#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare matched, orientation-normalized TSWV L/M/S whole-segment inputs.

Uses the study-aware strict isolate table from 07L. All three segment FASTAs
contain exactly the same paired_id labels. Whole records are oriented using a
coding anchor (L: RdRp; M: GnGc; S: N), so cRNA/genomic-RNA submissions are
placed in a consistent direction before alignment.
"""
from __future__ import annotations
import csv,os,re,statistics,sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT=Path(os.environ.get('ANALYSIS_ROOT','')).resolve()
if not (ROOT/'06_TSWV_evolution').is_dir():
    raise SystemExit('ANALYSIS_ROOT must be set to the analysis root (the directory containing 06_TSWV_evolution/); got: %s' % ROOT)
BASE=ROOT/'06_TSWV_evolution'
PAIR=BASE/'tables'/'11_TSWV_LMS_strict_isolates.tsv'
GB=BASE/'work'/'TSWV_formal_current.gb'
OUT=BASE/'pipeline_results'/'06_segment_concordance'/'inputs'
OUT.mkdir(parents=True,exist_ok=True)
if not PAIR.exists(): raise SystemExit(f'Missing study-aware strict pairs: {PAIR}')
if not GB.exists(): raise SystemExit(f'Missing GenBank: {GB}')

MIN_REL=float(os.environ.get('TSWV_SEGMENT_MIN_REL_LENGTH','0.90'))
MAX_REL=float(os.environ.get('TSWV_SEGMENT_MAX_REL_LENGTH','1.10'))
MAX_AMB=float(os.environ.get('TSWV_SEGMENT_MAX_AMBIGUOUS_FRACTION','0.02'))

GENES=['RdRp','GnGc','NSm','N','NSs']
def norm(s): return re.sub(r'[^a-z0-9]+',' ',(s or '').lower()).strip()
def classify(feat):
    q=feat.qualifiers
    text=' '.join(sum([q.get(k,[]) for k in ['gene','product','note','function','label']],[]))
    t=norm(text);scores={g:0 for g in GENES}
    if re.search(r'\brdrp\b|rna dependent rna polymerase|large polymerase|polymerase l|\bl protein\b',t): scores['RdRp']+=10
    if re.search(r'gngc|gn gc|glycoprotein precursor|glycoprotein polyprotein|envelope glycoprotein',t): scores['GnGc']+=10
    if re.search(r'\bnsm\b|nonstructural protein m|movement protein',t): scores['NSm']+=10
    if re.search(r'\bnss\b|nonstructural protein s|silencing suppressor',t): scores['NSs']+=10
    if re.search(r'\bnucleocapsid\b|\bnucleoprotein\b|\bn protein\b',t) and 'nonstructural' not in t: scores['N']+=10
    gv=norm(' '.join(q.get('gene',[])))
    if gv in {'l','rdrp'}: scores['RdRp']+=8
    if gv in {'g','gngc','gn gc','glycoprotein'}: scores['GnGc']+=8
    if gv=='nsm': scores['NSm']+=8
    if gv=='nss': scores['NSs']+=8
    if gv=='n': scores['N']+=8
    best=max(scores.values());win=[k for k,v in scores.items() if v==best and v>0]
    return win[0] if len(win)==1 else ''

def orient_record(rec,segment):
    anchor={'L':'RdRp','M':'GnGc','S':'N'}[segment]
    cand=[]
    for idx,f in enumerate(rec.features):
        if f.type!='CDS': continue
        g=classify(f)
        if g==anchor:
            strand=getattr(f.location,'strand',None)
            length=len(f.extract(rec.seq))
            cand.append((length,idx,strand))
    if not cand:
        return None,'anchor_not_found'
    cand.sort(reverse=True);strand=cand[0][2]
    if strand==1: return str(rec.seq).upper().replace('U','T'),'forward_anchor'
    if strand==-1: return str(rec.seq.reverse_complement()).upper().replace('U','T'),'reverse_complemented_anchor'
    return None,'anchor_strand_unresolved'

def write_tsv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)

pairs=list(csv.DictReader(PAIR.open(encoding='utf-8'),delimiter='\t'))
recmap={r.id:r for r in SeqIO.parse(str(GB),'genbank')}
raw=[]
for p in pairs:
    z={'paired_id':p['paired_id'],'isolate_name':p['isolate_name'],'country':p['country'],'collection_date':p['collection_date'],'host':p['host']}
    ok=True
    for seg in 'LMS':
        acc=p[f'{seg}_accession']; rec=recmap.get(acc)
        if rec is None:
            z[f'{seg}_status']='missing_record';ok=False;continue
        seq,orient=orient_record(rec,seg)
        z[f'{seg}_accession']=acc;z[f'{seg}_orientation']=orient
        if seq is None:
            z[f'{seg}_status']=orient;ok=False;continue
        amb=sum(c not in 'ACGT' for c in seq)/len(seq) if seq else 1
        z[f'{seg}_length']=len(seq);z[f'{seg}_ambiguous_fraction']=amb;z[f'{seg}_seq']=seq;z[f'{seg}_status']='candidate'
    z['initial_all_three_ok']=ok;raw.append(z)

# Segment-specific robust medians from orientation-resolved records.
med={}
for seg in 'LMS':
    vals=[r.get(f'{seg}_length',0) for r in raw if r.get(f'{seg}_status')=='candidate' and r.get(f'{seg}_length',0)>0]
    if not vals: raise SystemExit(f'No orientation-resolved {seg} records')
    med[seg]=statistics.median(vals)

kept=[];audit=[]
for r in raw:
    reasons=[]
    for seg in 'LMS':
        if r.get(f'{seg}_status')!='candidate': reasons.append(f'{seg}:{r.get(f"{seg}_status", "missing")}');continue
        L=r[f'{seg}_length']; amb=r[f'{seg}_ambiguous_fraction'];m=med[seg]
        if not (MIN_REL*m<=L<=MAX_REL*m): reasons.append(f'{seg}:outside_relative_length_window')
        if amb>MAX_AMB: reasons.append(f'{seg}:ambiguous_fraction_gt_{MAX_AMB}')
    keep=not reasons
    a={k:v for k,v in r.items() if not k.endswith('_seq')}
    a.update({'L_median_length':med['L'],'M_median_length':med['M'],'S_median_length':med['S'],'strict_segment_pass':str(keep).upper(),'failure_reasons':';'.join(reasons)})
    audit.append(a)
    if keep: kept.append(r)

if len(kept)<50: raise SystemExit(f'Only {len(kept)} matched isolates pass whole-segment QC; stop for review')
# FASTAs with EXACT same ordered paired_id set.
kept=sorted(kept,key=lambda r:r['paired_id'])
for seg in 'LMS':
    with (OUT/f'{seg}.matched_oriented.fna').open('w') as f:
        for r in kept: f.write(f'>{r["paired_id"]}\n{r[f"{seg}_seq"]}\n')

meta=[]
for r in kept:
    meta.append({
        'paired_id':r['paired_id'],'isolate_name':r['isolate_name'],'country':r['country'],'collection_date':r['collection_date'],'host':r['host'],
        'L_accession':r['L_accession'],'M_accession':r['M_accession'],'S_accession':r['S_accession'],
        'L_length':r['L_length'],'M_length':r['M_length'],'S_length':r['S_length'],
        'L_orientation':r['L_orientation'],'M_orientation':r['M_orientation'],'S_orientation':r['S_orientation']
    })
write_tsv(OUT/'LMS_matched_isolate_metadata.tsv',meta,list(meta[0].keys()))
write_tsv(OUT/'LMS_whole_segment_QC.tsv',audit,list(audit[0].keys()))
print('TSWV LMS segment input preparation PASS')
print('study_aware_pair_n',len(pairs),'whole_segment_QC_pass_n',len(kept),'median_lengths',med)
