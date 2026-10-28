#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Study-aware validation of candidate TSWV L/M/S isolate pairings.

Input (from 07K):
  tables/09_TSWV_segment_pairing_audit.tsv
  work/TSWV_formal_current.gb

This script does NOT infer reassortment. It only establishes a defensible set
of isolate-matched L/M/S records for subsequent segment phylogenies.

Strict-pass logic:
  - starts from LMS_unique_high_confidence candidates from 07K;
  - one L, one M, one S accession, complete coding evidence, no country/date conflict;
  - no explicit GenBank segment-qualifier contradiction;
  - requires independent contextual linkage across the three records:
      * shared primary reference, OR
      * shared direct-submission signature, OR
      * coherent accession-series + host/country context, OR
      * for non-generic isolate names, coherent host/country/date context.
    Generic isolate names (e.g. "1", "G", "100/23") require one of the first
    three stronger linkages and are not accepted from the isolate string alone.
"""
from __future__ import annotations
import csv,re,sys
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def find_root():
    p=Path(__file__).resolve()
    for x in [p.parents[2],Path.cwd().resolve(),*Path.cwd().resolve().parents]:
        if (x/'06_TSWV_evolution').is_dir(): return x
    raise SystemExit('Cannot locate analysis root')

ROOT=find_root(); BASE=ROOT/'06_TSWV_evolution'
AUD=BASE/'tables'/'09_TSWV_segment_pairing_audit.tsv'
GB=BASE/'work'/'TSWV_formal_current.gb'
if not AUD.exists(): raise SystemExit(f'Missing 07K pairing audit: {AUD}')
if not GB.exists(): raise SystemExit(f'Missing TSWV GenBank: {GB}')

def txt(x): return re.sub(r'\s+',' ',str(x or '')).strip()
def norm(x): return re.sub(r'[^a-z0-9]+',' ',txt(x).lower()).strip()
def host_core(x):
    z=norm(x)
    toks=z.split()
    stop={'cultivar','cv','hybrid','var','variety','strain'}
    out=[]
    for t in toks:
        if t in stop: break
        out.append(t)
        if len(out)>=2: break
    return ' '.join(out)
def country_core(x): return norm(txt(x).split(':',1)[0])
def is_generic_isolate(x):
    s=txt(x); compact=re.sub(r'[^A-Za-z0-9]','',s)
    if not compact: return True
    if re.fullmatch(r'\d+(?:[._/-]?\d+)*',s): return True
    if len(compact)<=2: return True
    return False

def first_source(rec):
    for f in rec.features:
        if f.type=='source': return f.qualifiers
    return {}
def first(q,keys):
    for k in keys:
        v=q.get(k,[])
        if v: return txt(v[0])
    return ''
def ref_signature(rec,direct=False):
    refs=rec.annotations.get('references',[]) or []
    pool=[]
    for r in refs:
        title=txt(getattr(r,'title','')); authors=txt(getattr(r,'authors','')); journal=txt(getattr(r,'journal','')); pmid=txt(getattr(r,'pubmed_id',''))
        isdir=title.lower().startswith('direct submission') or journal.lower().startswith('submitted')
        if isdir==direct:
            # PMIDs are strongest; otherwise title+authors for primary refs and authors+journal for direct submissions.
            if pmid: return 'pmid:'+pmid
            if direct: pool.append(norm(authors)+'|'+norm(journal))
            else: pool.append(norm(title)+'|'+norm(authors))
    return pool[0] if pool else ''
def acc_series(accessions):
    vals=[]; prefixes=[]
    for a in accessions:
        a=a.split('.')[0]
        m=re.fullmatch(r'([A-Za-z]+)(\d+)',a)
        if not m: return False,'',None
        prefixes.append(m.group(1));vals.append(int(m.group(2)))
    if len(set(prefixes))!=1: return False,prefixes[0] if prefixes else '',None
    span=max(vals)-min(vals)
    return span<=30,prefixes[0],span

def write_tsv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)

records={r.id:r for r in SeqIO.parse(str(GB),'genbank')}
cands=[r for r in csv.DictReader(AUD.open(encoding='utf-8'),delimiter='\t') if r['pairing_status']=='LMS_unique_high_confidence']
out=[]
for r in cands:
    acc={'L':r['L_accessions'],'M':r['M_accessions'],'S':r['S_accessions']}
    if any(';' in acc[s] or not acc[s] for s in 'LMS'):
        continue
    recs={s:records.get(acc[s]) for s in 'LMS'}
    if any(v is None for v in recs.values()):
        status='REVIEW_missing_GenBank_record'; evidence=[]
    else:
        src={s:first_source(recs[s]) for s in 'LMS'}
        segq={s:first(src[s],['segment']).upper() for s in 'LMS'}
        def norm_seg(v):
            q=re.sub(r'[^A-Z0-9]+',' ',v.upper()).strip()
            mp={'L':'L','RNA L':'L','L RNA':'L','RNA1':'L','RNA 1':'L','1':'L',
                'M':'M','RNA M':'M','M RNA':'M','RNA2':'M','RNA 2':'M','2':'M',
                'S':'S','RNA S':'S','S RNA':'S','RNA3':'S','RNA 3':'S','3':'S'}
            return mp.get(q,q)
        segment_mismatch=[s for s in 'LMS' if segq[s] and norm_seg(segq[s])!=s]
        hosts=[host_core(first(src[s],['host','lab_host'])) for s in 'LMS']
        countries=[country_core(first(src[s],['geo_loc_name','country'])) for s in 'LMS']
        dates=[norm(first(src[s],['collection_date'])) for s in 'LMS']
        host_compatible=len({x for x in hosts if x})<=1
        country_compatible=len({x for x in countries if x})<=1
        date_compatible=len({x for x in dates if x})<=1
        prim=[ref_signature(recs[s],False) for s in 'LMS']
        direct=[ref_signature(recs[s],True) for s in 'LMS']
        shared_primary=bool(prim[0] and len(set(prim))==1)
        shared_direct=bool(direct[0] and len(set(direct))==1)
        series_ok,series_prefix,series_span=acc_series([acc[s] for s in 'LMS'])
        generic=is_generic_isolate(r['isolate_raw_names'])
        strong=[]
        if shared_primary: strong.append('shared_primary_reference')
        if shared_direct: strong.append('shared_direct_submission')
        if series_ok and host_compatible and country_compatible: strong.append('coherent_accession_series_context')
        moderate=(not generic and host_compatible and country_compatible and date_compatible)
        if moderate: strong.append('non_generic_isolate_context')
        if segment_mismatch:
            status='REVIEW_segment_qualifier_mismatch'
        elif generic and not any(x in strong for x in ['shared_primary_reference','shared_direct_submission','coherent_accession_series_context']):
            status='REVIEW_generic_isolate_without_strong_context'
        elif strong:
            status='STRICT_PASS'
        else:
            status='REVIEW_insufficient_independent_context'
        evidence=strong

    row={**r,
        'L_segment_qualifier':segq.get('L','') if all(v is not None for v in recs.values()) else '',
        'M_segment_qualifier':segq.get('M','') if all(v is not None for v in recs.values()) else '',
        'S_segment_qualifier':segq.get('S','') if all(v is not None for v in recs.values()) else '',
        'segment_qualifier_mismatch_n':len(segment_mismatch) if all(v is not None for v in recs.values()) else '',
        'host_context_compatible':str(host_compatible).upper() if all(v is not None for v in recs.values()) else '',
        'country_context_compatible':str(country_compatible).upper() if all(v is not None for v in recs.values()) else '',
        'date_context_compatible':str(date_compatible).upper() if all(v is not None for v in recs.values()) else '',
        'shared_primary_reference':str(shared_primary).upper() if all(v is not None for v in recs.values()) else '',
        'shared_direct_submission':str(shared_direct).upper() if all(v is not None for v in recs.values()) else '',
        'coherent_accession_series':str(series_ok).upper() if all(v is not None for v in recs.values()) else '',
        'accession_series_prefix':series_prefix if all(v is not None for v in recs.values()) else '',
        'accession_series_span':series_span if all(v is not None for v in recs.values()) else '',
        'generic_isolate_name':str(generic).upper() if all(v is not None for v in recs.values()) else '',
        'independent_linkage_evidence':';'.join(evidence),
        'study_aware_status':status,
    }
    out.append(row)

out_path=BASE/'tables'/'10_TSWV_segment_pairing_study_aware_audit.tsv'
fields=list(out[0].keys()) if out else []
write_tsv(out_path,out,fields)
strict=[r for r in out if r['study_aware_status']=='STRICT_PASS']
strict=sorted(strict,key=lambda r:(r['isolate_key'],r['L_accessions']))
strict_rows=[]
for i,r in enumerate(strict,1):
    strict_rows.append({
        'paired_id':f'TSWVISO_{i:04d}','isolate_key':r['isolate_key'],'isolate_name':r['isolate_raw_names'],
        'L_accession':r['L_accessions'],'M_accession':r['M_accessions'],'S_accession':r['S_accessions'],
        'country':r['country_set'],'collection_date':r['collection_date_set'],'host':r['host_set'],
        'generic_isolate_name':r['generic_isolate_name'],'linkage_evidence':r['independent_linkage_evidence']
    })
strict_path=BASE/'tables'/'11_TSWV_LMS_strict_isolates.tsv'
write_tsv(strict_path,strict_rows,list(strict_rows[0].keys()) if strict_rows else ['paired_id'])

from collections import Counter
cnt=Counter(r['study_aware_status'] for r in out)
qc=[
    {'metric':'07K_LMS_unique_high_confidence_candidate_n','value':len(cands)},
    {'metric':'study_aware_STRICT_PASS_n','value':len(strict)},
    {'metric':'study_aware_review_n','value':len(out)-len(strict)},
]
for k,v in sorted(cnt.items()): qc.append({'metric':'status::'+k,'value':v})
qcp=BASE/'qc'/'05_TSWV_segment_pairing_study_aware_QC.tsv';write_tsv(qcp,qc,['metric','value'])
print('TSWV study-aware segment pairing PASS')
for z in qc: print(f"{z['metric']}: {z['value']}")
print('Audit :',out_path)
print('Strict:',strict_path)
print('QC    :',qcp)
if len(strict)<50:
    print('WARNING: <50 strict paired isolates; reassortment analysis may be underpowered.',file=sys.stderr)
