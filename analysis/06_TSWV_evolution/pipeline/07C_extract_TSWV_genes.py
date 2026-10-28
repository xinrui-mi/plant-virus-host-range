#!/usr/bin/env python3
import os,re,csv,statistics
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
ROOT=Path(os.environ['ANALYSIS_ROOT']);base=ROOT/'06_TSWV_evolution';gb=base/'work'/'TSWV_formal_current.gb'
out=base/'pipeline_results'/'01_gene_extraction'; genes=out/'genes'; meta=out/'metadata'; genes.mkdir(parents=True,exist_ok=True);meta.mkdir(parents=True,exist_ok=True)
if not gb.exists() or gb.stat().st_size==0: raise SystemExit('Fresh TSWV GenBank subset is missing; resolve membership gate first')
GENES=['RdRp','GnGc','NSm','N','NSs']
def norm(s): return re.sub(r'[^a-z0-9]+',' ',(s or '').lower()).strip()
def classify(feat):
    q=feat.qualifiers
    text=' '.join(sum([q.get(k,[]) for k in ['gene','product','note','function','label']],[]))
    t=norm(text)
    # High-specificity rules first; N is intentionally last to avoid "nonstructural" false matches.
    scores={g:0 for g in GENES}
    if re.search(r'\brdrp\b|rna dependent rna polymerase|large polymerase|polymerase l|\bl protein\b',t): scores['RdRp']+=10
    if re.search(r'gngc|gn gc|glycoprotein precursor|glycoprotein polyprotein|envelope glycoprotein',t): scores['GnGc']+=10
    if re.search(r'\bnsm\b|nonstructural protein m|movement protein',t): scores['NSm']+=10
    if re.search(r'\bnss\b|nonstructural protein s|silencing suppressor',t): scores['NSs']+=10
    if re.search(r'\bnucleocapsid\b|\bnucleoprotein\b|\bn protein\b',t) and 'nonstructural' not in t: scores['N']+=10
    # gene qualifier exact-ish clues
    gv=' '.join(q.get('gene',[]));g=norm(gv)
    if g in {'l','rdrp'}: scores['RdRp']+=8
    if g in {'g','gngc','gn gc','glycoprotein'}: scores['GnGc']+=8
    if g=='nsm': scores['NSm']+=8
    if g=='nss': scores['NSs']+=8
    if g=='n': scores['N']+=8
    best=max(scores.values())
    winners=[k for k,v in scores.items() if v==best and v>0]
    return (winners[0] if len(winners)==1 else ''),text,scores
records=[]; candidates=defaultdict(list); ambiguous=[]
for rec in SeqIO.parse(str(gb),'genbank'):
    for idx,feat in enumerate(rec.features,1):
        if feat.type!='CDS': continue
        gene,label,scores=classify(feat)
        nt=str(feat.extract(rec.seq)).upper().replace('U','T')
        aa=(feat.qualifiers.get('translation',[''])[0] or '').replace(' ','').replace('\n','').rstrip('*')
        if not aa and nt:
            try: aa=str(feat.extract(rec.seq).translate(to_stop=False)).rstrip('*')
            except Exception: aa=''
        row={'sequence_id':rec.id,'feature_index':idx,'assigned_gene':gene,'annotation_text':label,'nt_length':len(nt),'aa_length':len(aa),'nt':nt,'aa':aa,'score_vector':';'.join(f'{k}:{scores[k]}' for k in GENES)}
        records.append(row)
        if gene:candidates[(rec.id,gene)].append(row)
        elif max(scores.values())>0: ambiguous.append(row)
# exactly one CDS per accession x gene; if duplicate, prefer score then length and record audit
selected=[]; duplicate=[]
for key,z in candidates.items():
    gene=key[1]
    z=sorted(z,key=lambda r:(-int(dict(x.split(':') for x in r['score_vector'].split(';'))[gene]),-r['aa_length'],r['feature_index']))
    selected.append(z[0])
    if len(z)>1: duplicate.extend(z[1:])
# strict QC by basic CDS validity and robust relative-to-median length within each gene
med={g:statistics.median([r['nt_length'] for r in selected if r['assigned_gene']==g and r['nt_length']>0]) for g in GENES if any(r['assigned_gene']==g and r['nt_length']>0 for r in selected)}
strict=[]; audit=[]
lo=float(os.environ.get('TSWV_CDS_MIN_RELATIVE_LENGTH','0.90'));hi=float(os.environ.get('TSWV_CDS_MAX_RELATIVE_LENGTH','1.10'))
for r in selected:
    g=r['assigned_gene']; reasons=[]; nt=r['nt']; aa=r['aa']; m=med.get(g,0)
    if len(nt)%3!=0: reasons.append('CDS_length_not_multiple_of_3')
    if not aa: reasons.append('missing_translation')
    if '*' in aa: reasons.append('internal_stop')
    if m and not (lo*m <= len(nt) <= hi*m): reasons.append('outside_relative_length_window')
    # preserve ambiguous bases but flag excessive ambiguity
    amb=sum(c not in 'ACGT' for c in nt)
    if len(nt) and amb/len(nt)>0.02: reasons.append('ambiguous_nt_fraction_gt_0.02')
    ok=not reasons
    audit.append({**{k:r[k] for k in ['sequence_id','assigned_gene','feature_index','annotation_text','nt_length','aa_length']},'median_gene_nt_length':m,'strict_pass':str(ok).upper(),'failure_reasons':';'.join(reasons)})
    if ok: strict.append(r)
for g in GENES:
    zs=[r for r in strict if r['assigned_gene']==g]
    with (genes/f'{g}.strict.fna').open('w') as fn,(genes/f'{g}.strict.faa').open('w') as fa:
        for r in sorted(zs,key=lambda x:x['sequence_id']):fn.write(f">{r['sequence_id']}\n{r['nt']}\n");fa.write(f">{r['sequence_id']}\n{r['aa']}\n")
with (meta/'strict_qc.tsv').open('w',newline='') as f:
    fields=['sequence_id','assigned_gene','feature_index','annotation_text','nt_length','aa_length','median_gene_nt_length','strict_pass','failure_reasons'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(audit)
with (meta/'gene_summary.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['gene','selected_CDS_n','strict_pass_n','median_nt_length'])
    for g in GENES:w.writerow([g,sum(r['assigned_gene']==g for r in selected),sum(r['assigned_gene']==g for r in strict),med.get(g,'')])
with (meta/'annotation_ambiguity.tsv').open('w',newline='') as f:
    fields=['sequence_id','feature_index','annotation_text','score_vector','nt_length','aa_length'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows([{k:r[k] for k in fields} for r in ambiguous+duplicate])
print('TSWV gene extraction PASS', {g:sum(r['assigned_gene']==g for r in strict) for g in GENES})
