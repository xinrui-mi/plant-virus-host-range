#!/usr/bin/env python3
import os,gzip,csv
from pathlib import Path
REPO=Path(os.environ.get("REPO_ROOT", Path(os.environ.get("ANALYSIS_ROOT",".")).parent))
ROOT=Path(os.environ["ANALYSIS_ROOT"])
seqfile=Path(os.environ.get("SEQ_FINAL", REPO/"data"/"analysis_input"/"sequence_analysis_input.tsv.gz"))
rdrp=Path(os.environ.get("RDRP_CORE_FASTA", REPO/"external"/"08_RdRP_core.faa"))
out=ROOT/"02_host_relatedness_and_phylogeny"/"work"/"marker_inputs"; out.mkdir(parents=True,exist_ok=True)
rows={}
_opener = gzip.open if str(seqfile).endswith(".gz") else open
with _opener(seqfile,"rt",encoding='utf-8-sig') as f:
    for r in csv.DictReader(f,delimiter='\t'): rows[r['sequence_id']]=r
def current_family_eligible(r):
    return (r.get('virus_ictv_match_rank','').strip().lower() in {'species','genus','family'} and r.get('virus_family','').strip())
families=['Potyviridae','Tombusviridae']; targets={x:set() for x in families}
for sid,r in rows.items():
    if current_family_eligible(r) and r.get('virus_family') in targets and r.get('rdrp_operational_cluster_id','').strip(): targets[r['virus_family']].add(sid)
for fam in ['Geminiviridae','Caulimoviridae']:
    p=out/f"{fam}.target_sequence_ids.txt"
    ids=sorted(sid for sid,r in rows.items() if current_family_eligible(r) and r.get('virus_family')==fam)
    p.write_text('\n'.join(ids)+('\n' if ids else ''))

def fasta_iter(path):
    h=None; seq=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if h is not None: yield h,''.join(seq)
                h=line[1:].strip();seq=[]
            else: seq.append(line.strip())
        if h is not None: yield h,''.join(seq)
writers={fam:(out/f"{fam}.RdRP_core.raw.faa").open('w') for fam in families}; found={fam:set() for fam in families}; unmatched=[]
allids=set(rows)
for h,s in fasta_iter(rdrp):
    tok=h.split()[0]; cand=[tok,tok.split('|')[0]]
    sid=next((x for x in cand if x in allids),None)
    if sid is None:
        unmatched.append(tok);continue
    fam=rows[sid].get('virus_family','')
    if fam in writers and sid in targets[fam]:
        writers[fam].write(f">{sid}\n{s}\n");found[fam].add(sid)
for x in writers.values(): x.close()
with (out/'03D_marker_input_QC.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['family','target_RdRP_positive_sequence_n','recovered_RdRP_core_n','missing_n'])
    for fam in families:w.writerow([fam,len(targets[fam]),len(found[fam]),len(targets[fam]-found[fam])])

for fam in families:
    missing=targets[fam]-found[fam]
    if missing:
        raise SystemExit(f'{fam}: {len(missing)} RdRP-positive target sequences missing from frozen 08_RdRP_core.faa')
print('03D marker input prep PASS')
