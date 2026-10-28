#!/usr/bin/env python3
from __future__ import annotations
import csv, os
from pathlib import Path

ROOT=Path(os.environ["ANALYSIS_ROOT"])
WORK=ROOT/"02_host_relatedness_and_phylogeny"/"work"/"marker_inputs"
EVALUE=float(os.environ.get("FIG3_RDRP_EVALUE","1e-5"))
RULES={"Potyviridae":(220,120),"Tombusviridae":(80,180)}


def fasta(path):
    d={};h=None;s=[]
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if h is not None:d[h]="".join(s)
                h=line[1:].split()[0];s=[]
            else:s.append(line.strip())
        if h is not None:d[h]="".join(s)
    return d


def loadmap(path):
    d={}
    with path.open() as f:
        for r in csv.DictReader(f,delimiter="\t"):d[r["orf_id"]]=r
    return d


def dom(path):
    out=[]
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():continue
            p=line.split()
            if len(p)<22:continue
            try:
                h=dict(target=p[0],tlen=int(p[2]),query=p[3],qlen=int(p[5]),full_e=float(p[6]),full_score=float(p[7]),ie=float(p[12]),dom_score=float(p[13]),env_from=int(p[19]),env_to=int(p[20]))
            except ValueError:continue
            if h["ie"]<=EVALUE:out.append(h)
    return out

qc=[]
for fam,(up,down) in RULES.items():
    seqs=fasta(WORK/f"{fam}.all_ORFs.faa"); mp=loadmap(WORK/f"{fam}.all_ORFs_map.tsv"); hits=dom(WORK/f"{fam}.RdRP.domtbl")
    best={}
    for h in hits:
        if h["target"] not in mp or h["target"] not in seqs:continue
        sid=mp[h["target"]]["sequence_id"]
        key=(h["dom_score"],-h["ie"])
        if sid not in best or key>(best[sid][0],-best[sid][1]["ie"]):best[sid]=(h["dom_score"],h)
    rows=[]
    for sid in sorted(best):
        _,h=best[sid]; oid=h["target"]; s=seqs[oid]
        start=max(1,h["env_from"]-up);end=min(len(s),h["env_to"]+down); marker=s[start-1:end]
        if len(marker)<120:continue
        rows.append(dict(sequence_id=sid,orf_id=oid,hmm_model=h["query"],i_evalue=h["ie"],domain_score=h["dom_score"],hit_from=h["env_from"],hit_to=h["env_to"],up=up,down=down,marker_start=start,marker_end=end,marker_length=len(marker),marker_seq=marker))
    with (WORK/f"{fam}.RdRP_expanded.raw.faa").open("w") as o:
        for r in rows:o.write(f">{r['sequence_id']}\n{r['marker_seq']}\n")
    with (WORK/f"{fam}.RdRP_marker_detail.tsv").open("w",newline="") as f:
        fields=["sequence_id","orf_id","hmm_model","i_evalue","domain_score","hit_from","hit_to","up","down","marker_start","marker_end","marker_length"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();
        for r in rows:w.writerow({k:r[k] for k in fields})
    qc.append([fam,len(seqs),len(hits),len(best),len(rows),up,down])
    print(fam,"ORFs",len(seqs),"qualifying_hits",len(hits),"sequence_best_hits",len(best),"markers",len(rows),"UP",up,"DOWN",down)
with (WORK/"03E_RNA_marker_selection_QC.tsv").open("w",newline="") as f:
    w=csv.writer(f,delimiter="\t");w.writerow(["family","orf_n","qualifying_domain_hit_n","sequence_with_best_hit_n","marker_sequence_n","up_aa","down_aa"]);w.writerows(qc)
print("RNA marker selection PASS")
