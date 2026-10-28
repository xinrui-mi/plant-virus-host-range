#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, re

FAMS=[
"Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
"Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
"Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
]

def fasta_n(p):
    return sum(1 for x in open(p) if x.startswith(">"))

def tree_tip_n(p):
    # Alignment count is authoritative for expected tips.
    # Treefile nonempty + report success are checked separately.
    return None

def parse_report(path):
    txt=path.read_text(errors="replace")
    model=""
    logl=""
    m=re.search(r"Best-fit model according to BIC:\s+(\S+)",txt)
    if m:model=m.group(1)
    if not model:
        m=re.search(r"Best-fit model:\s+(\S+)",txt)
        if m:model=m.group(1)
    m=re.search(r"Log-likelihood of the tree:\s+([-\d.]+)",txt)
    if m:logl=m.group(1)
    return model,logl

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    a=ap.parse_args()
    root=Path(a.project_root).resolve()
    base=root/"three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    droot=base/"pipeline_results/03K_family_specific_markers/08_final_12family_trees"
    out=droot/"trees"

    rows=[]
    for fam in FAMS:
        aln=droot/"inputs"/f"{fam}.final_alignment.faa"
        D=out/fam
        tree=D/f"{fam}.treefile"
        rep=D/f"{fam}.iqtree"
        status="PASS" if tree.exists() and tree.stat().st_size>0 and rep.exists() and rep.stat().st_size>0 else "FAIL"
        model="";logl=""
        if rep.exists():
            model,logl=parse_report(rep)
        rows.append({
            "virus_family":fam,
            "alignment_tip_n":fasta_n(aln) if aln.exists() else "NA",
            "treefile_present":"YES" if tree.exists() and tree.stat().st_size>0 else "NO",
            "iqtree_report_present":"YES" if rep.exists() and rep.stat().st_size>0 else "NO",
            "best_fit_model":model,
            "log_likelihood":logl,
            "status":status,
        })

    p=droot/"03K_D_tree_QC.tsv"
    fields=list(rows[0].keys())
    with open(p,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader();w.writerows(rows)

    print(f"tree_QC={p}")
    print(f"PASS={sum(r['status']=='PASS' for r in rows)}/12")

if __name__=="__main__":
    main()
