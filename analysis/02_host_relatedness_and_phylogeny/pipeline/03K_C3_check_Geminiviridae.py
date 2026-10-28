#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse, csv, statistics, math

def read_fasta(path):
    rows=[]; h=None; b=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if h is not None: rows.append((h,"".join(b)))
                h=line[1:].split()[0]; b=[]
            else: b.append(line)
    if h is not None: rows.append((h,"".join(b)))
    return rows

def q(v,p):
    v=sorted(v)
    if not v:return None
    x=(len(v)-1)*p
    lo=math.floor(x);hi=math.ceil(x)
    if lo==hi:return v[lo]
    return v[lo]*(hi-x)+v[hi]*(x-lo)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    a=ap.parse_args()
    root=Path(a.project_root).resolve()
    base=root/"three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    out=base/"pipeline_results/03K_family_specific_markers/07_Geminiviridae_Rep_conserved"
    d=out/"AA95_alignment_QC"

    for p in [
        out/"Geminiviridae.Rep_conserved.summary.tsv",
        out/"Geminiviridae.Rep_conserved.03K_C3_summary.tsv",
        d/"Geminiviridae.Rep_conserved.AA95.mafft.trim.faa",
        d/"Geminiviridae.Rep_conserved.AA95.mafft.trim.filtered.faa",
        d/"Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv",
    ]:
        if not p.exists():
            raise SystemExit(f"MISSING {p}")

    print("===== EXTRACTION SUMMARY =====")
    print((out/"Geminiviridae.Rep_conserved.summary.tsv").read_text())

    print("===== PIPELINE SUMMARY =====")
    print((out/"Geminiviridae.Rep_conserved.03K_C3_summary.tsv").read_text())

    print("===== TRIMMED ALIGNMENT DISTRIBUTION =====")
    rows=read_fasta(d/"Geminiviridae.Rep_conserved.AA95.mafft.trim.faa")
    alen=[len(s) for _,s in rows]
    eff=[sum(c not in "-?." for c in s) for _,s in rows]
    gap=[s.count("-")/len(s) if s else 1 for _,s in rows]
    print("seq_n",len(rows))
    print("alignment_length", min(alen), statistics.median(alen), max(alen))
    print("effective_aa", min(eff), f"{q(eff,.05):.1f}", statistics.median(eff), f"{q(eff,.95):.1f}", max(eff))
    print("gap_fraction", f"{min(gap):.3f}", f"{q(gap,.05):.3f}", f"{statistics.median(gap):.3f}", f"{q(gap,.95):.3f}", f"{max(gap):.3f}")

    qc=list(csv.DictReader(open(d/"Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv"),delimiter="\t"))
    keep=sum(r["action"]=="keep" for r in qc)
    drop=len(qc)-keep
    print("QC keep",keep,"drop",drop)
    print("03K_C3_CHECK=PASS")

if __name__=="__main__":
    main()
