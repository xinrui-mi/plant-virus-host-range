#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, statistics, math
from collections import Counter

def read_fasta(path):
    rows = []
    h = None
    buf = []
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if h is not None:
                    rows.append((h, "".join(buf)))
                h=line[1:].split()[0]
                buf=[]
            else:
                buf.append(line)
    if h is not None:
        rows.append((h, "".join(buf)))
    return rows

def q(vals, p):
    vals=sorted(vals)
    if not vals: return None
    if len(vals)==1: return vals[0]
    x=(len(vals)-1)*p
    lo=math.floor(x); hi=math.ceil(x)
    if lo==hi: return vals[lo]
    return vals[lo]*(hi-x)+vals[hi]*(x-lo)

def stats(vals):
    if not vals:
        return "NA"
    return (
        f"n={len(vals)} min={min(vals):.3f} q05={q(vals,.05):.3f} "
        f"median={statistics.median(vals):.3f} q95={q(vals,.95):.3f} max={max(vals):.3f}"
    )

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    a=ap.parse_args()

    root=Path(a.project_root).resolve()
    base=root/"three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    d=base/"pipeline_results/03K_family_specific_markers/06_AA95_alignment_QC/Geminiviridae"

    raw=base/"pipeline_results/03K_family_specific_markers/03_final_marker_recovery/Geminiviridae.marker.raw.faa"
    aa95=d/"Geminiviridae.AA95.faa"
    aln=d/"Geminiviridae.AA95.mafft.faa"
    trim=d/"Geminiviridae.AA95.mafft.trim.faa"
    qc=d/"Geminiviridae.alignment_sequence_QC.tsv"

    required=[raw,aa95,aln,trim,qc]
    for p in required:
        if not p.exists():
            raise SystemExit(f"Missing: {p}")

    print("===== FASTA COUNTS / LENGTHS =====")
    for label,p in [("raw_marker",raw),("AA95",aa95),("MAFFT",aln),("trimAl",trim)]:
        rows=read_fasta(p)
        lens=[len(s) for _,s in rows]
        eff=[sum(c not in "-?." for c in s) for _,s in rows]
        gaps=[s.count("-")/len(s) if s else 1 for _,s in rows]
        print(f"{label}: seq_n={len(rows)}")
        print(f"  total/alignment length: {stats(lens)}")
        print(f"  effective aa: {stats(eff)}")
        print(f"  gap fraction: {stats(gaps)}")

    print()
    print("===== FROZEN QC REPORT =====")
    rows=[]
    with open(qc) as f:
        rd=csv.DictReader(f, delimiter="\t")
        for r in rd:
            rows.append(r)

    print("qc_rows", len(rows))
    actions=Counter(r["action"] for r in rows)
    print("action_counts", dict(actions))

    eff=[int(float(r["effective_aa"])) for r in rows]
    gap=[float(r["gap_fraction"]) for r in rows]
    alen=[int(float(r["alignment_length"])) for r in rows]
    print("alignment_length:", stats(alen))
    print("effective_aa:", stats(eff))
    print("gap_fraction:", stats(gap))

    reason=Counter()
    for r in rows:
        e=int(float(r["effective_aa"]))
        g=float(r["gap_fraction"])
        if g>0.50 and e<100:
            reason["gap_and_effective"]+=1
        elif g>0.50:
            reason["gap_only"]+=1
        elif e<100:
            reason["effective_aa_only"]+=1
        else:
            reason["passes_thresholds"]+=1
    print("derived_drop_reason_counts", dict(reason))

    print()
    print("===== FIRST 20 QC ROWS =====")
    for r in rows[:20]:
        print(
            r["tip"],
            "alignment_length="+r["alignment_length"],
            "effective_aa="+r["effective_aa"],
            "gap_fraction="+r["gap_fraction"],
            "action="+r["action"]
        )

    print()
    print("===== INTERPRETATION FLAG =====")
    n=len(rows)
    effective_fail=reason["effective_aa_only"]+reason["gap_and_effective"]
    gap_fail=reason["gap_only"]+reason["gap_and_effective"]
    if n and effective_fail==n and gap_fail<n:
        print("PRIMARY_FAILURE=effective_length")
    elif n and gap_fail==n and effective_fail<n:
        print("PRIMARY_FAILURE=gap_fraction")
    elif n and gap_fail==n and effective_fail==n:
        print("PRIMARY_FAILURE=both_gap_and_effective_for_all")
    elif reason["passes_thresholds"]==0:
        print("PRIMARY_FAILURE=mixed_but_all_fail")
    else:
        print("PRIMARY_FAILURE=unexpected_some_should_pass")

    print("03K_C2_GEMINI_AUDIT=PASS")

if __name__=="__main__":
    main()
