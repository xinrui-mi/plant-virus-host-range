#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path

def read_fasta(path):
    out=[];h=None;s=[]
    for line in path.open():
        if line.startswith(">"):
            if h is not None:out.append((h,"".join(s)))
            h=line[1:].strip();s=[]
        else:s.append(line.strip())
    if h is not None:out.append((h,"".join(s)))
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument("input");ap.add_argument("output");ap.add_argument("--max-gap",type=float,default=0.5);ap.add_argument("--min-effective-aa",type=int,default=100);ap.add_argument("--report",required=True);a=ap.parse_args()
    rows=read_fasta(Path(a.input));kept=[];report=[]
    for h,s in rows:
        if not s:report.append((h,0,0,1,"empty"));continue
        gap=s.count("-");eff=sum(ch not in "-?." for ch in s);gf=gap/len(s)
        keep=(gf<=a.max_gap and eff>=a.min_effective_aa)
        report.append((h,len(s),eff,gf,"keep" if keep else "drop"))
        if keep:kept.append((h,s))
    with Path(a.output).open("w") as o:
        for h,s in kept:o.write(f">{h}\n{s}\n")
    with Path(a.report).open("w") as o:
        o.write("tip\talignment_length\teffective_aa\tgap_fraction\taction\n")
        for r in report:o.write("\t".join(map(str,r))+"\n")
    print("Alignment filter input",len(rows),"kept",len(kept),"dropped",len(rows)-len(kept))
if __name__=="__main__":main()
