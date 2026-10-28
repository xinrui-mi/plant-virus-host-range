#!/usr/bin/env python3
from __future__ import annotations
import csv,re,sys
clstr,fam,out=sys.argv[1:4]
rows=[];cid=None;members=[]

def flush():
    global rows,members,cid
    if not members:return
    reps=[x for x in members if x[1]]
    if len(reps)!=1:raise SystemExit(f"{fam} {cid}: expected exactly one CD-HIT representative, got {len(reps)}")
    rep=reps[0][0]
    for seq,isrep,pct in members:rows.append([fam,cid,rep,seq,int(isrep),pct])

for line in open(clstr):
    if line.startswith(">Cluster"):
        flush();cid=line.strip()[1:].replace(" ","_");members=[];continue
    m=re.search(r">(.+?)\.\.\.",line)
    if not m:continue
    seq=m.group(1);isrep=line.rstrip().endswith("*");mi=re.search(r"at\s+([0-9.]+)%",line);pct="100" if isrep else (mi.group(1) if mi else "")
    members.append((seq,isrep,pct))
flush()
with open(out,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t");w.writerow(["virus_family","cdhit95_cluster_id","representative_sequence_id","member_sequence_id","is_representative","identity_to_representative_percent"]);w.writerows(rows)
print(f"CD-HIT membership PASS {fam}: clusters={len(set(r[1] for r in rows))} members={len(rows)}")
