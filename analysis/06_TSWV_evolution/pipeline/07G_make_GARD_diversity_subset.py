#!/usr/bin/env python3
import sys,csv
from pathlib import Path
from Bio import SeqIO
import numpy as np
inp=Path(sys.argv[1]); out=Path(sys.argv[2]); audit=Path(sys.argv[3]); maxn=int(sys.argv[4])
rec=list(SeqIO.parse(str(inp),'fasta'))
if not rec: raise SystemExit('empty GARD input alignment')
if len(rec)<=maxn:
    SeqIO.write(rec,str(out),'fasta');chosen=list(range(len(rec)))
else:
    seq=[str(r.seq).upper() for r in rec]; L=len(seq[0])
    if any(len(x)!=L for x in seq): raise SystemExit('unaligned input')
    trans=np.full(256,-1,dtype=np.int8)
    for i,c in enumerate(b'ACGT'): trans[c]=i
    raw=np.frombuffer(''.join(seq).encode('ascii'),dtype=np.uint8).reshape(len(seq),L)
    arr=trans[raw]
    valid=arr>=0
    ids=[r.id for r in rec]
    first=min(range(len(rec)),key=lambda i:ids[i]);chosen=[first]
    def all_dist_to(j):
        both=valid & valid[j]
        comp=both.sum(axis=1)
        diff=((arr!=arr[j]) & both).sum(axis=1)
        return np.divide(diff,comp,out=np.zeros(len(rec),dtype=float),where=comp>0)
    mind=all_dist_to(first);mind[first]=-1.0; chosen_set={first}
    while len(chosen)<maxn:
        candidates=np.array([i for i in range(len(rec)) if i not in chosen_set],dtype=int)
        mx=np.max(mind[candidates]); tied=candidates[np.isclose(mind[candidates],mx,rtol=0,atol=1e-15)]
        best=min(tied,key=lambda i:ids[i])
        chosen.append(int(best));chosen_set.add(int(best));d=all_dist_to(best);mind=np.minimum(mind,d)
        for i in chosen_set: mind[i]=-1.0
    SeqIO.write([rec[i] for i in chosen],str(out),'fasta')
sel=set(chosen)
with audit.open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sequence_id','selected_for_GARD','input_n','target_n'])
    for i,r in enumerate(rec):w.writerow([r.id,str(i in sel).upper(),len(rec),min(len(rec),maxn)])
print('GARD subset',len(rec),'->',len(chosen))
