#!/usr/bin/env python3
import os,csv
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
ROOT=Path(os.environ['ANALYSIS_ROOT']);base=ROOT/'06_TSWV_evolution';indir=base/'pipeline_results'/'01_gene_extraction'/'genes';out=base/'pipeline_results'/'02_exact_haplotypes';out.mkdir(parents=True,exist_ok=True)
for gene in ['RdRp','GnGc','NSm','N','NSs']:
    fna=indir/f'{gene}.strict.fna';faa=indir/f'{gene}.strict.faa'
    nt={r.id:str(r.seq).upper() for r in SeqIO.parse(str(fna),'fasta')}; aa={r.id:str(r.seq) for r in SeqIO.parse(str(faa),'fasta')}
    groups=defaultdict(list)
    for sid,s in nt.items(): groups[s].append(sid)
    reps=[]; rows=[]
    for i,(s,members) in enumerate(sorted(groups.items(),key=lambda z:sorted(z[1])[0]),1):
        rep=sorted(members)[0]; hid=f'{gene}_H{i:04d}';reps.append((rep,s,aa[rep]));
        for m in sorted(members): rows.append([gene,hid,rep,m,str(m==rep).upper(),len(members)])
    with (out/f'{gene}.haplotypes.fna').open('w') as fn,(out/f'{gene}.haplotypes.faa').open('w') as fa:
        for rep,ntseq,aaseq in reps:fn.write(f'>{rep}\n{ntseq}\n');fa.write(f'>{rep}\n{aaseq}\n')
    with (out/f'{gene}.haplotype_map.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['gene','haplotype_id','representative_sequence_id','member_sequence_id','is_representative','haplotype_member_n']);w.writerows(rows)
    print(gene,'strict',len(nt),'exact_haplotypes',len(reps))
