#!/usr/bin/env python3
from __future__ import annotations
import csv, json, os, re, sys
from pathlib import Path

GENES = ["RdRp", "GnGc", "NSm", "N", "NSs"]
VALID = {
    "FINAL_GARD",
    "PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1",
    "PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1",
}


def root_from_env_or_file() -> Path:
    if os.environ.get("ANALYSIS_ROOT"):
        return Path(os.environ["ANALYSIS_ROOT"]).resolve()
    return Path(__file__).resolve().parents[2]


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def read_fasta(path: Path):
    out=[]; name=None; chunks=[]
    with path.open() as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: out.append((name,''.join(chunks)))
                name=line[1:].split()[0]; chunks=[]
            else: chunks.append(line)
        if name is not None: out.append((name,''.join(chunks)))
    if not out: raise ValueError(f"No sequences: {path}")
    lens={len(s) for _,s in out}
    if len(lens)!=1: raise ValueError(f"Non-rectangular FASTA: {path}")
    return out, next(iter(lens))


def numeric_values(obj):
    if isinstance(obj,(int,float)) and not isinstance(obj,bool): yield float(obj)
    elif isinstance(obj,list):
        for x in obj: yield from numeric_values(x)
    elif isinstance(obj,dict):
        for k in sorted(obj,key=str): yield from numeric_values(obj[k])


def _partition_key(k):
    try:return (0,int(k))
    except:return (1,str(k))


def sorted_partition_data(js):
    """Parse a completed multi-partition GARD breakpointData object.

    This parser is intentionally used only when the final GARD model has >=1
    breakpoint. Zero-breakpoint GARD outputs are heterogeneous across HyPhy
    versions/checkpoint states (some encode a whole-alignment partition with a
    one-element `bps` object), so they are handled separately by
    single_tree_partition().
    """
    bd=js.get('breakpointData')
    if not isinstance(bd,dict) or not bd: raise ValueError('breakpointData missing/empty')
    parts=[]
    for k in sorted(bd,key=_partition_key):
        v=bd[k]
        if not isinstance(v,dict):
            raise ValueError(f'Invalid breakpointData partition {k}')
        nums=list(numeric_values(v.get('bps',[])))
        if len(nums)<2: raise ValueError(f'Cannot parse bps for partition {k}')
        tree=v.get('tree','')
        if not isinstance(tree,str) or '(' not in tree: raise ValueError(f'Missing tree for partition {k}')
        parts.append({'json_partition':str(k),'raw_start':int(round(nums[0])),'raw_end':int(round(nums[1])),'tree':tree.strip()})
    return parts


def single_tree_partition(js,total_len):
    """Return a whole-alignment partition and a GARD-matched tree.

    For a final no-breakpoint model we do NOT need breakpoint coordinates; we
    only need the GARD tree fitted to the exact diversity-subset taxa.  This is
    safer than using the original full-haplotype IQ-TREE tree, whose tip set can
    be larger than the GARD subset.
    """
    bd=js.get('breakpointData')
    if isinstance(bd,dict):
        for k in sorted(bd,key=_partition_key):
            v=bd[k]
            if isinstance(v,dict):
                tree=v.get('tree','')
                if isinstance(tree,str) and '(' in tree:
                    return {'json_partition':str(k),'raw_start':1,'raw_end':int(total_len),
                            'tree':tree.strip()}
    # Defensive fallback for HyPhy JSON variants where the fitted tree is stored
    # under `trees` rather than in breakpointData.
    trees=js.get('trees')
    candidates=[]
    if isinstance(trees,str): candidates=[trees]
    elif isinstance(trees,list): candidates=trees
    elif isinstance(trees,dict): candidates=[trees[k] for k in sorted(trees,key=_partition_key)]
    for obj in candidates:
        if isinstance(obj,str) and '(' in obj:
            return {'json_partition':'single_tree','raw_start':1,'raw_end':int(total_len),
                    'tree':obj.strip()}
        if isinstance(obj,dict):
            for key in ('tree','newick','string'):
                tree=obj.get(key,'')
                if isinstance(tree,str) and '(' in tree:
                    return {'json_partition':'single_tree','raw_start':1,'raw_end':int(total_len),
                            'tree':tree.strip()}
    raise ValueError('Cannot recover a GARD-matched single tree for zero-breakpoint model')


def snap_boundary_to_codon(bp:int, prev:int, total:int)->int:
    candidates=[]
    lo=(bp//3)*3
    hi=lo+3
    for c in {lo,hi}:
        if c>prev and c<total and c%3==0:
            candidates.append(c)
    if not candidates:
        raise ValueError(f'No valid codon boundary near breakpoint {bp}')
    return min(candidates,key=lambda c:(abs(c-bp),c))


def quote_nexus_name(s:str)->str:
    return "'"+s.replace("'","''")+"'"


def build_nexus(out:Path, seqs, spans, trees):
    nchar=len(seqs[0][1]); ntax=len(seqs)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf-8') as f:
        f.write('#NEXUS\n\nBEGIN TAXA;\n')
        f.write(f'  DIMENSIONS NTAX={ntax};\n  TAXLABELS\n')
        for name,_ in seqs: f.write(f'    {quote_nexus_name(name)}\n')
        f.write('  ;\nEND;\n\nBEGIN CHARACTERS;\n')
        f.write(f'  DIMENSIONS NCHAR={nchar};\n  FORMAT DATATYPE=DNA GAP=- MISSING=?;\n  MATRIX\n')
        for name,seq in seqs: f.write(f'    {quote_nexus_name(name)} {seq}\n')
        f.write('  ;\nEND;\n\nBEGIN ASSUMPTIONS;\n')
        for i,(a,b) in enumerate(spans,1): f.write(f'  CHARSET span_{i} = {a}-{b};\n')
        f.write('END;\n\nBEGIN TREES;\n')
        for i,t in enumerate(trees,1):
            t=t.rstrip(';')+';'
            f.write(f'  TREE tree_{i} = {t}\n')
        f.write('END;\n')


def main():
    root=root_from_env_or_file(); base=root/'06_TSWV_evolution'
    status_path=base/'tables'/'16_TSWV_GARD_final_status.tsv'
    if not status_path.exists():
        print('Missing GARD final status table; run 07O first.',file=sys.stderr); return 2
    status={r['gene']:r for r in read_tsv(status_path)}
    audits=[]; qc=[]
    outroot=base/'pipeline_results'/'05_selection'/'00_inputs'

    for gene in GENES:
        r=status.get(gene)
        if not r or r.get('status') not in VALID:
            raise SystemExit(f'{gene}: GARD status is not selection-safe: {r}')
        gd=base/'pipeline_results'/'04_GARD'/gene
        fasta=gd/f'{gene}.GARD_input.fna'; jp=gd/f'{gene}.GARD.json'
        seqs,L=read_fasta(fasta)
        if L%3!=0: raise SystemExit(f'{gene}: GARD input length {L} is not divisible by 3')
        with jp.open() as f: js=json.load(f)
        expected_bp=int(r['final_breakpoint_n'])

        # IMPORTANT: do not parse breakpoint coordinates for a zero-breakpoint
        # model. HyPhy GARD JSON variants differ here; some encode the whole
        # alignment with a one-element `bps` object.  The scientifically correct
        # selection input is simply the complete codon alignment plus the
        # GARD-fitted tree for the exact diversity-subset taxa.
        if expected_bp==0:
            parts=[single_tree_partition(js,L)]
            adjusted_breaks=[]
            use_trees=[parts[0]['tree']]
            raw_breaks=[]
            tree_source_label='GARD single-tree/baseline tree'
        else:
            parts=sorted_partition_data(js)
            if len(parts)-1 != expected_bp:
                raise SystemExit(f'{gene}: expected {expected_bp} breakpoints but JSON has {len(parts)-1}')
            raw_breaks=[p['raw_end'] for p in parts[:-1]]
            adjusted_breaks=[]; prev=0
            for bp in raw_breaks:
                adj=snap_boundary_to_codon(bp,prev,L)
                if adjusted_breaks and adj<=adjusted_breaks[-1]:
                    raise SystemExit(f'{gene}: codon-snapped breakpoints are not strictly increasing')
                adjusted_breaks.append(adj); prev=adj
            use_trees=[p['tree'] for p in parts]
            tree_source_label='GARD breakpointData partition tree'

        starts=[1]+[x+1 for x in adjusted_breaks]
        ends=adjusted_breaks+[L]
        spans=list(zip(starts,ends))
        if len(spans)!=len(use_trees):
            raise SystemExit(f'{gene}: partitions ({len(spans)}) != trees ({len(use_trees)})')
        for a,b in spans:
            if (b-a+1)%3!=0:
                raise SystemExit(f'{gene}: non-codon partition {a}-{b}')

        od=outroot/gene; nex=od/f'{gene}.selection_codon_safe.nex'
        build_nexus(nex,seqs,spans,use_trees)

        for i,(a,b) in enumerate(spans,1):
            raw=parts[i-1]
            raw_a=raw['raw_start']; raw_b=raw['raw_end']
            raw_boundary='' if i==len(spans) else raw_b
            adj_boundary='' if i==len(spans) else b
            audits.append({
                'gene':gene,'gard_status':r['status'],'partition_index':i,
                'raw_start_nt':raw_a,'raw_end_nt':raw_b,
                'adjusted_start_nt':a,'adjusted_end_nt':b,
                'codon_start':(a+2)//3,'codon_end':b//3,
                'raw_breakpoint_nt':raw_boundary,'codon_safe_breakpoint_nt':adj_boundary,
                'breakpoint_shift_nt': '' if raw_boundary=='' else int(adj_boundary)-int(raw_boundary),
                'tree_source':tree_source_label,
                'selection_nexus':str(nex),
            })
        qc.append({'gene':gene,'status':'PASS','sequence_n':len(seqs),'alignment_nt':L,
                   'partition_n':len(spans),'breakpoint_n':len(spans)-1,
                   'max_abs_breakpoint_shift_nt':max([abs(int(x['breakpoint_shift_nt'])) for x in audits if x['gene']==gene and x['breakpoint_shift_nt']!=''] or [0]),
                   'selection_nexus':str(nex)})
        print(gene,'selection input PASS',f'n={len(seqs)}',f'nt={L}',f'partitions={len(spans)}',nex)

    write_tsv(base/'tables'/'17_TSWV_selection_partition_audit.tsv',audits,list(audits[0].keys()))
    qc_path=base/'qc'/'08_TSWV_selection_input_QC.tsv'
    write_tsv(qc_path,qc,list(qc[0].keys()))
    print('Fig7 07P codon-safe selection input preparation PASS')
    return 0

if __name__=='__main__': raise SystemExit(main())
