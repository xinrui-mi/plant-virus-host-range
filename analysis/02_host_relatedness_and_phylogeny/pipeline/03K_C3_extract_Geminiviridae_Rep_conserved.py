#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse, csv, math, statistics

IE = 1e-5
FLANK = 10

def read_fasta(path):
    d={}
    h=None; b=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if h is not None: d[h]="".join(b)
                h=line[1:].split()[0]; b=[]
            else:
                b.append(line)
    if h is not None: d[h]="".join(b)
    return d

def read_tsv(path):
    with open(path) as f:
        yield from csv.DictReader(f, delimiter="\t")

def parse_domtbl(path):
    hits=defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            p=line.split()
            if len(p)<22: continue
            hits[p[3]].append({
                "model":p[0],
                "i_e":float(p[12]),
                "hmm_from":int(p[15]),
                "hmm_to":int(p[16]),
                "ali_from":int(p[17]),
                "ali_to":int(p[18]),
            })
    return hits

def q(vals,p):
    vals=sorted(vals)
    if not vals: return None
    if len(vals)==1: return float(vals[0])
    x=(len(vals)-1)*p
    lo=math.floor(x); hi=math.ceil(x)
    if lo==hi: return float(vals[lo])
    return vals[lo]*(hi-x)+vals[hi]*(x-lo)

def best_pairs(hs):
    a=[h for h in hs if h["model"]=="Gemini_AL1" and h["i_e"]<=IE]
    b=[h for h in hs if h["model"]=="Gemini_AL1_M" and h["i_e"]<=IE]
    out=[]
    for x in a:
        for y in b:
            # Require the canonical ordered N-terminal Rep architecture.
            # Allow a small overlap but reject obviously unrelated repeated hits.
            if y["ali_from"] < x["ali_from"]:
                continue
            if y["ali_from"] - x["ali_to"] > 120:
                continue
            span = y["ali_to"] - x["ali_from"] + 1
            if 140 <= span <= 360:
                score = -math.log10(max(x["i_e"],1e-300)) - math.log10(max(y["i_e"],1e-300))
                out.append((score, x, y))
    return sorted(out, reverse=True, key=lambda z:z[0])

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args=ap.parse_args()

    root=Path(args.project_root).resolve()
    base=root/"three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    kroot=base/"pipeline_results/03K_family_specific_markers"

    faa=kroot/"01_orf_recovery/Geminiviridae.candidate_orfs.faa"
    meta=kroot/"01_orf_recovery/Geminiviridae.candidate_orfs.tsv"
    dom=kroot/"02_pfam/Geminiviridae.Pfam.domtblout"
    seqmaster=root/"final_three_source_pipeline_v1/results/20_three_source_sequence_master_FINAL.tsv"
    assoc=root/"final_three_source_pipeline_v1/results/20_three_source_association_master_FINAL.tsv"

    outdir=kroot/"07_Geminiviridae_Rep_conserved"
    outdir.mkdir(parents=True, exist_ok=True)

    seqs=read_fasta(faa)
    m={r["orf_id"]:r for r in read_tsv(meta)}
    hits=parse_domtbl(dom)
    sm={r["sequence_id"]:r for r in read_tsv(seqmaster)}

    hostfam=defaultdict(set)
    with open(assoc) as f:
        rd=csv.DictReader(f, delimiter="\t")
        for r in rd:
            hf=(r.get("host_family") or "").strip()
            if hf and hf.lower() not in {"na","nan","none","unclassified"}:
                hostfam[r["sequence_id"]].add(hf)

    by_sid=defaultdict(list)
    for oid,row in m.items():
        sid=row["sequence_id"]
        seq=seqs.get(oid)
        if not seq: continue
        pairs=best_pairs(hits.get(oid,[]))
        if not pairs: continue
        score,a1,am=pairs[0]
        start=max(1,a1["ali_from"]-FLANK)
        end=min(len(seq),am["ali_to"]+FLANK)
        marker=seq[start-1:end]
        if not (160 <= len(marker) <= 380):
            continue
        by_sid[sid].append((score,oid,start,end,marker,a1,am))

    selected=[]
    for sid,cands in by_sid.items():
        cands.sort(reverse=True,key=lambda z:z[0])
        score,oid,start,end,marker,a1,am=cands[0]
        x=sm.get(sid,{})
        selected.append({
            "sequence_id":sid,
            "orf_id":oid,
            "marker_start_aa":start,
            "marker_end_aa":end,
            "marker_length_aa":len(marker),
            "Gemini_AL1_iE":f"{a1['i_e']:.6g}",
            "Gemini_AL1_M_iE":f"{am['i_e']:.6g}",
            "virus_taxid":(x.get("virus_taxid") or "").strip(),
            "virus_genus":(x.get("virus_genus") or "").strip(),
            "virus_species":(x.get("virus_species") or "").strip(),
            "host_annotated":"YES" if hostfam.get(sid) else "NO",
            "host_family_n":len(hostfam.get(sid,set())),
            "_marker":marker
        })

    fasta=outdir/"Geminiviridae.Rep_conserved.raw.faa"
    with open(fasta,"w") as f:
        for r in selected:
            f.write(f">{r['sequence_id']}|marker=Rep_conserved_AL1_to_AL1M|orf={r['orf_id']}\n")
            s=r["_marker"]
            for i in range(0,len(s),80): f.write(s[i:i+80]+"\n")

    fields=[k for k in selected[0] if k!="_marker"] if selected else [
        "sequence_id","orf_id","marker_start_aa","marker_end_aa","marker_length_aa",
        "Gemini_AL1_iE","Gemini_AL1_M_iE","virus_taxid","virus_genus","virus_species",
        "host_annotated","host_family_n"
    ]
    with open(outdir/"Geminiviridae.Rep_conserved.map.tsv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader()
        for r in selected:
            w.writerow({k:r[k] for k in fields})

    species={(r["virus_species"] or "").strip() for r in selected if (r["virus_species"] or "").strip()}
    genera={(r["virus_genus"] or "").strip() for r in selected if (r["virus_genus"] or "").strip()}
    lengths=[r["marker_length_aa"] for r in selected]
    hf=set()
    for r in selected: hf |= hostfam.get(r["sequence_id"],set())

    with open(outdir/"Geminiviridae.Rep_conserved.summary.tsv","w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["metric","value"])
        w.writerow(["marker_definition","Gemini_AL1_to_Gemini_AL1_M_continuous_region_plus_10aa_flanks"])
        w.writerow(["marker_sequence_n",len(selected)])
        w.writerow(["species_level_label_n",len(species)])
        w.writerow(["named_genus_n",len(genera)])
        w.writerow(["host_family_n",len(hf)])
        w.writerow(["host_annotated_sequence_n",sum(r["host_annotated"]=="YES" for r in selected)])
        if lengths:
            w.writerow(["marker_length_min",min(lengths)])
            w.writerow(["marker_length_q05",f"{q(lengths,.05):.1f}"])
            w.writerow(["marker_length_median",f"{statistics.median(lengths):.1f}"])
            w.writerow(["marker_length_q95",f"{q(lengths,.95):.1f}"])
            w.writerow(["marker_length_max",max(lengths)])

    print("03K-C3 Geminiviridae Rep-conserved extraction PASS")
    print(f"marker_sequences={len(selected)}")
    print(f"species_level_labels={len(species)}")
    print(f"named_genera={len(genera)}")
    print(f"host_families={len(hf)}")
    if lengths:
        print(f"length=min:{min(lengths)} q05:{q(lengths,.05):.1f} median:{statistics.median(lengths):.1f} q95:{q(lengths,.95):.1f} max:{max(lengths)}")

if __name__=="__main__":
    main()
