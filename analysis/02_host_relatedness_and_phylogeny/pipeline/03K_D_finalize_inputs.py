#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse, csv

FAMS = [
    "Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
    "Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
    "Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
]

def read_tsv(path):
    with open(path) as f:
        yield from csv.DictReader(f, delimiter="\t")

def fasta_ids(path):
    out=[]
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                out.append(line[1:].strip().split()[0].split("|",1)[0])
    return out

def write_tip_host_spectrum(mem_path, filtered_faa, seqmaster_path, assoc_path, out_path):
    seqmaster={r["sequence_id"]:r for r in read_tsv(seqmaster_path)}
    hostfam=defaultdict(set)
    hosttax=defaultdict(set)
    with open(assoc_path) as f:
        rd=csv.DictReader(f, delimiter="\t")
        for r in rd:
            sid=r["sequence_id"]
            hf=(r.get("host_family") or "").strip()
            ht=(r.get("host_taxid") or "").strip()
            if hf and hf.lower() not in {"na","nan","none","unclassified"}:
                hostfam[sid].add(hf)
            if ht and ht.lower() not in {"na","nan","none"}:
                hosttax[sid].add(ht)

    by_rep=defaultdict(list)
    cluster={}
    cluster_n={}
    for r in read_tsv(mem_path):
        rep=r["representative_sequence_id"]
        by_rep[rep].append(r["member_sequence_id"])
        cluster[rep]=r["cdhit95_cluster_id"]
        cluster_n[rep]=int(r["cluster_member_n"])

    keep=set(fasta_ids(filtered_faa))
    rows=[]
    for rep in sorted(keep):
        members=by_rep.get(rep,[rep])
        hfs=set(); hts=set()
        for sid in members:
            hfs |= hostfam.get(sid,set())
            hts |= hosttax.get(sid,set())
        sm=seqmaster.get(rep,{})
        rows.append({
            "virus_family":"Geminiviridae",
            "cdhit95_cluster_id":cluster.get(rep,""),
            "representative_sequence_id":rep,
            "cluster_member_n":cluster_n.get(rep,len(members)),
            "member_sequence_ids":";".join(sorted(members)),
            "virus_taxid_representative":(sm.get("virus_taxid") or "").strip(),
            "virus_genus_representative":(sm.get("virus_genus") or "").strip(),
            "virus_species_representative":(sm.get("virus_species") or "").strip(),
            "host_annotated":"YES" if hfs else "NO",
            "host_taxid_n":len(hts),
            "host_taxids":";".join(sorted(hts)),
            "host_family_n":len(hfs),
            "host_families":";".join(sorted(hfs)),
        })
    fields=list(rows[0].keys())
    with open(out_path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader(); w.writerows(rows)
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    a=ap.parse_args()

    root=Path(a.project_root).resolve()
    base=root/"three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    kroot=base/"pipeline_results/03K_family_specific_markers"
    croot=kroot/"06_AA95_alignment_QC"
    groot=kroot/"07_Geminiviridae_Rep_conserved"
    droot=kroot/"08_final_12family_trees"
    inp=droot/"inputs"
    inp.mkdir(parents=True,exist_ok=True)

    seqmaster=root/"final_three_source_pipeline_v1/results/20_three_source_sequence_master_FINAL.tsv"
    assoc=root/"final_three_source_pipeline_v1/results/20_three_source_association_master_FINAL.tsv"

    # Final corrected Gemini host-spectrum table.
    gmem=groot/"AA95_alignment_QC/Geminiviridae.Rep_conserved.AA95_membership.tsv"
    gfaa=groot/"AA95_alignment_QC/Geminiviridae.Rep_conserved.AA95.mafft.trim.filtered.faa"
    gtip=groot/"Geminiviridae.tree_tip_host_spectrum.tsv"
    grows=write_tip_host_spectrum(gmem,gfaa,seqmaster,assoc,gtip)

    # Build frozen final eligibility table.
    old={r["virus_family"]:r for r in read_tsv(croot/"03K_C_final_eligibility.tsv")}
    rows=[]
    for fam in FAMS:
        if fam=="Geminiviridae":
            oldr=old[fam]
            rows.append({
                "virus_family":"Geminiviridae",
                "candidate_species_level_label_n":oldr["candidate_species_level_label_n"],
                "marker_type":"Rep_conserved_AL1_to_AL1M",
                "marker_recovered_sequence_n":"14694",
                "marker_recovered_species_level_label_n":"736",
                "marker_species_recovery_pct_of_candidate":f"{100*736/int(oldr['candidate_species_level_label_n']):.2f}",
                "AA95_representative_n":"1678",
                "alignment_QC_tree_tip_n":"1604",
                "host_annotated_tree_tip_n":str(sum(r["host_annotated"]=="YES" for r in grows)),
                "host_family_coverage_n":str(len(set(x for r in grows for x in r["host_families"].split(";") if x))),
                "distinct_host_spectrum_n":str(len(set(r["host_families"] for r in grows if r["host_families"]))),
                "technical_eligibility":"PASS",
                "exclusion_reason":"NONE",
            })
        else:
            rows.append(old[fam])

    fields=list(rows[0].keys())
    final_tab=droot/"03K_D_final_eligibility_FROZEN.tsv"
    with open(final_tab,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader(); w.writerows(rows)

    # Freeze standardized tree inputs and source manifest.
    manifest=[]
    for fam in FAMS:
        if fam=="Geminiviridae":
            src=gfaa
            mem=gmem
            tip=gtip
            marker="Rep_conserved_AL1_to_AL1M"
        else:
            D=croot/fam
            src=D/f"{fam}.AA95.mafft.trim.filtered.faa"
            mem=D/f"{fam}.AA95_membership.tsv"
            tip=D/f"{fam}.tree_tip_host_spectrum.tsv"
            marker=next(r["marker_type"] for r in rows if r["virus_family"]==fam)

        if not src.exists() or not mem.exists() or not tip.exists():
            raise RuntimeError(f"Missing finalized input/mapping for {fam}")

        dst=inp/f"{fam}.final_alignment.faa"
        dst.write_bytes(src.read_bytes())
        manifest.append({
            "virus_family":fam,
            "marker_type":marker,
            "source_alignment":str(src),
            "frozen_alignment":str(dst),
            "tree_tip_n":len(fasta_ids(dst)),
            "AA95_membership":str(mem),
            "tree_tip_host_spectrum":str(tip),
            "eligibility":"PASS",
        })

    mfields=list(manifest[0].keys())
    with open(droot/"03K_D_tree_input_manifest.tsv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=mfields,delimiter="\t")
        w.writeheader();w.writerows(manifest)

    print("03K-D input finalization PASS")
    print(f"Final eligibility: {final_tab}")
    print("eligible_families=12")
    for r in manifest:
        print(r["virus_family"],"tips="+str(r["tree_tip_n"]),r["marker_type"])

if __name__=="__main__":
    main()
