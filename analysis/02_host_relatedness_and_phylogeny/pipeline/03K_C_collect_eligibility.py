#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse, csv

FAMS = [
    "Geminiviridae","Caulimoviridae","Alphaflexiviridae",
    "Betaflexiviridae","Closteroviridae","Mitoviridae",
    "Partitiviridae","Potyviridae","Rhabdoviridae",
    "Secoviridae","Solemoviridae","Tombusviridae"
]

def read_tsv(path):
    with open(path) as f:
        yield from csv.DictReader(f, delimiter="\t")

def read_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0].split("|",1)[0])
    return ids

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root).resolve()
    base = root / "three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    kroot = base / "pipeline_results/03K_family_specific_markers"
    outroot = kroot / "06_AA95_alignment_QC"
    plan_path = base / "pipeline_results/03J_family_specific_marker_audit/01_family_marker_plan.tsv"
    seqmaster_path = root / "final_three_source_pipeline_v1/results/20_three_source_sequence_master_FINAL.tsv"
    assoc_path = root / "final_three_source_pipeline_v1/results/20_three_source_association_master_FINAL.tsv"

    plan = {r["virus_family"]: r for r in read_tsv(plan_path)}
    seqmaster = {r["sequence_id"]: r for r in read_tsv(seqmaster_path)}

    host_fam = defaultdict(set)
    host_taxid = defaultdict(set)
    with open(assoc_path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            sid = r["sequence_id"]
            hf = (r.get("host_family") or "").strip()
            ht = (r.get("host_taxid") or "").strip()
            if hf and hf.lower() not in {"na","nan","none","unclassified"}:
                host_fam[sid].add(hf)
            if ht and ht.lower() not in {"na","nan","none"}:
                host_taxid[sid].add(ht)

    # marker maps: use corrected Secoviridae map, standard 03K-B maps for others.
    marker_species = {}
    marker_seq_n = {}
    marker_type = {}

    for fam in FAMS:
        if fam == "Secoviridae":
            mp = kroot / "05_Secoviridae_corrected/Secoviridae.ProPol_CG_to_GDD.corrected.tsv"
            rows = list(read_tsv(mp))
            marker_seq_n[fam] = len(rows)
            marker_species[fam] = len({
                (r.get("virus_species") or "").strip()
                for r in rows
                if (r.get("virus_species") or "").strip() not in {"","NA","nan","None"}
            })
            marker_type[fam] = "Pro-Pol CG-to-GDD"
        else:
            mp = kroot / f"03_final_marker_recovery/{fam}.marker_map.tsv"
            rows = list(read_tsv(mp))
            marker_seq_n[fam] = len(rows)
            marker_species[fam] = len({
                (r.get("virus_species") or "").strip()
                for r in rows
                if (r.get("virus_species") or "").strip() not in {"","NA","nan","None"}
            })
            vals = sorted({
                (r.get("marker_type") or "").strip()
                for r in rows if (r.get("marker_type") or "").strip()
            })
            marker_type[fam] = ";".join(vals)

    eligibility = []

    for fam in FAMS:
        D = outroot / fam
        mem_path = D / f"{fam}.AA95_membership.tsv"
        qc_faa = D / f"{fam}.AA95.mafft.trim.filtered.faa"
        aa95_faa = D / f"{fam}.AA95.faa"
        qc_report = D / f"{fam}.alignment_sequence_QC.tsv"

        candidate_species_n = int(plan[fam]["distinct_virus_species_n"])

        if not (mem_path.exists() and qc_faa.exists() and aa95_faa.exists()):
            eligibility.append({
                "virus_family": fam,
                "candidate_species_level_label_n": candidate_species_n,
                "marker_type": marker_type.get(fam,"NA"),
                "marker_recovered_sequence_n": marker_seq_n.get(fam,0),
                "marker_recovered_species_level_label_n": marker_species.get(fam,0),
                "marker_species_recovery_pct_of_candidate": (
                    f"{100*marker_species.get(fam,0)/candidate_species_n:.2f}"
                    if candidate_species_n else "NA"
                ),
                "AA95_representative_n": "NA",
                "alignment_QC_tree_tip_n": "NA",
                "host_annotated_tree_tip_n": "NA",
                "host_family_coverage_n": "NA",
                "distinct_host_spectrum_n": "NA",
                "technical_eligibility": "FAIL",
                "exclusion_reason": "03K_C_family_output_missing_or_job_failed",
            })
            continue

        members = list(read_tsv(mem_path))
        by_rep = defaultdict(list)
        cluster_id = {}
        cluster_n = {}
        for r in members:
            rep = r["representative_sequence_id"]
            by_rep[rep].append(r["member_sequence_id"])
            cluster_id[rep] = r["cdhit95_cluster_id"]
            cluster_n[rep] = int(r["cluster_member_n"])

        aa95_ids = set(read_fasta_ids(aa95_faa))
        kept_ids = set(read_fasta_ids(qc_faa))

        # Tree-tip host spectrum is aggregated across every member of the AA95 cluster.
        tip_rows = []
        all_hf = set()
        host_tip_n = 0
        spectra = set()

        for rep in sorted(kept_ids):
            mids = by_rep.get(rep, [rep])
            hfs = set()
            hts = set()
            for sid in mids:
                hfs |= host_fam.get(sid, set())
                hts |= host_taxid.get(sid, set())

            if hfs:
                host_tip_n += 1
                all_hf |= hfs
                spectra.add(tuple(sorted(hfs)))

            sm = seqmaster.get(rep, {})
            tip_rows.append({
                "virus_family": fam,
                "cdhit95_cluster_id": cluster_id.get(rep,""),
                "representative_sequence_id": rep,
                "cluster_member_n": cluster_n.get(rep, len(mids)),
                "member_sequence_ids": ";".join(sorted(mids)),
                "virus_taxid_representative": (sm.get("virus_taxid") or "").strip(),
                "virus_genus_representative": (sm.get("virus_genus") or "").strip(),
                "virus_species_representative": (sm.get("virus_species") or "").strip(),
                "host_annotated": "YES" if hfs else "NO",
                "host_taxid_n": len(hts),
                "host_taxids": ";".join(sorted(hts)),
                "host_family_n": len(hfs),
                "host_families": ";".join(sorted(hfs)),
            })

        tip_path = D / f"{fam}.tree_tip_host_spectrum.tsv"
        fields = [
            "virus_family","cdhit95_cluster_id","representative_sequence_id",
            "cluster_member_n","member_sequence_ids",
            "virus_taxid_representative","virus_genus_representative",
            "virus_species_representative","host_annotated",
            "host_taxid_n","host_taxids","host_family_n","host_families"
        ]
        with open(tip_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(tip_rows)

        reasons = []
        # Freeze the original data-rich criterion at the marker-recovered level.
        if marker_species.get(fam,0) < 100:
            reasons.append("marker_recovered_species_level_labels_lt100")
        # Reuse the old frozen technical minimum for phylogeny construction.
        if len(aa95_ids) < 10:
            reasons.append("AA95_representatives_lt10")
        if len(kept_ids) < 10:
            reasons.append("alignment_QC_tree_tips_lt10")
        # Minimal host-spectrum analyzability requirements.
        if host_tip_n < 10:
            reasons.append("host_annotated_tree_tips_lt10")
        if len(all_hf) < 2:
            reasons.append("host_family_coverage_lt2")
        if len(spectra) < 2:
            reasons.append("distinct_host_spectra_lt2")

        eligibility.append({
            "virus_family": fam,
            "candidate_species_level_label_n": candidate_species_n,
            "marker_type": marker_type.get(fam,"NA"),
            "marker_recovered_sequence_n": marker_seq_n.get(fam,0),
            "marker_recovered_species_level_label_n": marker_species.get(fam,0),
            "marker_species_recovery_pct_of_candidate": (
                f"{100*marker_species.get(fam,0)/candidate_species_n:.2f}"
                if candidate_species_n else "NA"
            ),
            "AA95_representative_n": len(aa95_ids),
            "alignment_QC_tree_tip_n": len(kept_ids),
            "host_annotated_tree_tip_n": host_tip_n,
            "host_family_coverage_n": len(all_hf),
            "distinct_host_spectrum_n": len(spectra),
            "technical_eligibility": "PASS" if not reasons else "FAIL",
            "exclusion_reason": "NONE" if not reasons else ";".join(reasons),
        })

    out = outroot / "03K_C_final_eligibility.tsv"
    fields = list(eligibility[0].keys())
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(eligibility)

    print(f"Eligibility table: {out}")
    print(f"PASS={sum(r['technical_eligibility']=='PASS' for r in eligibility)}/{len(eligibility)}")

if __name__ == "__main__":
    main()
