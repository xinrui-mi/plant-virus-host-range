#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Collect 03I RNA family tree outputs, membership, and final QC."""

from __future__ import annotations
import argparse, csv
from pathlib import Path

FAMILIES = [
    "Potyviridae",
    "Betaflexiviridae",
    "Secoviridae",
    "Tombusviridae",
    "Rhabdoviridae",
    "Solemoviridae",
    "Partitiviridae",
    "Closteroviridae",
    "Mitoviridae",
    "Alphaflexiviridae",
]

def fasta_n(path: Path) -> int:
    if not path.is_file():
        return 0
    n = 0
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n

def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    args = ap.parse_args()
    root = Path(args.root).resolve()

    fig3 = root / "three_source_manuscript_analysis_v1_1" / "Fig3_host_phylogenetic_constraint"
    work = fig3 / "work" / "03I_RNA_standardized_RdRP"
    outroot = fig3 / "pipeline_results" / "03I_RNA_standardized_RdRP"
    tables = fig3 / "tables"
    qcdir = fig3 / "qc"

    summary = []
    membership_all = []
    failed = []

    for fam in FAMILIES:
        d = outroot / fam
        raw = work / f"{fam}.RdRP_core.raw.faa"
        nr95 = d / f"{fam}.RdRP_core.nr95.faa"
        filtered = d / f"{fam}.RdRP_core.nr95.trim.filtered.faa"
        tree = d / f"{fam}.RdRP_core.treefile"
        membership = d / f"{fam}.cluster95_membership.tsv"

        raw_n = fasta_n(raw)
        cluster_n = fasta_n(nr95)
        tip_n = fasta_n(filtered)
        tree_ok = tree.is_file() and tree.stat().st_size > 0
        mem_ok = membership.is_file() and membership.stat().st_size > 0
        status = "PASS" if raw_n > 0 and cluster_n > 0 and tip_n >= 4 and tree_ok and mem_ok else "FAIL"
        if status != "PASS":
            failed.append(fam)

        summary.append({
            "virus_family": fam,
            "marker": "standardized Stage-8 RdRP conserved core",
            "raw_RdRP_core_n": raw_n,
            "cluster95_n": cluster_n,
            "post_alignment_QC_tip_n": tip_n,
            "treefile": str(tree),
            "status": status,
        })

        if mem_ok:
            membership_all.extend(read_tsv(membership))

    write_tsv(
        outroot / "03I_RNA_tree_pipeline_manifest.tsv",
        [
            "virus_family", "marker", "raw_RdRP_core_n", "cluster95_n",
            "post_alignment_QC_tip_n", "treefile", "status"
        ],
        summary
    )

    mem_fields = [
        "virus_family", "cluster95_id", "representative_sequence_id",
        "member_sequence_id", "is_representative",
        "identity_to_representative_pct"
    ]
    write_tsv(
        tables / "41_03I_RNA_cluster95_membership_all.tsv",
        mem_fields,
        membership_all
    )

    write_tsv(
        tables / "42_03I_RNA_tree_tip_summary.tsv",
        [
            "virus_family", "marker", "raw_RdRP_core_n", "cluster95_n",
            "post_alignment_QC_tip_n", "treefile", "status"
        ],
        summary
    )

    qc = [
        {
            "metric": "eligible_RNA_family_n",
            "value": len(FAMILIES),
            "status": "PASS",
            "note": ",".join(FAMILIES),
        },
        {
            "metric": "families_with_PASS_tree",
            "value": sum(r["status"] == "PASS" for r in summary),
            "status": "PASS" if not failed else "FAIL",
            "note": "",
        },
        {
            "metric": "families_failed",
            "value": len(failed),
            "status": "PASS" if not failed else "FAIL",
            "note": ",".join(failed),
        },
        {
            "metric": "combined_cluster95_membership_rows",
            "value": len(membership_all),
            "status": "PASS" if membership_all else "FAIL",
            "note": str(tables / "41_03I_RNA_cluster95_membership_all.tsv"),
        },
    ]
    write_tsv(
        qcdir / "08_Fig3_03I_RNA_tree_QC.tsv",
        ["metric", "value", "status", "note"],
        qc
    )

    print("03I collector")
    for r in summary:
        print(
            r["virus_family"],
            "raw", r["raw_RdRP_core_n"],
            "nr95", r["cluster95_n"],
            "tips", r["post_alignment_QC_tip_n"],
            r["status"]
        )
    if failed:
        raise SystemExit("03I failed families: " + ", ".join(failed))
    print("Fig3 03I standardized-RdRP ten-family phylogeny pipeline PASS")

if __name__ == "__main__":
    main()
