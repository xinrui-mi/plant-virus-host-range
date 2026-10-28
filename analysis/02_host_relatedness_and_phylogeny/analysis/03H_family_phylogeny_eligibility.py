#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03H_family_phylogeny_eligibility.py

Audit virus-family eligibility for expanded family-level phylogeny–host-spectrum analyses.

Design principles
-----------------
1. Read only the frozen Stage-20 three-source masters and the frozen Stage-8 RdRP core FASTA.
2. Do NOT rerun HMM searches and do NOT extend RdRP regions upstream/downstream.
3. Count formal ICTV species separately from sequence abundance.
4. For RNA families, quantify how many sequences already have a standardized Stage-8 RdRP core
   and how many of those sequences have plant-host annotations.
5. Do NOT decide final tree inclusion here. The script only produces an objective audit table.

Outputs
-------
Fig3_host_phylogenetic_constraint/tables/
  30_03H_family_phylogeny_eligibility_all.tsv
  31_03H_family_phylogeny_candidates_ge100_species.tsv
  32_03H_family_source_composition.tsv
  33_03H_RNA_family_RdRP_core_summary.tsv
Fig3_host_phylogenetic_constraint/qc/
  05_Fig3_03H_eligibility_QC.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple


MISSING_STRINGS = {"", "na", "n/a", "nan", "none", "null", "."}
EXISTING_FOUR = {"Geminiviridae", "Caulimoviridae", "Potyviridae", "Tombusviridae"}


def present(x: object) -> bool:
    if x is None:
        return False
    s = str(x).strip()
    return bool(s) and s.lower() not in MISSING_STRINGS


def norm(x: object) -> str:
    return str(x).strip() if present(x) else ""


def safe_pct(num: int, den: int) -> str:
    if den <= 0:
        return "NA"
    return f"{100.0 * num / den:.4f}"


def open_tsv(path: Path):
    return path.open("r", encoding="utf-8-sig", newline="")


def fasta_headers(path: Path) -> Iterator[str]:
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip()


def resolve_core_sequence_id(header: str, all_sequence_ids: Set[str]) -> Optional[str]:
    """Match Stage-8 RdRP FASTA headers using the same conservative rule as frozen 03D."""
    token = header.split()[0] if header.split() else ""
    candidates = [token]
    if "|" in token:
        candidates.append(token.split("|", 1)[0])
    for candidate in candidates:
        if candidate in all_sequence_ids:
            return candidate
    return None


def write_rows(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def add_set(d: Dict[str, Set[str]], key: str, value: str) -> None:
    if present(value):
        d[key].add(value)


def classify_marker_route(rna_n: int, dna_n: int, other_n: int) -> str:
    if rna_n > 0 and dna_n == 0 and other_n == 0:
        return "RNA_standardized_RdRP_core"
    if dna_n > 0 and rna_n == 0 and other_n == 0:
        return "non_RNA_family_specific_marker_required"
    if rna_n > 0 and (dna_n > 0 or other_n > 0):
        return "mixed_or_special_family_manual_review"
    return "non_RNA_or_unresolved_manual_review"


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit family-level phylogeny eligibility from frozen three-source data.")
    parser.add_argument(
        "--root",
        default=None,
        help="final_three_source root. If omitted, infer from this script location.",
    )
    parser.add_argument(
        "--min-species",
        type=int,
        default=100,
        help="Distinct virus_species threshold used only to define the teacher-requested candidate review set (default: 100).",
    )
    args = parser.parse_args()

    script = Path(__file__).resolve()
    if args.root:
        root = Path(args.root).resolve()
    else:
        # <root>/three_source_manuscript_analysis_v1_1/Fig3.../analysis/03H...
        try:
            root = script.parents[3]
        except IndexError:
            raise SystemExit("Cannot infer project root; supply --root explicitly.")

    pipeline_results = root / "final_three_source_pipeline_v1" / "results"
    analysis_root = root / "three_source_manuscript_analysis_v1_1"
    fig3 = analysis_root / "Fig3_host_phylogenetic_constraint"
    table_dir = fig3 / "tables"
    qc_dir = fig3 / "qc"
    table_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    seq_path = pipeline_results / "20_three_source_sequence_master_FINAL.tsv"
    assoc_path = pipeline_results / "20_three_source_association_master_FINAL.tsv"
    core_fasta = pipeline_results / "08_RdRP_core.faa"

    for p in (seq_path, assoc_path, core_fasta):
        if not p.is_file():
            raise SystemExit(f"Missing required frozen input: {p}")

    # ------------------------------------------------------------------
    # 1) Sequence master: family composition, formal species, source, etc.
    # ------------------------------------------------------------------
    seq_meta: Dict[str, Dict[str, str]] = {}
    family_sequence_ids: Dict[str, Set[str]] = defaultdict(set)
    family_species: Dict[str, Set[str]] = defaultdict(set)
    family_species_rank_matched: Dict[str, Set[str]] = defaultdict(set)
    family_taxids: Dict[str, Set[str]] = defaultdict(set)
    family_genera: Dict[str, Set[str]] = defaultdict(set)
    family_types: Dict[str, Set[str]] = defaultdict(set)
    family_genomes: Dict[str, Set[str]] = defaultdict(set)
    family_sources: Dict[str, Counter] = defaultdict(Counter)
    family_rna_ids: Dict[str, Set[str]] = defaultdict(set)
    family_dna_ids: Dict[str, Set[str]] = defaultdict(set)
    family_other_type_ids: Dict[str, Set[str]] = defaultdict(set)
    master_rdrp_cluster_ids: Set[str] = set()

    sequence_row_n = 0
    duplicate_sequence_id_n = 0

    with open_tsv(seq_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sequence_id", "virus_family", "virus_species", "virus_taxid", "virus_type", "genome", "rdrp_operational_cluster_id"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"Sequence master missing required columns: {', '.join(missing)}")
        source_col = "source_dataset" if "source_dataset" in (reader.fieldnames or []) else (
            "source_label_4way" if "source_label_4way" in (reader.fieldnames or []) else None
        )
        rank_col = "virus_ictv_match_rank" if "virus_ictv_match_rank" in (reader.fieldnames or []) else None
        genus_col = "virus_genus" if "virus_genus" in (reader.fieldnames or []) else None

        for row in reader:
            sequence_row_n += 1
            sid = norm(row.get("sequence_id"))
            fam = norm(row.get("virus_family"))
            if not sid:
                continue
            if sid in seq_meta:
                duplicate_sequence_id_n += 1
                continue

            source = norm(row.get(source_col)) if source_col else ""
            vtype = norm(row.get("virus_type"))
            genome = norm(row.get("genome"))
            species = norm(row.get("virus_species"))
            taxid = norm(row.get("virus_taxid"))
            genus = norm(row.get(genus_col)) if genus_col else ""
            rank = norm(row.get(rank_col)).lower() if rank_col else ""
            rdrp_cluster = norm(row.get("rdrp_operational_cluster_id"))

            seq_meta[sid] = {
                "family": fam,
                "source": source,
                "virus_type": vtype,
                "genome": genome,
                "virus_species": species,
                "virus_taxid": taxid,
                "virus_genus": genus,
                "ictv_rank": rank,
                "rdrp_cluster": rdrp_cluster,
            }

            if not fam:
                continue

            family_sequence_ids[fam].add(sid)
            if species:
                # Main species count follows the final frozen virus_species field.
                # This directly answers the review question: how many distinct virus species
                # are represented in each family. The stricter ICTV match-rank count is
                # retained separately as an audit metric.
                family_species[fam].add(species)
                if not rank_col or rank == "species":
                    family_species_rank_matched[fam].add(species)
            if taxid:
                family_taxids[fam].add(taxid)
            if genus:
                family_genera[fam].add(genus)
            if vtype:
                family_types[fam].add(vtype)
            if genome:
                family_genomes[fam].add(genome)
            family_sources[fam][source if source else "UNKNOWN"] += 1

            vt = vtype.upper()
            if vt == "RNA" or "RNA" in genome.upper():
                family_rna_ids[fam].add(sid)
            elif vt == "DNA" or "DNA" in genome.upper():
                family_dna_ids[fam].add(sid)
            else:
                family_other_type_ids[fam].add(sid)

            if rdrp_cluster:
                master_rdrp_cluster_ids.add(sid)

    if duplicate_sequence_id_n:
        raise SystemExit(f"Sequence master contains {duplicate_sequence_id_n} duplicate sequence_id values; aborting audit.")

    all_sequence_ids = set(seq_meta)

    # ---------------------------------------------------------------
    # 2) Frozen Stage-8 standardized RdRP core FASTA; no HMM rerun.
    # ---------------------------------------------------------------
    rdrp_core_ids: Set[str] = set()
    rdrp_core_duplicate_id_n = 0
    rdrp_core_unmatched_headers: List[str] = []
    rdrp_core_fasta_record_n = 0

    for header in fasta_headers(core_fasta):
        rdrp_core_fasta_record_n += 1
        sid = resolve_core_sequence_id(header, all_sequence_ids)
        if sid is None:
            rdrp_core_unmatched_headers.append(header.split()[0] if header.split() else header)
            continue
        if sid in rdrp_core_ids:
            rdrp_core_duplicate_id_n += 1
        rdrp_core_ids.add(sid)

    family_core_ids: Dict[str, Set[str]] = defaultdict(set)
    family_core_sources: Dict[str, Counter] = defaultdict(Counter)
    for sid in rdrp_core_ids:
        meta = seq_meta.get(sid)
        if not meta:
            continue
        fam = meta["family"]
        if not fam:
            continue
        family_core_ids[fam].add(sid)
        family_core_sources[fam][meta["source"] if meta["source"] else "UNKNOWN"] += 1

    # ---------------------------------------------------------------
    # 3) Association master: host coverage overall and within RdRP core.
    # ---------------------------------------------------------------
    family_assoc_seq: Dict[str, Set[str]] = defaultdict(set)
    family_host_annot_seq: Dict[str, Set[str]] = defaultdict(set)
    family_host_taxids: Dict[str, Set[str]] = defaultdict(set)
    family_host_families: Dict[str, Set[str]] = defaultdict(set)

    family_core_host_annot_seq: Dict[str, Set[str]] = defaultdict(set)
    family_core_host_taxids: Dict[str, Set[str]] = defaultdict(set)
    family_core_host_families: Dict[str, Set[str]] = defaultdict(set)

    association_row_n = 0
    association_sequence_missing_master_n = 0

    with open_tsv(assoc_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sequence_id", "virus_family", "host_taxid", "host_family"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"Association master missing required columns: {', '.join(missing)}")

        for row in reader:
            association_row_n += 1
            sid = norm(row.get("sequence_id"))
            if not sid:
                continue
            meta = seq_meta.get(sid)
            if meta is None:
                association_sequence_missing_master_n += 1
                continue
            fam = meta["family"] or norm(row.get("virus_family"))
            if not fam:
                continue

            family_assoc_seq[fam].add(sid)
            htax = norm(row.get("host_taxid"))
            hfam = norm(row.get("host_family"))
            if hfam:
                family_host_annot_seq[fam].add(sid)
                family_host_families[fam].add(hfam)
            if htax:
                family_host_taxids[fam].add(htax)

            if sid in rdrp_core_ids:
                if hfam:
                    family_core_host_annot_seq[fam].add(sid)
                    family_core_host_families[fam].add(hfam)
                if htax:
                    family_core_host_taxids[fam].add(htax)

    # ---------------------------------------------------------------
    # 4) Family audit tables.
    # ---------------------------------------------------------------
    all_families = sorted(family_sequence_ids)
    audit_rows: List[Dict[str, object]] = []

    for fam in all_families:
        total_n = len(family_sequence_ids[fam])
        species_n = len(family_species[fam])
        species_rank_matched_n = len(family_species_rank_matched[fam])
        taxid_n = len(family_taxids[fam])
        genus_n = len(family_genera[fam])
        rna_n = len(family_rna_ids[fam])
        dna_n = len(family_dna_ids[fam])
        other_n = len(family_other_type_ids[fam])
        core_n = len(family_core_ids[fam])
        core_host_n = len(family_core_host_annot_seq[fam])
        cluster_n = sum(1 for sid in family_sequence_ids[fam] if sid in master_rdrp_cluster_ids)

        ge100 = species_n >= args.min_species
        marker_route = classify_marker_route(rna_n, dna_n, other_n)
        if ge100:
            if marker_route == "RNA_standardized_RdRP_core":
                review_group = "GE100_RNA_RdRP_CORE_REVIEW"
            elif marker_route == "non_RNA_family_specific_marker_required":
                review_group = "GE100_NON_RNA_MARKER_REVIEW"
            else:
                review_group = "GE100_SPECIAL_MANUAL_REVIEW"
        else:
            review_group = "BELOW_SPECIES_THRESHOLD"

        audit_rows.append({
            "virus_family": fam,
            "existing_four_family": "YES" if fam in EXISTING_FOUR else "NO",
            "distinct_virus_species_n": species_n,
            "ICTV_species_rank_matched_n": species_rank_matched_n,
            "formal_virus_taxid_n": taxid_n,
            "virus_genus_n": genus_n,
            "total_sequence_n": total_n,
            "RNA_sequence_n": rna_n,
            "DNA_sequence_n": dna_n,
            "other_or_unresolved_type_sequence_n": other_n,
            "virus_type_set": ";".join(sorted(family_types[fam])),
            "genome_set": ";".join(sorted(family_genomes[fam])),
            "source_dataset_n": len(family_sources[fam]),
            "association_sequence_n": len(family_assoc_seq[fam]),
            "host_annotated_sequence_n": len(family_host_annot_seq[fam]),
            "host_taxid_n": len(family_host_taxids[fam]),
            "host_family_n": len(family_host_families[fam]),
            "rdrp_cluster_assigned_sequence_n": cluster_n,
            "rdrp_core_fasta_sequence_n": core_n,
            "rdrp_core_recovery_pct_of_RNA_sequences": safe_pct(core_n, rna_n),
            "rdrp_core_host_annotated_sequence_n": core_host_n,
            "rdrp_core_host_annotation_pct": safe_pct(core_host_n, core_n),
            "rdrp_core_host_taxid_n": len(family_core_host_taxids[fam]),
            "rdrp_core_host_family_n": len(family_core_host_families[fam]),
            f"virus_species_ge_{args.min_species}": "YES" if ge100 else "NO",
            "marker_route_for_review": marker_route,
            "eligibility_review_group": review_group,
        })

    audit_rows.sort(key=lambda r: (-int(r["distinct_virus_species_n"]), -int(r["total_sequence_n"]), str(r["virus_family"])))

    fields = [
        "virus_family", "existing_four_family",
        "distinct_virus_species_n", "ICTV_species_rank_matched_n", "formal_virus_taxid_n", "virus_genus_n",
        "total_sequence_n", "RNA_sequence_n", "DNA_sequence_n", "other_or_unresolved_type_sequence_n",
        "virus_type_set", "genome_set", "source_dataset_n",
        "association_sequence_n", "host_annotated_sequence_n", "host_taxid_n", "host_family_n",
        "rdrp_cluster_assigned_sequence_n", "rdrp_core_fasta_sequence_n", "rdrp_core_recovery_pct_of_RNA_sequences",
        "rdrp_core_host_annotated_sequence_n", "rdrp_core_host_annotation_pct",
        "rdrp_core_host_taxid_n", "rdrp_core_host_family_n",
        f"virus_species_ge_{args.min_species}", "marker_route_for_review", "eligibility_review_group",
    ]

    out_all = table_dir / "30_03H_family_phylogeny_eligibility_all.tsv"
    out_candidates = table_dir / "31_03H_family_phylogeny_candidates_ge100_species.tsv"
    out_sources = table_dir / "32_03H_family_source_composition.tsv"
    out_rna = table_dir / "33_03H_RNA_family_RdRP_core_summary.tsv"
    out_qc = qc_dir / "05_Fig3_03H_eligibility_QC.tsv"

    write_rows(out_all, fields, audit_rows)

    candidate_rows = [r for r in audit_rows if r[f"virus_species_ge_{args.min_species}"] == "YES"]
    write_rows(out_candidates, fields, candidate_rows)

    source_rows: List[Dict[str, object]] = []
    for fam in all_families:
        sources = sorted(set(family_sources[fam]) | set(family_core_sources[fam]))
        for src in sources:
            source_rows.append({
                "virus_family": fam,
                "source_dataset": src,
                "all_sequence_n": family_sources[fam][src],
                "rdrp_core_sequence_n": family_core_sources[fam][src],
            })
    source_rows.sort(key=lambda r: (str(r["virus_family"]), str(r["source_dataset"])))
    write_rows(out_sources, ["virus_family", "source_dataset", "all_sequence_n", "rdrp_core_sequence_n"], source_rows)

    rna_fields = [
        "virus_family", "existing_four_family", "distinct_virus_species_n", "formal_virus_taxid_n",
        "RNA_sequence_n", "rdrp_cluster_assigned_sequence_n", "rdrp_core_fasta_sequence_n",
        "rdrp_core_recovery_pct_of_RNA_sequences", "rdrp_core_host_annotated_sequence_n",
        "rdrp_core_host_annotation_pct", "rdrp_core_host_taxid_n", "rdrp_core_host_family_n",
        f"virus_species_ge_{args.min_species}", "eligibility_review_group",
    ]
    rna_rows = [r for r in audit_rows if int(r["RNA_sequence_n"]) > 0]
    write_rows(out_rna, rna_fields, rna_rows)

    # ---------------------------------------------------------------
    # 5) QC: cross-check Stage-8 core FASTA vs final RdRP65 assignment.
    # ---------------------------------------------------------------
    master_missing_core = master_rdrp_cluster_ids - rdrp_core_ids
    core_without_master_cluster = rdrp_core_ids - master_rdrp_cluster_ids

    qc_rows = [
        {"metric": "sequence_master_rows", "value": sequence_row_n, "status": "PASS", "note": str(seq_path)},
        {"metric": "sequence_master_unique_sequence_ids", "value": len(all_sequence_ids), "status": "PASS", "note": ""},
        {"metric": "association_master_rows", "value": association_row_n, "status": "PASS", "note": str(assoc_path)},
        {"metric": "association_sequence_ids_missing_from_sequence_master", "value": association_sequence_missing_master_n,
         "status": "PASS" if association_sequence_missing_master_n == 0 else "WARN", "note": "Expected 0 for frozen Stage-20 masters."},
        {"metric": "rdrp_core_fasta_records", "value": rdrp_core_fasta_record_n, "status": "PASS", "note": str(core_fasta)},
        {"metric": "rdrp_core_mapped_unique_sequence_ids", "value": len(rdrp_core_ids), "status": "PASS", "note": "No HMM rerun; frozen Stage-8 core only."},
        {"metric": "rdrp_core_unmatched_fasta_headers", "value": len(rdrp_core_unmatched_headers),
         "status": "PASS" if len(rdrp_core_unmatched_headers) == 0 else "WARN",
         "note": ";".join(rdrp_core_unmatched_headers[:10])},
        {"metric": "rdrp_core_duplicate_mapped_sequence_ids", "value": rdrp_core_duplicate_id_n,
         "status": "PASS" if rdrp_core_duplicate_id_n == 0 else "WARN", "note": ""},
        {"metric": "final_master_rdrp65_assigned_sequence_ids", "value": len(master_rdrp_cluster_ids), "status": "PASS", "note": "rdrp_operational_cluster_id present"},
        {"metric": "rdrp65_assigned_missing_from_core_fasta", "value": len(master_missing_core),
         "status": "PASS" if len(master_missing_core) == 0 else "WARN", "note": "Cross-check only; investigate if nonzero."},
        {"metric": "core_fasta_ids_without_final_rdrp65_assignment", "value": len(core_without_master_cluster),
         "status": "PASS" if len(core_without_master_cluster) == 0 else "WARN", "note": "Cross-check only; investigate if nonzero."},
        {"metric": f"families_with_virus_species_ge_{args.min_species}", "value": len(candidate_rows), "status": "PASS", "note": "Candidate review set, not final eligibility."},
        {"metric": f"RNA_families_with_virus_species_ge_{args.min_species}",
         "value": sum(1 for r in candidate_rows if int(r["RNA_sequence_n"]) > 0), "status": "PASS", "note": ""},
        {"metric": f"non_RNA_or_special_families_with_virus_species_ge_{args.min_species}",
         "value": sum(1 for r in candidate_rows if int(r["RNA_sequence_n"]) == 0), "status": "PASS", "note": "Need family-specific marker review."},
    ]
    write_rows(out_qc, ["metric", "value", "status", "note"], qc_rows)

    print("03H family phylogeny eligibility audit PASS")
    print(f"Root: {root}")
    print(f"Families audited: {len(audit_rows)}")
    print(f"Families with >= {args.min_species} distinct virus species: {len(candidate_rows)}")
    print(f"Frozen RdRP core FASTA records: {rdrp_core_fasta_record_n}")
    print(f"Mapped unique RdRP core sequence IDs: {len(rdrp_core_ids)}")
    print("\nCandidate families:")
    for r in candidate_rows:
        print(
            f"  {r['virus_family']}: species={r['distinct_virus_species_n']}, "
            f"seq={r['total_sequence_n']}, RNA={r['RNA_sequence_n']}, "
            f"RdRP_core={r['rdrp_core_fasta_sequence_n']}, "
            f"RdRP_core_host_annot={r['rdrp_core_host_annotated_sequence_n']}, "
            f"host_families(core)={r['rdrp_core_host_family_n']}, "
            f"review={r['eligibility_review_group']}"
        )
    print("\nOutputs:")
    for p in (out_all, out_candidates, out_sources, out_rna, out_qc):
        print(f"  {p}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
