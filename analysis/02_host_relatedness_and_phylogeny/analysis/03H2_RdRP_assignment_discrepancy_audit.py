#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03H2_RdRP_assignment_discrepancy_audit.py

Purpose
-------
Audit the 1,199 sequences that have a frozen Stage-8 RdRP core but lack a
non-empty rdrp_operational_cluster_id in the final Stage-20 sequence master.

This script:
1. Reuses the frozen Stage-8 RdRP core FASTA; it does NOT rerun HMM.
2. Maps all Stage-8 RdRP core records back to the Stage-20 final sequence master.
3. Identifies core-positive sequences lacking final rdrp_operational_cluster_id.
4. Summarizes those sequences by source dataset and virus family.
5. Cross-checks whether they are present in:
   - 11_RdRP_selected_membership.tsv
   - 09_RdRP_id65_membership.tsv
6. Produces a small QC report to distinguish a downstream field-propagation
   issue from an earlier RdRP detection / clustering issue.

No upstream result is modified.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


MISSING = {"", "na", "n/a", "nan", "none", "null", "."}


def present(x) -> bool:
    if x is None:
        return False
    s = str(x).strip()
    return bool(s) and s.lower() not in MISSING


def norm(x) -> str:
    return str(x).strip() if present(x) else ""


def open_tsv(path: Path):
    return path.open("r", encoding="utf-8-sig", newline="")


def fasta_headers(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip()


def resolve_core_sequence_id(header: str, master_ids: Set[str]) -> Optional[str]:
    token = header.split()[0] if header.split() else ""
    candidates = [token]
    if "|" in token:
        candidates.append(token.split("|", 1)[0])
    for x in candidates:
        if x in master_ids:
            return x
    return None


def write_rows(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def detect_sequence_id_column(path: Path, valid_ids: Set[str]) -> Tuple[str, Dict[str, str], int]:
    """
    Detect which column of a membership table contains sequence IDs.
    Returns: (column_name, sequence_id -> raw row key/cluster-ish value map, matched_n)
    The map value is a compact text representation of non-ID fields for audit only.
    """
    with open_tsv(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields:
            raise RuntimeError(f"No header detected in {path}")

        rows = list(reader)

    preferred = [
        "sequence_id", "seq_id", "member_id", "member", "accession",
        "accession_id", "query_id", "query", "id"
    ]
    candidates = [c for c in preferred if c in fields] + [c for c in fields if c not in preferred]

    best_col = ""
    best_n = -1
    for col in candidates:
        n = 0
        for row in rows:
            v = norm(row.get(col))
            if v in valid_ids:
                n += 1
        if n > best_n:
            best_n = n
            best_col = col

    if best_n <= 0:
        raise RuntimeError(
            f"Could not detect a sequence-ID column in {path}. Columns: {fields}"
        )

    out = {}
    for row in rows:
        sid = norm(row.get(best_col))
        if sid not in valid_ids:
            continue
        # Keep a compact record of other fields; useful if we later need to inspect
        # how the selected membership was encoded.
        extras = []
        for col in fields:
            if col == best_col:
                continue
            v = norm(row.get(col))
            if v:
                extras.append(f"{col}={v}")
        out[sid] = ";".join(extras[:8])

    return best_col, out, best_n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=None,
                    help="final_three_source root; inferred from script location if omitted")
    args = ap.parse_args()

    script = Path(__file__).resolve()
    if args.root:
        root = Path(args.root).resolve()
    else:
        root = script.parents[3]

    res = root / "final_three_source_pipeline_v1" / "results"
    fig3 = root / "three_source_manuscript_analysis_v1_1" / "Fig3_host_phylogenetic_constraint"
    tables = fig3 / "tables"
    qcdir = fig3 / "qc"
    tables.mkdir(parents=True, exist_ok=True)
    qcdir.mkdir(parents=True, exist_ok=True)

    seq_master = res / "20_three_source_sequence_master_FINAL.tsv"
    core_faa = res / "08_RdRP_core.faa"
    selected_membership = res / "11_RdRP_selected_membership.tsv"
    id65_membership = res / "09_RdRP_id65_membership.tsv"
    threshold_file = res / "10_RdRP_threshold_selected.txt"

    required = [seq_master, core_faa]
    for p in required:
        if not p.is_file():
            raise SystemExit(f"Missing required input: {p}")

    # ------------------------------------------------------------
    # 1. Final master
    # ------------------------------------------------------------
    meta: Dict[str, Dict[str, str]] = {}
    with open_tsv(seq_master) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or [])
        req = {"sequence_id", "source_dataset", "virus_family",
               "virus_genus", "virus_species",
               "rdrp_operational_cluster_id"}
        miss = sorted(req - fields)
        if miss:
            raise SystemExit("Final sequence master missing columns: " + ", ".join(miss))

        for row in reader:
            sid = norm(row.get("sequence_id"))
            if not sid:
                continue
            meta[sid] = {
                "source_dataset": norm(row.get("source_dataset")),
                "virus_family": norm(row.get("virus_family")),
                "virus_genus": norm(row.get("virus_genus")),
                "virus_species": norm(row.get("virus_species")),
                "virus_taxid": norm(row.get("virus_taxid")),
                "virus_type": norm(row.get("virus_type")),
                "genome": norm(row.get("genome")),
                "rdrp_operational_cluster_id": norm(row.get("rdrp_operational_cluster_id")),
            }

    master_ids = set(meta)

    # ------------------------------------------------------------
    # 2. Stage-8 core FASTA mapping
    # ------------------------------------------------------------
    core_ids: Set[str] = set()
    unmatched_headers: List[str] = []
    core_record_n = 0

    for h in fasta_headers(core_faa):
        core_record_n += 1
        sid = resolve_core_sequence_id(h, master_ids)
        if sid is None:
            unmatched_headers.append(h)
        else:
            core_ids.add(sid)

    final_cluster_ids = {
        sid for sid, m in meta.items()
        if present(m["rdrp_operational_cluster_id"])
    }

    discrepancy_ids = sorted(core_ids - final_cluster_ids)

    # ------------------------------------------------------------
    # 3. Cross-check frozen membership tables
    # ------------------------------------------------------------
    selected_map = {}
    selected_id_col = ""
    selected_matched_n = 0
    if selected_membership.is_file():
        selected_id_col, selected_map, selected_matched_n = detect_sequence_id_column(
            selected_membership, master_ids
        )

    id65_map = {}
    id65_id_col = ""
    id65_matched_n = 0
    if id65_membership.is_file():
        id65_id_col, id65_map, id65_matched_n = detect_sequence_id_column(
            id65_membership, master_ids
        )

    selected_threshold = ""
    if threshold_file.is_file():
        selected_threshold = threshold_file.read_text(
            encoding="utf-8", errors="replace"
        ).strip()

    # ------------------------------------------------------------
    # 4. Per-sequence audit
    # ------------------------------------------------------------
    detail_rows = []
    for sid in discrepancy_ids:
        m = meta[sid]
        detail_rows.append({
            "sequence_id": sid,
            "source_dataset": m["source_dataset"],
            "virus_family": m["virus_family"],
            "virus_genus": m["virus_genus"],
            "virus_species": m["virus_species"],
            "virus_taxid": m["virus_taxid"],
            "virus_type": m["virus_type"],
            "genome": m["genome"],
            "stage8_RdRP_core_present": "YES",
            "final_rdrp_operational_cluster_id": m["rdrp_operational_cluster_id"],
            "in_11_RdRP_selected_membership": "YES" if sid in selected_map else "NO",
            "in_09_RdRP_id65_membership": "YES" if sid in id65_map else "NO",
            "selected_membership_fields": selected_map.get(sid, ""),
            "id65_membership_fields": id65_map.get(sid, ""),
        })

    detail_fields = [
        "sequence_id", "source_dataset", "virus_family", "virus_genus",
        "virus_species", "virus_taxid", "virus_type", "genome",
        "stage8_RdRP_core_present", "final_rdrp_operational_cluster_id",
        "in_11_RdRP_selected_membership", "in_09_RdRP_id65_membership",
        "selected_membership_fields", "id65_membership_fields"
    ]

    # ------------------------------------------------------------
    # 5. Summaries
    # ------------------------------------------------------------
    by_source = Counter()
    by_family = Counter()
    by_source_family = Counter()
    membership_pattern = Counter()

    for sid in discrepancy_ids:
        m = meta[sid]
        src = m["source_dataset"] or "NA"
        fam = m["virus_family"] or "NA"
        by_source[src] += 1
        by_family[fam] += 1
        by_source_family[(src, fam)] += 1
        membership_pattern[
            (
                "YES" if sid in selected_map else "NO",
                "YES" if sid in id65_map else "NO"
            )
        ] += 1

    source_rows = [
        {"source_dataset": k, "sequence_n": v,
         "percent_of_discrepancy": f"{100*v/len(discrepancy_ids):.4f}" if discrepancy_ids else "NA"}
        for k, v in sorted(by_source.items(), key=lambda x: (-x[1], x[0]))
    ]

    family_rows = [
        {"virus_family": k, "sequence_n": v,
         "percent_of_discrepancy": f"{100*v/len(discrepancy_ids):.4f}" if discrepancy_ids else "NA"}
        for k, v in sorted(by_family.items(), key=lambda x: (-x[1], x[0]))
    ]

    sf_rows = [
        {"source_dataset": src, "virus_family": fam, "sequence_n": n,
         "percent_of_discrepancy": f"{100*n/len(discrepancy_ids):.4f}" if discrepancy_ids else "NA"}
        for (src, fam), n in sorted(by_source_family.items(), key=lambda x: (-x[1], x[0][0], x[0][1]))
    ]

    membership_rows = [
        {
            "in_11_RdRP_selected_membership": a,
            "in_09_RdRP_id65_membership": b,
            "sequence_n": n,
            "percent_of_discrepancy": f"{100*n/len(discrepancy_ids):.4f}" if discrepancy_ids else "NA"
        }
        for (a, b), n in sorted(membership_pattern.items(), key=lambda x: -x[1])
    ]

    # ------------------------------------------------------------
    # 6. QC
    # ------------------------------------------------------------
    all_discrepancy_in_selected = (
        bool(discrepancy_ids) and
        len(selected_map) > 0 and
        all(sid in selected_map for sid in discrepancy_ids)
    )
    all_discrepancy_in_id65 = (
        bool(discrepancy_ids) and
        len(id65_map) > 0 and
        all(sid in id65_map for sid in discrepancy_ids)
    )

    if all_discrepancy_in_selected and all_discrepancy_in_id65:
        interpretation = (
            "All Stage-8-core / Stage-20-cluster discrepancies are still present in both "
            "the selected RdRP membership and the 65% membership table. This supports a "
            "downstream field-propagation/join issue rather than loss at HMM detection or RdRP65 clustering."
        )
    elif all_discrepancy_in_selected:
        interpretation = (
            "All discrepancy sequences remain in the selected RdRP membership, but not all "
            "are recovered in the explicit 65% membership cross-check. Inspect membership encoding "
            "before assigning cause."
        )
    else:
        interpretation = (
            "At least some discrepancy sequences are absent from the selected RdRP membership. "
            "The difference cannot yet be attributed solely to downstream Stage-20 field propagation."
        )

    qc_rows = [
        {"metric": "stage20_sequence_master_rows", "value": len(master_ids), "status": "PASS", "note": str(seq_master)},
        {"metric": "stage8_RdRP_core_fasta_records", "value": core_record_n, "status": "PASS", "note": str(core_faa)},
        {"metric": "stage8_RdRP_core_mapped_unique_ids", "value": len(core_ids),
         "status": "PASS" if not unmatched_headers else "WARN",
         "note": f"unmatched_headers={len(unmatched_headers)}"},
        {"metric": "stage20_final_rdrp_cluster_nonempty_ids", "value": len(final_cluster_ids), "status": "PASS", "note": ""},
        {"metric": "core_positive_without_stage20_cluster_id", "value": len(discrepancy_ids),
         "status": "WARN" if discrepancy_ids else "PASS", "note": ""},
        {"metric": "11_selected_membership_exists", "value": str(selected_membership.is_file()),
         "status": "PASS" if selected_membership.is_file() else "WARN",
         "note": f"id_column={selected_id_col}; matched_master_ids={selected_matched_n}"},
        {"metric": "09_id65_membership_exists", "value": str(id65_membership.is_file()),
         "status": "PASS" if id65_membership.is_file() else "WARN",
         "note": f"id_column={id65_id_col}; matched_master_ids={id65_matched_n}"},
        {"metric": "selected_threshold_text", "value": selected_threshold or "NA", "status": "PASS", "note": str(threshold_file)},
        {"metric": "all_discrepancy_ids_in_11_selected_membership",
         "value": str(all_discrepancy_in_selected),
         "status": "PASS" if all_discrepancy_in_selected else "WARN", "note": ""},
        {"metric": "all_discrepancy_ids_in_09_id65_membership",
         "value": str(all_discrepancy_in_id65),
         "status": "PASS" if all_discrepancy_in_id65 else "WARN", "note": ""},
        {"metric": "interpretation", "value": interpretation, "status": "INFO", "note": ""},
    ]

    # ------------------------------------------------------------
    # 7. Write
    # ------------------------------------------------------------
    p_detail = tables / "34_03H2_core_without_final_rdrp65_assignment.tsv"
    p_source = tables / "35_03H2_discrepancy_by_source.tsv"
    p_family = tables / "36_03H2_discrepancy_by_family.tsv"
    p_sf = tables / "37_03H2_discrepancy_by_source_family.tsv"
    p_mem = tables / "38_03H2_membership_crosscheck.tsv"
    p_qc = qcdir / "06_Fig3_03H2_RdRP_assignment_discrepancy_QC.tsv"

    write_rows(p_detail, detail_fields, detail_rows)
    write_rows(p_source, ["source_dataset", "sequence_n", "percent_of_discrepancy"], source_rows)
    write_rows(p_family, ["virus_family", "sequence_n", "percent_of_discrepancy"], family_rows)
    write_rows(p_sf, ["source_dataset", "virus_family", "sequence_n", "percent_of_discrepancy"], sf_rows)
    write_rows(p_mem, [
        "in_11_RdRP_selected_membership", "in_09_RdRP_id65_membership",
        "sequence_n", "percent_of_discrepancy"
    ], membership_rows)
    write_rows(p_qc, ["metric", "value", "status", "note"], qc_rows)

    print("03H2 RdRP assignment discrepancy audit finished")
    print(f"Stage-8 core mapped IDs: {len(core_ids)}")
    print(f"Stage-20 non-empty rdrp_operational_cluster_id: {len(final_cluster_ids)}")
    print(f"Discrepancy IDs: {len(discrepancy_ids)}")
    print("\nBy source:")
    for r in source_rows:
        print(f"  {r['source_dataset']}: {r['sequence_n']} ({r['percent_of_discrepancy']}%)")
    print("\nTop families:")
    for r in family_rows[:20]:
        print(f"  {r['virus_family']}: {r['sequence_n']}")
    print("\nMembership cross-check:")
    for r in membership_rows:
        print(
            f"  selected={r['in_11_RdRP_selected_membership']}, "
            f"id65={r['in_09_RdRP_id65_membership']}: {r['sequence_n']}"
        )
    print("\nInterpretation:")
    print(" ", interpretation)
    print("\nOutputs:")
    for p in [p_detail, p_source, p_family, p_sf, p_mem, p_qc]:
        print(" ", p)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
