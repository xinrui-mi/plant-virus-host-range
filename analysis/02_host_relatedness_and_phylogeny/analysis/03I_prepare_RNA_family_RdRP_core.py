#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03I_prepare_RNA_family_RdRP_core.py

Prepare standardized Stage-8 RdRP conserved-core FASTA files for the 10
predefined eligible RNA virus families.

Important:
- Reads the frozen Stage-8 08_RdRP_core.faa.
- Does NOT rerun HMM.
- Does NOT extend upstream/downstream of the HMM-defined conserved core.
- Family assignment is taken from the frozen Stage-20 final sequence master.
- FASTA headers are normalized to sequence_id only, which makes downstream
  CD-HIT membership and tree-tip mapping deterministic.
"""

from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
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

MISSING = {"", "na", "n/a", "nan", "none", "null", "."}

def present(x):
    if x is None:
        return False
    s = str(x).strip()
    return bool(s) and s.lower() not in MISSING

def norm(x):
    return str(x).strip() if present(x) else ""

def read_fasta(path: Path):
    header = None
    seq = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)

def resolve_sequence_id(header: str, master_ids: set[str]) -> str | None:
    token = header.split()[0] if header.split() else ""
    candidates = [token]
    if "|" in token:
        candidates.append(token.split("|", 1)[0])
    for x in candidates:
        if x in master_ids:
            return x
    return None

def write_tsv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    args = ap.parse_args()
    root = Path(args.root).resolve()

    results = root / "final_three_source_pipeline_v1" / "results"
    fig3 = root / "three_source_manuscript_analysis_v1_1" / "Fig3_host_phylogenetic_constraint"
    work = fig3 / "work" / "03I_RNA_standardized_RdRP"
    tables = fig3 / "tables"
    qcdir = fig3 / "qc"
    work.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    qcdir.mkdir(parents=True, exist_ok=True)

    seq_master = results / "20_three_source_sequence_master_FINAL.tsv"
    core_faa = results / "08_RdRP_core.faa"
    if not seq_master.is_file():
        raise SystemExit(f"Missing {seq_master}")
    if not core_faa.is_file():
        raise SystemExit(f"Missing {core_faa}")

    meta = {}
    with seq_master.open("r", encoding="utf-8-sig", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        required = {"sequence_id", "source_dataset", "virus_family", "virus_type", "genome"}
        miss = required - set(rd.fieldnames or [])
        if miss:
            raise SystemExit("Missing columns in final sequence master: " + ", ".join(sorted(miss)))
        for row in rd:
            sid = norm(row.get("sequence_id"))
            if not sid:
                continue
            meta[sid] = {
                "source_dataset": norm(row.get("source_dataset")),
                "virus_family": norm(row.get("virus_family")),
                "virus_type": norm(row.get("virus_type")),
                "genome": norm(row.get("genome")),
                "virus_genus": norm(row.get("virus_genus")),
                "virus_species": norm(row.get("virus_species")),
                "virus_taxid": norm(row.get("virus_taxid")),
            }

    master_ids = set(meta)
    handles = {}
    for fam in FAMILIES:
        p = work / f"{fam}.RdRP_core.raw.faa"
        handles[fam] = p.open("w", encoding="utf-8")

    counts = Counter()
    source_counts = Counter()
    unmatched = []
    duplicated_ids = Counter()
    seen = set()
    total_core_records = 0

    try:
        for header, seq in read_fasta(core_faa):
            total_core_records += 1
            sid = resolve_sequence_id(header, master_ids)
            if sid is None:
                unmatched.append(header)
                continue
            duplicated_ids[sid] += 1
            if sid in seen:
                continue
            seen.add(sid)

            fam = meta[sid]["virus_family"]
            if fam not in handles:
                continue
            handles[fam].write(f">{sid}\n{seq}\n")
            counts[fam] += 1
            source_counts[(fam, meta[sid]["source_dataset"] or "NA")] += 1
    finally:
        for h in handles.values():
            h.close()

    summary = []
    for fam in FAMILIES:
        summary.append({
            "virus_family": fam,
            "stage8_RdRP_core_sequence_n": counts[fam],
            "raw_fasta": str(work / f"{fam}.RdRP_core.raw.faa"),
        })
    write_tsv(
        tables / "39_03I_RNA_family_RdRP_core_input_summary.tsv",
        ["virus_family", "stage8_RdRP_core_sequence_n", "raw_fasta"],
        summary
    )

    src_rows = []
    for (fam, src), n in sorted(source_counts.items()):
        src_rows.append({"virus_family": fam, "source_dataset": src, "sequence_n": n})
    write_tsv(
        tables / "40_03I_RNA_family_RdRP_core_source_composition.tsv",
        ["virus_family", "source_dataset", "sequence_n"],
        src_rows
    )

    qc_rows = [
        {
            "metric": "stage8_RdRP_core_records",
            "value": total_core_records,
            "status": "PASS",
            "note": str(core_faa),
        },
        {
            "metric": "stage8_core_unmatched_headers",
            "value": len(unmatched),
            "status": "PASS" if not unmatched else "FAIL",
            "note": "Expected 0",
        },
        {
            "metric": "stage8_core_duplicate_sequence_ids",
            "value": sum(1 for v in duplicated_ids.values() if v > 1),
            "status": "PASS" if all(v == 1 for v in duplicated_ids.values()) else "FAIL",
            "note": "Expected 0",
        },
        {
            "metric": "eligible_RNA_family_n",
            "value": len(FAMILIES),
            "status": "PASS",
            "note": ",".join(FAMILIES),
        },
        {
            "metric": "eligible_RNA_core_sequence_total",
            "value": sum(counts.values()),
            "status": "PASS",
            "note": "Sum across 10 eligible RNA families",
        },
        {
            "metric": "families_with_zero_core",
            "value": sum(1 for fam in FAMILIES if counts[fam] == 0),
            "status": "PASS" if all(counts[fam] > 0 for fam in FAMILIES) else "FAIL",
            "note": ",".join(fam for fam in FAMILIES if counts[fam] == 0),
        },
    ]
    write_tsv(
        qcdir / "07_Fig3_03I_RNA_input_QC.tsv",
        ["metric", "value", "status", "note"],
        qc_rows
    )

    if unmatched:
        (work / "unmatched_stage8_headers.txt").write_text(
            "\n".join(unmatched) + "\n", encoding="utf-8"
        )

    print("03I RNA input preparation PASS")
    for fam in FAMILIES:
        print(f"{fam}\t{counts[fam]}")
    print("TOTAL", sum(counts.values()))

if __name__ == "__main__":
    main()
