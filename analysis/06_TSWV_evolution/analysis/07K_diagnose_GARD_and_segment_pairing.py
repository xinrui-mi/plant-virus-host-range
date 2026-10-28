#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07K_diagnose_GARD_and_segment_pairing.py

One-shot diagnostic for the TSWV case study (Results 3.7):
1) Diagnose the latest per-gene GARD failure from internal HyPhy logs/results.
2) Audit whether TSWV L/M/S segments can be reliably paired by isolate for
   the segment-tree comparison.

Run from anywhere inside the analysis tree:
    python analysis/06_TSWV_evolution/analysis/07K_diagnose_GARD_and_segment_pairing.py

Outputs:
    06_TSWV_evolution/tables/08_GARD_failure_diagnosis.tsv
    06_TSWV_evolution/tables/09_TSWV_segment_pairing_audit.tsv
    06_TSWV_evolution/qc/04_TSWV_segment_pairing_QC.tsv
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import defaultdict, Counter
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("ERROR: Biopython is required (Bio.SeqIO). Install it in the same Python environment used to run the released analysis scripts (see environment/python_version.txt).")

GENES = ["RdRp", "GnGc", "NSm", "N", "NSs"]
GENE_TO_SEGMENT = {"RdRp": "L", "GnGc": "M", "NSm": "M", "N": "S", "NSs": "S"}
EXPECTED_GENES = {"L": {"RdRp"}, "M": {"GnGc", "NSm"}, "S": {"N", "NSs"}}


def find_analysis_root() -> Path:
    # Script is expected at ROOT/06_TSWV_evolution/analysis/
    p = Path(__file__).resolve()
    candidate = p.parents[2]
    if (candidate / "06_TSWV_evolution").is_dir():
        return candidate
    # fallback: walk upward
    for parent in [Path.cwd().resolve(), *Path.cwd().resolve().parents]:
        if (parent / "06_TSWV_evolution").is_dir():
            return parent
    raise SystemExit("ERROR: Cannot locate the analysis root (the directory containing 06_TSWV_evolution/).")


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def norm_text(x: str) -> str:
    x = re.sub(r"\s+", " ", (x or "").strip())
    return x


def norm_isolate(x: str) -> str:
    x = norm_text(x)
    x = re.sub(r"^(isolate|strain)\s*[:=]\s*", "", x, flags=re.I)
    return x.casefold()


def first(q, keys):
    for k in keys:
        v = q.get(k, [])
        if v:
            return norm_text(str(v[0]))
    return ""


def diagnose_gard(base: Path):
    logdir = base / "logs"
    gardroot = base / "pipeline_results" / "04_GARD"
    script_path = base / "pipeline" / "07H_run_GARD_one_gene.sh"
    script_text = script_path.read_text(encoding="utf-8", errors="replace") if script_path.exists() else ""
    tolerance_in_script = "TOLERATE_NUMERICAL_ERRORS=1" in script_text

    rows = []
    for gene in GENES:
        stderr_p = logdir / f"{gene}.GARD.stderr.log"
        stdout_p = logdir / f"{gene}.GARD.stdout.md"
        outdir = gardroot / gene
        json_p = outdir / f"{gene}.GARD.json"
        nex_p = outdir / f"{gene}.best-gard.nex"
        sub_p = outdir / f"{gene}.GARD_input.fna"

        err = stderr_p.read_text(encoding="utf-8", errors="replace") if stderr_p.exists() else ""
        out = stdout_p.read_text(encoding="utf-8", errors="replace") if stdout_p.exists() else ""
        text = err + "\n" + out

        flags = []
        primary = "unknown_or_no_internal_log"
        if "ComputeBranchCache" in text or "numerical instability" in text.lower():
            primary = "hyphy_numerical_instability_ComputeBranchCache"
            flags.append("ComputeBranchCache")
        if "MPI_ABORT" in text:
            flags.append("MPI_ABORT_secondary")
        if "OpenFabrics" in text or "openib" in text.lower() or "mlx5_" in text:
            flags.append("OpenFabrics_warning")
        if "killed by signal 15" in text.lower() or "signal 15" in text.lower():
            flags.append("SIGTERM_secondary")
        if "unknown option" in text.lower() or "unrecognized" in text.lower():
            primary = "command_line_option_error"
            flags.append("CLI_option_error")
        if "TOLERATE_NUMERICAL_ERRORS" in text:
            flags.append("HyPhy_mentions_tolerance")
        if "GARD PASS" in out:
            primary = "PASS"

        json_ok = False
        json_size = json_p.stat().st_size if json_p.exists() else 0
        if json_size > 0:
            try:
                json.loads(json_p.read_text(encoding="utf-8", errors="replace"))
                json_ok = True
            except Exception:
                pass
        nex_size = nex_p.stat().st_size if nex_p.exists() else 0

        rows.append({
            "gene": gene,
            "diagnosis": primary,
            "secondary_flags": ";".join(flags),
            "tolerance_present_in_07H_script": str(tolerance_in_script).upper(),
            "stderr_bytes": stderr_p.stat().st_size if stderr_p.exists() else 0,
            "stdout_bytes": stdout_p.stat().st_size if stdout_p.exists() else 0,
            "subset_bytes": sub_p.stat().st_size if sub_p.exists() else 0,
            "GARD_json_bytes": json_size,
            "GARD_json_valid": str(json_ok).upper(),
            "best_gard_nex_bytes": nex_size,
            "stderr_path": str(stderr_p),
        })

    out_tsv = base / "tables" / "08_GARD_failure_diagnosis.tsv"
    write_tsv(out_tsv, rows, list(rows[0].keys()))

    print("\n===== GARD DIAGNOSIS =====")
    for r in rows:
        print(
            f"{r['gene']:4s}  {r['diagnosis']}"
            f"  [{r['secondary_flags']}]"
            f"  JSON={r['GARD_json_bytes']}  NEX={r['best_gard_nex_bytes']}"
        )
    print(f"07H script contains TOLERATE_NUMERICAL_ERRORS=1: {tolerance_in_script}")
    print(f"GARD diagnosis table: {out_tsv}")
    return rows


def load_strict_gene_map(base: Path):
    p = base / "pipeline_results" / "01_gene_extraction" / "metadata" / "strict_qc.tsv"
    if not p.exists():
        raise SystemExit(f"ERROR: Missing {p}")
    gene_by_acc = defaultdict(set)
    with p.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("strict_pass") or "").upper() == "TRUE":
                gene = r.get("assigned_gene", "")
                acc = r.get("sequence_id", "")
                if gene in GENE_TO_SEGMENT and acc:
                    gene_by_acc[acc].add(gene)
    return gene_by_acc


def segment_pairing_audit(base: Path):
    gene_by_acc = load_strict_gene_map(base)
    gb = base / "work" / "TSWV_formal_current.gb"
    if not gb.exists():
        raise SystemExit(f"ERROR: Missing {gb}")

    record_rows = []
    for rec in SeqIO.parse(str(gb), "genbank"):
        if rec.id not in gene_by_acc:
            continue

        source_q = {}
        for feat in rec.features:
            if feat.type == "source":
                source_q = feat.qualifiers
                break

        isolate = first(source_q, ["isolate"])
        strain = first(source_q, ["strain"])
        key_raw = isolate or strain
        key_source = "isolate" if isolate else ("strain" if strain else "")
        iso_key = norm_isolate(key_raw)

        country = first(source_q, ["geo_loc_name", "country"])
        date = first(source_q, ["collection_date"])
        host = first(source_q, ["host", "lab_host"])
        segment_qual = first(source_q, ["segment"])

        genes = set(gene_by_acc[rec.id])
        inferred_segments = sorted({GENE_TO_SEGMENT[g] for g in genes})

        record_rows.append({
            "sequence_id": rec.id,
            "isolate_raw": key_raw,
            "isolate_key": iso_key,
            "isolate_key_source": key_source,
            "country": country,
            "collection_date": date,
            "host": host,
            "segment_qualifier": segment_qual,
            "inferred_segments": ",".join(inferred_segments),
            "genes": ",".join(sorted(genes)),
            "description": rec.description or "",
        })

    by_iso = defaultdict(list)
    no_isolate = []
    for r in record_rows:
        if r["isolate_key"]:
            by_iso[r["isolate_key"]].append(r)
        else:
            no_isolate.append(r)

    audit = []
    for iso_key, grp in sorted(by_iso.items()):
        seg_to_acc = defaultdict(set)
        seg_to_genes = defaultdict(set)
        countries, dates, hosts, raw_names = set(), set(), set(), set()

        for r in grp:
            raw_names.add(r["isolate_raw"])
            if r["country"]:
                countries.add(r["country"])
            if r["collection_date"]:
                dates.add(r["collection_date"])
            if r["host"]:
                hosts.add(r["host"])
            genes = set(filter(None, r["genes"].split(",")))
            for gene in genes:
                seg = GENE_TO_SEGMENT[gene]
                seg_to_acc[seg].add(r["sequence_id"])
                seg_to_genes[seg].add(gene)

        segs = set(seg_to_acc)
        all_three = segs == {"L", "M", "S"}
        complete_L = EXPECTED_GENES["L"].issubset(seg_to_genes["L"])
        complete_M = EXPECTED_GENES["M"].issubset(seg_to_genes["M"])
        complete_S = EXPECTED_GENES["S"].issubset(seg_to_genes["S"])
        complete_coding_segments = complete_L and complete_M and complete_S

        multi_accession_segments = [s for s in "LMS" if len(seg_to_acc[s]) > 1]
        metadata_conflict = len(countries) > 1 or len(dates) > 1

        if not all_three:
            status = "incomplete_segments"
        elif metadata_conflict:
            status = "LMS_metadata_conflict_review"
        elif multi_accession_segments:
            status = "LMS_multiple_accessions_review"
        elif complete_coding_segments:
            status = "LMS_unique_high_confidence"
        else:
            status = "LMS_unique_but_gene_incomplete"

        audit.append({
            "isolate_key": iso_key,
            "isolate_raw_names": ";".join(sorted(raw_names)),
            "record_n": len(grp),
            "segment_n": len(segs),
            "segments": ",".join(sorted(segs)),
            "L_accession_n": len(seg_to_acc["L"]),
            "M_accession_n": len(seg_to_acc["M"]),
            "S_accession_n": len(seg_to_acc["S"]),
            "L_accessions": ";".join(sorted(seg_to_acc["L"])),
            "M_accessions": ";".join(sorted(seg_to_acc["M"])),
            "S_accessions": ";".join(sorted(seg_to_acc["S"])),
            "L_genes": ",".join(sorted(seg_to_genes["L"])),
            "M_genes": ",".join(sorted(seg_to_genes["M"])),
            "S_genes": ",".join(sorted(seg_to_genes["S"])),
            "complete_L_coding": str(complete_L).upper(),
            "complete_M_coding": str(complete_M).upper(),
            "complete_S_coding": str(complete_S).upper(),
            "complete_LMS_coding": str(complete_coding_segments).upper(),
            "country_set": ";".join(sorted(countries)),
            "collection_date_set": ";".join(sorted(dates)),
            "host_set": ";".join(sorted(hosts)),
            "metadata_conflict_country_or_date": str(metadata_conflict).upper(),
            "multiple_accession_segments": ",".join(multi_accession_segments),
            "pairing_status": status,
        })

    out_audit = base / "tables" / "09_TSWV_segment_pairing_audit.tsv"
    fields = list(audit[0].keys()) if audit else []
    if audit:
        write_tsv(out_audit, audit, fields)

    status_counts = Counter(r["pairing_status"] for r in audit)
    lms_n = sum(r["segments"] == "L,M,S" for r in audit)
    two_plus_n = sum(r["segment_n"] >= 2 for r in audit)
    high_conf_n = status_counts.get("LMS_unique_high_confidence", 0)
    lms_complete_coding_n = sum(
        r["segments"] == "L,M,S" and r["complete_LMS_coding"] == "TRUE" for r in audit
    )

    qc_rows = [
        {"metric": "strict_gene_bearing_accession_n", "value": len(record_rows)},
        {"metric": "accession_with_isolate_or_strain_n", "value": len(record_rows) - len(no_isolate)},
        {"metric": "accession_without_isolate_or_strain_n", "value": len(no_isolate)},
        {"metric": "unique_isolate_key_n", "value": len(audit)},
        {"metric": "isolate_with_at_least_2_segments_n", "value": two_plus_n},
        {"metric": "isolate_with_LMS_n", "value": lms_n},
        {"metric": "isolate_with_LMS_complete_coding_n", "value": lms_complete_coding_n},
        {"metric": "LMS_unique_high_confidence_n", "value": high_conf_n},
    ]
    for status, n in sorted(status_counts.items()):
        qc_rows.append({"metric": f"pairing_status::{status}", "value": n})

    out_qc = base / "qc" / "04_TSWV_segment_pairing_QC.tsv"
    write_tsv(out_qc, qc_rows, ["metric", "value"])

    print("\n===== SEGMENT PAIRING AUDIT =====")
    for r in qc_rows[:8]:
        print(f"{r['metric']}: {r['value']}")
    print("Pairing status counts:")
    for status, n in sorted(status_counts.items()):
        print(f"  {status}: {n}")

    # This is only a data-completeness decision, NOT evidence of reassortment.
    if high_conf_n >= 30:
        decision = "FEASIBLE: enough high-confidence L/M/S isolate pairs for a formal segment-concordance/reassortment analysis."
    elif lms_n >= 30:
        decision = "POTENTIALLY FEASIBLE: enough L/M/S isolate labels, but duplicate/metadata conflicts require resolution before phylogenetic analysis."
    else:
        decision = "NOT YET ROBUST: too few reliable L/M/S isolate pairs for a strong reassortment analysis."
    print("Decision:", decision)
    print(f"Segment pairing audit: {out_audit}")
    print(f"Segment pairing QC   : {out_qc}")
    return audit, qc_rows


def main():
    root = find_analysis_root()
    base = root / "06_TSWV_evolution"
    print("Analysis root:", root)
    print("TSWV base    :", base)

    diagnose_gard(base)
    segment_pairing_audit(base)

    print("\nIMPORTANT:")
    print("The segment-pairing result only tells whether a segment-tree comparison is technically supportable.")
    print("It does NOT itself demonstrate reassortment; establishing that requires testing topological discordance between the L/M/S segment phylogenies.")


if __name__ == "__main__":
    main()
