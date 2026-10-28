#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, os, re, sys
from pathlib import Path

GENES = ["RdRp", "GnGc", "NSm", "N", "NSs"]
VALID_STATUSES = {
    "FINAL_GARD",
    "PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1",
    "PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1",
}


def root_from_env_or_file() -> Path:
    if os.environ.get("ANALYSIS_ROOT"):
        return Path(os.environ["ANALYSIS_ROOT"]).resolve()
    return Path(__file__).resolve().parents[2]


def read_fasta(path: Path):
    seqs = []
    name = None
    chunks = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs.append((name, "".join(chunks)))
    if not seqs:
        raise ValueError(f"No sequences in {path}")
    lengths = {len(s) for _, s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Non-rectangular alignment {path}: lengths={sorted(lengths)[:10]}")
    return seqs, next(iter(lengths))


def numeric_values(obj):
    if isinstance(obj, (int, float)) and not isinstance(obj, bool):
        yield float(obj)
    elif isinstance(obj, list):
        for x in obj:
            yield from numeric_values(x)
    elif isinstance(obj, dict):
        for k in sorted(obj, key=lambda x: (str(x))):
            yield from numeric_values(obj[k])


def partition_rows_from_json(js):
    bd = js.get("breakpointData")
    if not isinstance(bd, dict) or not bd:
        return []
    def keyfun(k):
        try:
            return (0, int(k))
        except Exception:
            return (1, str(k))
    out = []
    for k in sorted(bd, key=keyfun):
        v = bd[k]
        if not isinstance(v, dict):
            continue
        nums = list(numeric_values(v.get("bps", [])))
        start = int(round(nums[0])) if len(nums) >= 1 else None
        end = int(round(nums[1])) if len(nums) >= 2 else None
        tree = v.get("tree", "")
        out.append((str(k), start, end, tree))
    return out


def parse_stdout(log_path: Path):
    text = log_path.read_text(errors="replace") if log_path.exists() else ""
    baseline = None
    # Baseline line: * Log(L) = ..., AIC-c = 16972.32 (...)
    m = re.search(r"AIC-c\s*=\s*([0-9eE+\-.]+)", text)
    if m:
        baseline = float(m.group(1))
    done_single = "Done with single breakpoint analysis." in text
    single_loc = None
    m = re.search(r"Best single break point location:\s*([0-9]+)", text)
    if m:
        single_loc = int(m.group(1))
    single_caic = None
    # Prefer c-AIC immediately after the 'Best single' report.
    m = re.search(
        r"Best single break point location:\s*[0-9]+\s*\n\s*c-AIC\s*=\s*([0-9eE+\-.]+)",
        text,
        flags=re.I,
    )
    if m:
        single_caic = float(m.group(1))
    return {
        "done_single": done_single,
        "baseline_log": baseline,
        "single_location_log": single_loc,
        "single_caic_log": single_caic,
    }


def analyze_gene(base: Path, gene: str):
    gd = base / "pipeline_results" / "04_GARD" / gene
    fasta = gd / f"{gene}.GARD_input.fna"
    js_path = gd / f"{gene}.GARD.json"
    nex_path = gd / f"{gene}.best-gard.nex"
    stdout = base / "logs" / f"{gene}.GARD.stdout.md"

    row = {
        "gene": gene,
        "status": "FAIL_MISSING_GARD_INPUT",
        "reason": "",
        "sequence_n": "",
        "alignment_nt": "",
        "min_partition_nt": "",
        "max_feasible_breakpoint_n": "",
        "masterList_present": "",
        "masterList_n": "",
        "final_breakpoint_n": "",
        "raw_breakpoint_positions_nt": "",
        "baseline_cAIC": "",
        "best_single_breakpoint_cAIC": "",
        "bestModelAICc": "",
        "singleTreeAICc": "",
        "selection_input_nexus": str(nex_path),
        "selection_input_status": "MISSING",
        "json_path": str(js_path),
        "stdout_path": str(stdout),
    }

    if not fasta.exists() or fasta.stat().st_size == 0:
        row["reason"] = "Missing GARD input FASTA"
        return row
    seqs, aln_len = read_fasta(fasta)
    nseq = len(seqs)
    min_part = 2 * nseq - 3
    max_feasible = max(0, aln_len // min_part - 1) if min_part > 0 else 0
    row.update({
        "sequence_n": nseq,
        "alignment_nt": aln_len,
        "min_partition_nt": min_part,
        "max_feasible_breakpoint_n": max_feasible,
    })

    if not js_path.exists() or js_path.stat().st_size == 0:
        row["status"] = "FAIL_MISSING_GARD_JSON"
        row["reason"] = "Missing/empty GARD JSON"
        return row
    try:
        with js_path.open() as f:
            js = json.load(f)
    except Exception as e:
        row["status"] = "FAIL_INVALID_GARD_JSON"
        row["reason"] = repr(e)
        return row

    log = parse_stdout(stdout)
    baseline = js.get("baselineScore", log["baseline_log"])
    if baseline is None:
        baseline = log["baseline_log"]
    single_caic = log["single_caic_log"]
    row["baseline_cAIC"] = baseline if baseline is not None else ""
    row["best_single_breakpoint_cAIC"] = single_caic if single_caic is not None else ""
    row["bestModelAICc"] = js.get("bestModelAICc", "")
    row["singleTreeAICc"] = js.get("singleTreeAICc", "")

    master_present = isinstance(js.get("masterList"), dict)
    row["masterList_present"] = str(master_present).upper()
    row["masterList_n"] = len(js.get("masterList", {})) if master_present else 0

    parts = partition_rows_from_json(js)
    bp_positions = [p[2] for p in parts[:-1] if p[2] is not None]
    row["raw_breakpoint_positions_nt"] = ",".join(map(str, bp_positions))
    inferred_bp_n = max(0, len(parts) - 1) if parts else ""

    if nex_path.exists() and nex_path.stat().st_size > 0:
        row["selection_input_status"] = "PRESENT"

    if not master_present:
        row["status"] = "FINAL_GARD"
        row["final_breakpoint_n"] = inferred_bp_n
        row["reason"] = "GARD completed normally; masterList absent from final JSON"
        if row["selection_input_status"] != "PRESENT":
            row["status"] = "FAIL_MISSING_SELECTION_SNAPSHOT"
            row["reason"] = "Final GARD JSON present but best-gard.nex missing/empty"
        return row

    # GARD v0.2 writes a resumable checkpoint after the exhaustive single-breakpoint step.
    # A partial JSON is scientifically complete only when >1 breakpoint is mathematically
    # impossible under GARD's own minPartitionSize heuristic, and the exhaustive single-BP
    # search itself has completed.
    if max_feasible <= 1 and log["done_single"] and baseline is not None and single_caic is not None:
        if single_caic < float(baseline):
            row["status"] = "PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1"
            row["final_breakpoint_n"] = 1
            row["reason"] = (
                "All feasible breakpoint models were exhausted; one breakpoint improves cAIC, "
                "while >=2 breakpoints are impossible under GARD minPartitionSize"
            )
        else:
            row["status"] = "PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1"
            row["final_breakpoint_n"] = 0
            row["reason"] = (
                "All feasible breakpoint models were exhausted; the best one-breakpoint model "
                "does not improve the no-breakpoint baseline, and >=2 breakpoints are impossible "
                "under GARD minPartitionSize"
            )
            # The checkpoint breakpointData corresponds to the baseline single-tree model.
            row["raw_breakpoint_positions_nt"] = ""
        if row["selection_input_status"] != "PRESENT":
            row["status"] = "FAIL_MISSING_SELECTION_SNAPSHOT"
            row["reason"] = "Feasible GARD search complete but best-gard.nex missing/empty"
        return row

    row["status"] = "FAIL_PARTIAL_GARD_REQUIRES_REVIEW"
    row["final_breakpoint_n"] = inferred_bp_n
    row["reason"] = (
        "GARD JSON still contains masterList and cannot be closed by the max-feasible-breakpoint rule"
    )
    return row


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", choices=GENES)
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args()

    root = root_from_env_or_file()
    base = root / "06_TSWV_evolution"
    genes = [args.gene] if args.gene else GENES
    rows = [analyze_gene(base, g) for g in genes]

    for r in rows:
        print(
            r["gene"], r["status"],
            f"n={r['sequence_n']}", f"nt={r['alignment_nt']}",
            f"minPart={r['min_partition_nt']}",
            f"maxBP={r['max_feasible_breakpoint_n']}",
            f"finalBP={r['final_breakpoint_n']}",
            f"baseline={r['baseline_cAIC']}",
            f"best1={r['best_single_breakpoint_cAIC']}"
        )

    bad = [r for r in rows if r["status"] not in VALID_STATUSES]

    if not args.no_write and not args.gene:
        fields = list(rows[0].keys())
        out = base / "tables" / "16_TSWV_GARD_final_status.tsv"
        write_tsv(out, rows, fields)
        qc_rows = [
            {"metric": "status", "value": "PASS" if not bad else "FAIL"},
            {"metric": "gene_n", "value": len(rows)},
            {"metric": "valid_gene_n", "value": len(rows) - len(bad)},
            {"metric": "FINAL_GARD_n", "value": sum(r["status"] == "FINAL_GARD" for r in rows)},
            {"metric": "PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1_n", "value": sum(r["status"] == "PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1" for r in rows)},
            {"metric": "PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1_n", "value": sum(r["status"] == "PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1" for r in rows)},
        ]
        write_tsv(base / "qc" / "07_TSWV_GARD_final_status_QC.tsv", qc_rows, ["metric", "value"])
        print("Status table:", out)

    if bad:
        print("GARD status gate FAIL:", ", ".join(f"{r['gene']}={r['status']}" for r in bad), file=sys.stderr)
        return 2
    print("Fig7 07O GARD final-status gate PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
