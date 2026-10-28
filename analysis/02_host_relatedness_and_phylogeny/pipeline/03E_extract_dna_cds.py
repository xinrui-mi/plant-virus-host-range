#!/usr/bin/env python3
from __future__ import annotations
import csv, os, re
from pathlib import Path
from Bio import SeqIO

ROOT = Path(os.environ["ANALYSIS_ROOT"])
GB = Path(os.environ["DNA_GENBANK"])
WORK = ROOT / "02_host_relatedness_and_phylogeny" / "work" / "marker_inputs"
WORK.mkdir(parents=True, exist_ok=True)

FAMILIES = ["Geminiviridae", "Caulimoviridae"]


def base_acc(x: str) -> str:
    return re.sub(r"\.[0-9]+$", "", (x or "").strip())


def load_targets(fam: str):
    vals = [x.strip() for x in (WORK / f"{fam}.target_sequence_ids.txt").read_text().splitlines() if x.strip()]
    exact = set(vals)
    base = {base_acc(x) for x in vals}
    return vals, exact, base


def clean_translation(x: str) -> str:
    x = re.sub(r"\s+", "", x or "").upper().rstrip("*")
    return x


def translate_feature(feat, rec):
    q = feat.qualifiers
    tr = clean_translation((q.get("translation") or [""])[0])
    if tr:
        return tr, "qualifier"
    try:
        nt = feat.extract(rec.seq)
        codon_start = int((q.get("codon_start") or ["1"])[0]) - 1
        nt = nt[codon_start:]
        usable = (len(nt) // 3) * 3
        if usable < 3:
            return "", "fallback_too_short"
        nt = nt[:usable]
        table = int((q.get("transl_table") or ["1"])[0])
        aa = clean_translation(str(nt.translate(table=table, to_stop=False)))
        if "*" in aa:
            return "", "fallback_internal_stop"
        return aa, "fallback"
    except Exception:
        return "", "fallback_failed"


def is_partial(feat):
    s = str(feat.location)
    return int("<" in s), int(">" in s)


def rec_keys(rec):
    vals = {rec.id, base_acc(rec.id)}
    for x in rec.annotations.get("accessions", []) or []:
        vals.add(str(x)); vals.add(base_acc(str(x)))
    return {x for x in vals if x}


qc_rows = []
for fam in FAMILIES:
    target_raw, target_exact, target_base = load_targets(fam)
    out_faa = WORK / f"{fam}.all_CDS.faa"
    out_map = WORK / f"{fam}.all_CDS_map.tsv"
    missing_out = WORK / f"{fam}.missing_from_DNA_GenBank.txt"

    found_base = set()
    cds_n = 0
    fallback_n = 0
    partial_n = 0
    internal_stop_skipped = 0
    translation_failed = 0

    with out_faa.open("w") as faa, out_map.open("w", newline="") as mf:
        fieldnames = [
            "protein_id", "sequence_id", "record_id", "gene", "product", "note",
            "aa_length", "translation_source", "partial_5p", "partial_3p"
        ]
        w = csv.DictWriter(mf, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for rec in SeqIO.parse(str(GB), "genbank"):
            keys = rec_keys(rec)
            hit_base = next((base_acc(x) for x in keys if base_acc(x) in target_base), None)
            if hit_base is None:
                continue
            # Prefer the released sequence ID spelling if it exists; otherwise record.id.
            sid = next((x for x in target_raw if base_acc(x) == hit_base), rec.id)
            found_base.add(hit_base)
            idx = 0
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                idx += 1
                aa, source = translate_feature(feat, rec)
                if source == "fallback_internal_stop":
                    internal_stop_skipped += 1
                    continue
                if not aa:
                    translation_failed += 1
                    continue
                if len(aa) < 40:
                    continue
                p5, p3 = is_partial(feat)
                partial_n += int(p5 or p3)
                fallback_n += int(source == "fallback")
                pid = f"{sid}__CDS{idx:03d}"
                q = feat.qualifiers
                gene = ";".join(q.get("gene", []))
                product = ";".join(q.get("product", []))
                note = " | ".join(q.get("note", []))
                faa.write(f">{pid}\n{aa}\n")
                w.writerow({
                    "protein_id": pid, "sequence_id": sid, "record_id": rec.id,
                    "gene": gene, "product": product, "note": note,
                    "aa_length": len(aa), "translation_source": source,
                    "partial_5p": p5, "partial_3p": p3,
                })
                cds_n += 1

    missing = [x for x in target_raw if base_acc(x) not in found_base]
    missing_out.write_text("\n".join(sorted(missing)) + ("\n" if missing else ""))
    qc_rows.append({
        "family": fam,
        "target_sequence_n": len(target_raw),
        "GenBank_record_found_n": len(found_base),
        "GenBank_record_missing_n": len(missing),
        "CDS_extracted_n": cds_n,
        "fallback_translation_n": fallback_n,
        "partial_CDS_kept_n": partial_n,
        "internal_stop_skipped_n": internal_stop_skipped,
        "translation_failed_n": translation_failed,
        "record_recovery_fraction": len(found_base) / len(target_raw) if target_raw else 0,
    })
    print(
        fam, "targets", len(target_raw), "records_found", len(found_base), "CDS", cds_n,
        "missing", len(missing), "fallback_translation", fallback_n,
        "partial_CDS_kept", partial_n, "internal_stop_skipped", internal_stop_skipped
    )

with (WORK / "03E_DNA_GenBank_recovery_QC.tsv").open("w", newline="") as f:
    fields = list(qc_rows[0].keys())
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(qc_rows)
