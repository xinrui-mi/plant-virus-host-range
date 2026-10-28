#!/usr/bin/env python3
from __future__ import annotations
import csv, os, re
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ["ANALYSIS_ROOT"])
WORK = ROOT / "02_host_relatedness_and_phylogeny" / "work" / "marker_inputs"


def norm(x: str) -> str:
    return re.sub(r"\s+", " ", (x or "").strip().lower())


def read_fasta(path: Path):
    seqs = {}
    h = None
    buf = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(buf)
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if h is not None:
            seqs[h] = "".join(buf)
    return seqs


def read_map(fam: str):
    rows = []
    with (WORK / f"{fam}.all_CDS_map.tsv").open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["aa_length"] = int(r.get("aa_length", "0") or 0)
            rows.append(r)
    return rows


def gem_score(r):
    gene = norm(r.get("gene", ""))
    product = norm(r.get("product", ""))
    note = norm(r.get("note", ""))
    text = " ".join([gene, product, note])
    L = r["aa_length"]
    if L < 80:
        return -1000, "too_short"
    excludes = [
        "coat protein", "capsid", "movement protein", "nuclear shuttle",
        "transport protein", "silencing suppressor", "transcriptional activator",
        "transactivator", "replication enhancer", "pathogenicity", "satellite",
        "beta c1", "betac1", "pre-coat", "precoat"
    ]
    if any(k in text for k in excludes):
        return -500, "non_rep_annotation"
    score = 0
    reasons = []
    if re.search(r"replication[- ]associated protein|replication initiator|rolling[- ]circle replication", text):
        score += 120; reasons.append("replication_associated")
    if re.search(r"\brep protein\b|\brepa\b", text):
        score += 100; reasons.append("rep_name")
    if gene in {"rep", "repa", "al1", "ac1"}:
        score += 140; reasons.append("rep_gene")
    if gene == "c1":
        score += 90 if re.search(r"rep|replication|initiator|rolling", text) else 25
        reasons.append("c1")
    if "replicase" in text:
        score += 45; reasons.append("replicase")
    if "replication" in text:
        score += 20
    if "helicase" in text:
        score += 15
    if 180 <= L <= 500:
        score += 20
    elif 120 <= L <= 700:
        score += 8
    return score, "+".join(reasons) if reasons else "no_rep_keyword"


def caul_score(r):
    gene = norm(r.get("gene", ""))
    product = norm(r.get("product", ""))
    note = norm(r.get("note", ""))
    text = " ".join([gene, product, note])
    L = r["aa_length"]
    if L < 180:
        return -1000, "too_short"
    excludes = ["coat protein", "capsid", "movement protein", "transactivator", "virion-associated", "matrix protein"]
    if any(k in text for k in excludes):
        return -500, "non_polymerase_annotation"
    score = 0
    reasons = []
    if "reverse transcriptase" in text:
        score += 140; reasons.append("reverse_transcriptase")
    if "rnase h" in text or "ribonuclease h" in text:
        score += 100; reasons.append("rnase_h")
    if "polymerase" in text:
        score += 70; reasons.append("polymerase")
    if "polyprotein" in text:
        score += 50; reasons.append("polyprotein")
    if "replicase" in text:
        score += 45; reasons.append("replicase")
    if re.search(r"\borf ?v\b|open reading frame v", text):
        score += 60; reasons.append("orf_v")
    if "protease" in text or "retropepsin" in text or "aspartic" in text:
        score += 15
    if L >= 300:
        score += 10
    if L >= 500:
        score += 10
    return score, "+".join(reasons) if reasons else "no_polymerase_keyword"


def choose_candidates(fam: str, scorer, fallback_min_len: int):
    seqs = read_fasta(WORK / f"{fam}.all_CDS.faa")
    rows = read_map(fam)
    by_sid = defaultdict(list)
    for r in rows:
        if r["protein_id"] in seqs:
            score, reason = scorer(r)
            rr = dict(r)
            rr["annotation_score"] = score
            rr["candidate_reason"] = reason
            by_sid[r["sequence_id"]].append(rr)

    selected = []
    audit = []
    for sid in sorted(by_sid):
        grp = by_sid[sid]
        positive = [r for r in grp if r["annotation_score"] > 0]
        if positive:
            positive.sort(key=lambda r: (r["annotation_score"], r["aa_length"]), reverse=True)
            # Keep up to two annotation-supported candidates; Pfam will arbitrate the final marker.
            keep = positive[:2]
            mode = "annotation_supported"
        else:
            plausible = [r for r in grp if r["aa_length"] >= fallback_min_len]
            plausible.sort(key=lambda r: r["aa_length"], reverse=True)
            keep = plausible[:2]
            mode = "fallback_longest_CDS"
        keep_ids = {r["protein_id"] for r in keep}
        for rank, r in enumerate(keep, 1):
            rr = dict(r)
            rr["selection_mode"] = mode
            rr["candidate_rank"] = rank
            selected.append(rr)
        audit.append({
            "family": fam,
            "sequence_id": sid,
            "CDS_n": len(grp),
            "positive_annotation_candidate_n": len(positive),
            "selected_candidate_n": len(keep),
            "selection_mode": mode if keep else "no_candidate",
            "selected_protein_ids": ";".join(r["protein_id"] for r in keep),
        })

    faa = WORK / ("Geminiviridae.Rep_candidates.faa" if fam == "Geminiviridae" else "Caulimoviridae.Pol_candidates.faa")
    tsv = WORK / ("Geminiviridae.Rep_candidates.tsv" if fam == "Geminiviridae" else "Caulimoviridae.Pol_candidates.tsv")
    with faa.open("w") as out:
        for r in selected:
            seq = seqs[r["protein_id"]]
            out.write(f">{r['protein_id']}\n{seq}\n")
    fields = [
        "protein_id","sequence_id","record_id","gene","product","note","aa_length",
        "translation_source","partial_5p","partial_3p","annotation_score","candidate_reason",
        "selection_mode","candidate_rank"
    ]
    with tsv.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(selected)
    return selected, audit


gem, aq1 = choose_candidates("Geminiviridae", gem_score, 120)
caul, aq2 = choose_candidates("Caulimoviridae", caul_score, 250)
with (WORK / "03E_DNA_candidate_selection_audit.tsv").open("w", newline="") as out:
    fields = ["family","sequence_id","CDS_n","positive_annotation_candidate_n","selected_candidate_n","selection_mode","selected_protein_ids"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(aq1 + aq2)

with (WORK / "03E_DNA_candidate_selection_QC.tsv").open("w", newline="") as out:
    fields = ["family","sequence_n","candidate_protein_n","annotation_supported_sequence_n","fallback_sequence_n","no_candidate_sequence_n"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); w.writeheader()
    for fam, audit, selected in [("Geminiviridae", aq1, gem), ("Caulimoviridae", aq2, caul)]:
        w.writerow({
            "family": fam,
            "sequence_n": len(audit),
            "candidate_protein_n": len(selected),
            "annotation_supported_sequence_n": sum(r["selection_mode"] == "annotation_supported" for r in audit),
            "fallback_sequence_n": sum(r["selection_mode"] == "fallback_longest_CDS" for r in audit),
            "no_candidate_sequence_n": sum(r["selection_mode"] == "no_candidate" for r in audit),
        })
print("DNA candidate selection PASS", len(gem), len(caul))
