#!/usr/bin/env python3
from __future__ import annotations
import csv, os, re
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ["ANALYSIS_ROOT"])
WORK = ROOT / "02_host_relatedness_and_phylogeny" / "work" / "marker_inputs"
EVALUE = float(os.environ.get("FIG3_MARKER_EVALUE", "1e-5"))
HMM_COV = float(os.environ.get("FIG3_MARKER_HMM_COV", "0.50"))
FLANK = int(os.environ.get("FIG3_DNA_MARKER_FLANK", "15"))


def read_fasta(path: Path):
    d = {}; h = None; buf = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None: d[h] = "".join(buf)
                h = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
        if h is not None: d[h] = "".join(buf)
    return d


def read_candidates(path: Path):
    d = {}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[r["protein_id"]] = r
    return d


def parse_domtbl(path: Path):
    out = []
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.strip().split()
            if len(p) < 23: continue
            try:
                h = {
                    "model_name": p[0], "model_acc": p[1], "hmm_len": int(p[2]),
                    "protein_id": p[3], "protein_len": int(p[5]),
                    "full_e": float(p[6]), "full_score": float(p[7]),
                    "i_evalue": float(p[12]), "domain_score": float(p[13]),
                    "hmm_from": int(p[15]), "hmm_to": int(p[16]),
                    "ali_from": int(p[17]), "ali_to": int(p[18]),
                    "env_from": int(p[19]), "env_to": int(p[20]),
                    "description": " ".join(p[22:]),
                }
            except ValueError:
                continue
            h["hmm_cov"] = (h["hmm_to"] - h["hmm_from"] + 1) / max(1, h["hmm_len"])
            if h["i_evalue"] <= EVALUE and h["hmm_cov"] >= HMM_COV:
                out.append(h)
    return out


def text(h):
    return f"{h['model_name']} {h['model_acc']} {h['description']}".lower()


def is_gem_catalytic(h):
    s = text(h)
    return (
        h["model_name"].lower() == "gemini_al1" or
        "geminivirus rep catalytic" in s or
        ("geminivirus" in s and "rep" in s and "catalytic" in s)
    )


def is_gem_extension(h):
    s = text(h)
    n = h["model_name"].lower()
    return n.startswith("gemini_al1") or ("geminivirus" in s and "rep" in s)


def is_rt(h):
    s = text(h)
    n = h["model_name"].lower()
    if "rnase h" in s or "ribonuclease h" in s or "rnase_h" in n:
        return False
    return n.startswith("rvt") or "reverse transcriptase" in s or "rna-dependent dna polymerase" in s


def is_rh(h):
    s = text(h)
    n = h["model_name"].lower()
    return "rnase h" in s or "ribonuclease h" in s or "rnase_h" in n


def write_fasta(path: Path, rows):
    with path.open("w") as out:
        for r in rows:
            out.write(f">{r['sequence_id']}\n{r['marker_seq']}\n")


def write_tsv(path: Path, rows, fields):
    with path.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

# -------- Geminiviridae --------
gseq = read_fasta(WORK / "Geminiviridae.Rep_candidates.faa")
gmap = read_candidates(WORK / "Geminiviridae.Rep_candidates.tsv")
ghits = parse_domtbl(WORK / "Geminiviridae.Pfam.domtbl")
by_pid = defaultdict(list)
for h in ghits:
    by_pid[h["protein_id"]].append(h)

g_candidates = []
for pid, hs in by_pid.items():
    if pid not in gseq or pid not in gmap: continue
    catalytic = [h for h in hs if is_gem_catalytic(h)]
    if not catalytic: continue
    selected = [h for h in hs if is_gem_extension(h)]
    if not selected: selected = catalytic
    start = max(1, min(h["env_from"] for h in selected) - FLANK)
    end = min(len(gseq[pid]), max(h["env_to"] for h in selected) + FLANK)
    marker = gseq[pid][start-1:end]
    if len(marker) < 100: continue
    sid = gmap[pid]["sequence_id"]
    g_candidates.append({
        "sequence_id": sid, "protein_id": pid, "marker_start": start, "marker_end": end,
        "marker_length": len(marker), "best_catalytic_iE": min(h["i_evalue"] for h in catalytic),
        "best_domain_score": max(h["domain_score"] for h in catalytic),
        "included_models": ";".join(sorted({h["model_name"] for h in selected})),
        "included_accessions": ";".join(sorted({h["model_acc"] for h in selected})),
        "marker_seq": marker,
    })

g_best = {}
for r in g_candidates:
    sid = r["sequence_id"]
    key = (r["best_domain_score"], r["marker_length"])
    if sid not in g_best or key > g_best[sid][0]: g_best[sid] = (key, r)
g_rows = [g_best[sid][1] for sid in sorted(g_best)]
write_fasta(WORK / "Geminiviridae.Rep_core.raw.faa", g_rows)
write_tsv(WORK / "03E_Geminiviridae_marker_detail.tsv", g_rows,
          ["sequence_id","protein_id","marker_start","marker_end","marker_length","best_catalytic_iE","best_domain_score","included_models","included_accessions"])

# -------- Caulimoviridae --------
cseq = read_fasta(WORK / "Caulimoviridae.Pol_candidates.faa")
cmap = read_candidates(WORK / "Caulimoviridae.Pol_candidates.tsv")
chits = parse_domtbl(WORK / "Caulimoviridae.Pfam.domtbl")
by_pid = defaultdict(list)
for h in chits:
    by_pid[h["protein_id"]].append(h)

c_candidates = []
same_cds_both = 0
ordered = 0
for pid, hs in by_pid.items():
    if pid not in cseq or pid not in cmap: continue
    rts = [h for h in hs if is_rt(h)]
    rhs = [h for h in hs if is_rh(h)]
    if not rts or not rhs: continue
    same_cds_both += 1
    best_pair = None
    for rt in rts:
        for rh in rhs:
            if rt["env_from"] >= rh["env_from"]: continue
            score = rt["domain_score"] + rh["domain_score"]
            if best_pair is None or score > best_pair[0]: best_pair = (score, rt, rh)
    if best_pair is None: continue
    ordered += 1
    score, rt, rh = best_pair
    start = max(1, min(rt["env_from"], rh["env_from"]) - FLANK)
    end = min(len(cseq[pid]), max(rt["env_to"], rh["env_to"]) + FLANK)
    marker = cseq[pid][start-1:end]
    if len(marker) < 120: continue
    sid = cmap[pid]["sequence_id"]
    c_candidates.append({
        "sequence_id": sid, "protein_id": pid,
        "marker_start": start, "marker_end": end, "marker_length": len(marker),
        "RT_model": rt["model_name"], "RT_acc": rt["model_acc"], "RT_iE": rt["i_evalue"],
        "RH_model": rh["model_name"], "RH_acc": rh["model_acc"], "RH_iE": rh["i_evalue"],
        "combined_domain_score": score, "marker_seq": marker,
    })

c_best = {}
for r in c_candidates:
    sid = r["sequence_id"]
    key = (r["combined_domain_score"], r["marker_length"])
    if sid not in c_best or key > c_best[sid][0]: c_best[sid] = (key, r)
c_rows = [c_best[sid][1] for sid in sorted(c_best)]
write_fasta(WORK / "Caulimoviridae.RT_RNaseH.raw.faa", c_rows)
write_tsv(WORK / "03E_Caulimoviridae_marker_detail.tsv", c_rows,
          ["sequence_id","protein_id","marker_start","marker_end","marker_length","RT_model","RT_acc","RT_iE","RH_model","RH_acc","RH_iE","combined_domain_score"])

with (WORK / "03E_DNA_marker_selection_QC.tsv").open("w", newline="") as out:
    fields = ["family","candidate_protein_n","qualifying_pfam_domain_hit_n","same_CDS_RT_RH_n","ordered_RT_RH_CDS_n","marker_sequence_n"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); w.writeheader()
    w.writerow({"family":"Geminiviridae","candidate_protein_n":len(gseq),"qualifying_pfam_domain_hit_n":len(ghits),"same_CDS_RT_RH_n":"","ordered_RT_RH_CDS_n":"","marker_sequence_n":len(g_rows)})
    w.writerow({"family":"Caulimoviridae","candidate_protein_n":len(cseq),"qualifying_pfam_domain_hit_n":len(chits),"same_CDS_RT_RH_n":same_cds_both,"ordered_RT_RH_CDS_n":ordered,"marker_sequence_n":len(c_rows)})

print("DNA full-Pfam marker selection PASS", len(g_rows), len(c_rows))
