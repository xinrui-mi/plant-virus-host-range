#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict, Counter
import csv, math, argparse

IE_CUTOFF = 1e-5

def read_fasta(path):
    d = {}
    h = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    d[h] = "".join(seq)
                h = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if h is not None:
        d[h] = "".join(seq)
    return d

def read_tsv(path):
    with open(path) as f:
        yield from csv.DictReader(f, delimiter="\t")

def parse_domtbl(path):
    hits = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 22:
                continue
            r = {
                "model": p[0],
                "model_len": int(p[2]),
                "query": p[3],
                "query_len": int(p[5]),
                "full_e": float(p[6]),
                "i_e": float(p[12]),
                "hmm_from": int(p[15]),
                "hmm_to": int(p[16]),
                "ali_from": int(p[17]),
                "ali_to": int(p[18]),
            }
            hits[r["query"]].append(r)
    return hits

def all_occurrences(seq, motif, start_1, end_1):
    # 1-based inclusive search window; return 1-based motif starts.
    start = max(1, start_1) - 1
    end = min(len(seq), end_1)
    out = []
    i = seq.find(motif, start, end)
    while i != -1:
        out.append(i + 1)
        i = seq.find(motif, i + 1, end)
    return out

def is_protease_like(model):
    m = model.lower()
    return (
        model.startswith("Peptidase_")
        or "protease" in m
        or "proteas" in m
        or model.startswith("RVP")
        or model.startswith("Trypsin")
    )

def best_rdrp(hs):
    z = [h for h in hs if h["i_e"] <= IE_CUTOFF and h["model"].startswith("RdRP_")]
    return min(z, key=lambda h: h["i_e"]) if z else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root).resolve()
    base = root / "three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    kroot = base / "pipeline_results/03K_family_specific_markers"

    faa = kroot / "01_orf_recovery/Secoviridae.candidate_orfs.faa"
    meta = kroot / "01_orf_recovery/Secoviridae.candidate_orfs.tsv"
    dom = kroot / "02_pfam/Secoviridae.Pfam.domtblout"
    seqmaster = root / "final_three_source_pipeline_v1/results/20_three_source_sequence_master_FINAL.tsv"

    outdir = kroot / "04_Secoviridae_marker_audit"
    outdir.mkdir(parents=True, exist_ok=True)

    seqs = read_fasta(faa)
    meta_rows = {r["orf_id"]: r for r in read_tsv(meta)}
    hits = parse_domtbl(dom)
    sm = {r["sequence_id"]: r for r in read_tsv(seqmaster)}

    detail = []
    protease_model_counter = Counter()
    protease_by_genus = Counter()

    for oid, m in meta_rows.items():
        hs = hits.get(oid, [])
        pol = best_rdrp(hs)
        if not pol:
            continue

        seq = seqs.get(oid, "")
        sid = m["sequence_id"]
        genus = (sm.get(sid, {}).get("virus_genus") or "").strip() or "NA"
        species = (sm.get(sid, {}).get("virus_species") or "").strip() or "NA"

        # GDD should be in/near the polymerase alignment, but use a generous window.
        gdd = all_occurrences(
            seq, "GDD",
            max(1, pol["ali_from"] - 120),
            min(len(seq), pol["ali_to"] + 180)
        )

        # CG motif of the 3C-like proteinase must be upstream of Pol.
        # Audit a broad window rather than assuming a particular Pfam peptidase model.
        cg = all_occurrences(
            seq, "CG",
            max(1, pol["ali_from"] - 1800),
            max(1, pol["ali_from"] - 1)
        )

        proteases = [
            h for h in hs
            if h["i_e"] <= IE_CUTOFF
            and is_protease_like(h["model"])
            and h["ali_from"] < pol["ali_from"]
        ]
        proteases.sort(key=lambda h: h["i_e"])

        for h in proteases:
            protease_model_counter[h["model"]] += 1
            protease_by_genus[(genus, h["model"])] += 1

        # CG positions close to any significant upstream protease-like domain.
        cg_near_pro = []
        for c in cg:
            for h in proteases:
                if h["ali_from"] - 120 <= c <= h["ali_to"] + 120:
                    cg_near_pro.append(c)
                    break

        # Evaluate possible CG->GDD regions without freezing a rule yet.
        candidate_pairs = []
        for c in cg:
            for g in gdd:
                end = g + 2
                L = end - c + 1
                if 200 <= L <= 1400:
                    candidate_pairs.append((c, g, L))

        near_pairs = []
        for c in cg_near_pro:
            for g in gdd:
                end = g + 2
                L = end - c + 1
                if 200 <= L <= 1400:
                    near_pairs.append((c, g, L))

        best_pair = None
        pair_source = ""
        if near_pairs:
            # Typical Pro-Pol is several hundred aa; choose a central plausible length
            # only for audit reporting, not final extraction.
            best_pair = min(near_pairs, key=lambda x: abs(x[2] - 600))
            pair_source = "CG_near_significant_protease"
        elif candidate_pairs:
            best_pair = min(candidate_pairs, key=lambda x: abs(x[2] - 600))
            pair_source = "CG_broad_upstream_window"

        detail.append({
            "sequence_id": sid,
            "virus_genus": genus,
            "virus_species": species,
            "orf_id": oid,
            "orf_aa_length": len(seq),
            "best_RdRP_model": pol["model"],
            "best_RdRP_iE": f"{pol['i_e']:.6g}",
            "RdRP_ali_from": pol["ali_from"],
            "RdRP_ali_to": pol["ali_to"],
            "significant_upstream_protease_models": ";".join(sorted({h["model"] for h in proteases})),
            "significant_upstream_protease_n": len(proteases),
            "CG_any_upstream_n": len(cg),
            "CG_near_protease_n": len(set(cg_near_pro)),
            "GDD_near_RdRP_n": len(gdd),
            "plausible_CG_GDD_pair": "YES" if best_pair else "NO",
            "candidate_pair_source": pair_source,
            "candidate_CG_pos": best_pair[0] if best_pair else "",
            "candidate_GDD_pos": best_pair[1] if best_pair else "",
            "candidate_ProPol_length_aa": best_pair[2] if best_pair else "",
        })

    # Keep one best audit ORF per sequence for genus/species coverage summaries.
    by_sid = defaultdict(list)
    for r in detail:
        by_sid[r["sequence_id"]].append(r)

    selected = []
    for sid, rows in by_sid.items():
        # Prefer plausible motif pair, then protease-supported pair, then best RdRP E.
        def key(r):
            plausible = 1 if r["plausible_CG_GDD_pair"] == "YES" else 0
            pro = 1 if r["CG_near_protease_n"] > 0 else 0
            ie = float(r["best_RdRP_iE"])
            return (plausible, pro, -math.log10(max(ie, 1e-300)))
        selected.append(max(rows, key=key))

    def summarize(rows, group_field):
        g = defaultdict(list)
        for r in rows:
            g[r[group_field]].append(r)
        out = []
        for group, rr in sorted(g.items()):
            sids = {r["sequence_id"] for r in rr}
            spp = {r["virus_species"] for r in rr if r["virus_species"] != "NA"}
            pro = {r["sequence_id"] for r in rr if r["significant_upstream_protease_n"] > 0}
            cg = {r["sequence_id"] for r in rr if r["CG_any_upstream_n"] > 0}
            cgpro = {r["sequence_id"] for r in rr if r["CG_near_protease_n"] > 0}
            gdd = {r["sequence_id"] for r in rr if r["GDD_near_RdRP_n"] > 0}
            pair = {r["sequence_id"] for r in rr if r["plausible_CG_GDD_pair"] == "YES"}
            out.append({
                group_field: group,
                "sequence_with_RdRP_n": len(sids),
                "species_with_RdRP_n": len(spp),
                "sequence_with_upstream_protease_HMM_n": len(pro),
                "sequence_with_any_upstream_CG_n": len(cg),
                "sequence_with_CG_near_protease_n": len(cgpro),
                "sequence_with_GDD_near_RdRP_n": len(gdd),
                "sequence_with_plausible_CG_GDD_pair_n": len(pair),
                "plausible_pair_pct_of_RdRP": f"{100*len(pair)/len(sids):.2f}" if sids else "NA",
            })
        return out

    # Write detail
    detail_path = outdir / "01_Secoviridae_RdRP_ORF_motif_audit.tsv"
    fields = list(detail[0].keys()) if detail else []
    with open(detail_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(detail)

    selected_path = outdir / "02_Secoviridae_best_ORF_per_sequence.tsv"
    with open(selected_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(selected)

    genus_rows = summarize(selected, "virus_genus")
    genus_path = outdir / "03_Secoviridae_genus_audit.tsv"
    with open(genus_path, "w", newline="") as f:
        gf = list(genus_rows[0].keys()) if genus_rows else []
        w = csv.DictWriter(f, fieldnames=gf, delimiter="\t")
        w.writeheader()
        w.writerows(genus_rows)

    model_path = outdir / "04_upstream_protease_Pfam_models.tsv"
    with open(model_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Pfam_model", "significant_hit_ORF_n"])
        for model, n in protease_model_counter.most_common():
            w.writerow([model, n])

    gbm_path = outdir / "05_genus_by_upstream_protease_Pfam.tsv"
    with open(gbm_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["virus_genus", "Pfam_model", "significant_hit_ORF_n"])
        for (genus, model), n in sorted(protease_by_genus.items()):
            w.writerow([genus, model, n])

    print("03K-B2 Secoviridae audit PASS")
    print(f"Output: {outdir}")
    print(f"ORFs_with_significant_RdRP={len(detail)}")
    print(f"sequences_with_significant_RdRP={len({r['sequence_id'] for r in selected})}")
    print(f"sequences_with_plausible_CG_GDD_pair={sum(r['plausible_CG_GDD_pair']=='YES' for r in selected)}")
    print("genus_summary:")
    for r in genus_rows:
        print(
            r["virus_genus"],
            f"RdRP_seq={r['sequence_with_RdRP_n']}",
            f"proteaseHMM={r['sequence_with_upstream_protease_HMM_n']}",
            f"CG={r['sequence_with_any_upstream_CG_n']}",
            f"GDD={r['sequence_with_GDD_near_RdRP_n']}",
            f"pair={r['sequence_with_plausible_CG_GDD_pair_n']}",
            f"pair_pct={r['plausible_pair_pct_of_RdRP']}",
        )

if __name__ == "__main__":
    main()
