#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse, csv, math, re, statistics

IE_CUTOFF = 1e-5
MIN_LEN = 200
MAX_LEN = 1200

def read_fasta(path):
    d = {}
    h = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    d[h] = "".join(buf)
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if h is not None:
        d[h] = "".join(buf)
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
            hits[p[3]].append({
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
            })
    return hits

def best_rdrp(hs):
    z = [h for h in hs if h["i_e"] <= IE_CUTOFF and h["model"].startswith("RdRP_")]
    return min(z, key=lambda x: x["i_e"]) if z else None

def is_protease_like(model):
    m = model.lower()
    return (
        model.startswith("Peptidase_")
        or "protease" in m
        or "proteas" in m
        or model.startswith("RVP")
        or model.startswith("Trypsin")
    )

def find_all(seq, motif):
    out = []
    start = 0
    while True:
        i = seq.find(motif, start)
        if i < 0:
            break
        out.append(i + 1)  # 1-based motif start
        start = i + 1
    return out

def choose_gdd(seq, pol):
    # First prefer GDD inside the RdRP alignment itself.
    candidates = []
    for p in find_all(seq, "GDD"):
        if pol["ali_from"] <= p <= pol["ali_to"]:
            candidates.append(p)
    if not candidates:
        # Then use a modest window around the HMM alignment.
        lo = max(1, pol["ali_from"] - 80)
        hi = min(len(seq), pol["ali_to"] + 150)
        candidates = [p for p in find_all(seq, "GDD") if lo <= p <= hi]
    if not candidates:
        return None

    # The catalytic GDD tends to sit in the later portion of the Pol domain.
    target = pol["ali_from"] + 0.75 * (pol["ali_to"] - pol["ali_from"])
    return min(candidates, key=lambda p: abs(p - target))

def gx_cg_starts(seq):
    # Return the C position (1-based) of GxCG motifs.
    # Example GSCG: G S C G -> C is the third residue.
    return [m.start() + 3 for m in re.finditer(r"G.CG", seq)]

def plain_cg_starts(seq):
    return [p for p in find_all(seq, "CG")]

def plausible(c_pos, gdd_pos):
    end = gdd_pos + 2
    L = end - c_pos + 1
    return MIN_LEN <= L <= MAX_LEN

def significant_upstream_proteases(hs, pol):
    z = [
        h for h in hs
        if h["i_e"] <= IE_CUTOFF
        and is_protease_like(h["model"])
        and h["ali_from"] < pol["ali_from"]
    ]
    return sorted(z, key=lambda x: x["i_e"])

def near_protease(c_pos, proteases):
    for h in proteases:
        if h["ali_from"] - 120 <= c_pos <= h["ali_to"] + 120:
            return True
    return False

def q(vals, p):
    vals = sorted(vals)
    if not vals:
        return None
    if len(vals) == 1:
        return float(vals[0])
    x = (len(vals)-1)*p
    lo, hi = math.floor(x), math.ceil(x)
    if lo == hi:
        return float(vals[lo])
    return vals[lo]*(hi-x) + vals[hi]*(x-lo)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root).resolve()
    base = root / "three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
    kroot = base / "pipeline_results/03K_family_specific_markers"

    faa_path = kroot / "01_orf_recovery/Secoviridae.candidate_orfs.faa"
    meta_path = kroot / "01_orf_recovery/Secoviridae.candidate_orfs.tsv"
    dom_path = kroot / "02_pfam/Secoviridae.Pfam.domtblout"
    seqmaster_path = root / "final_three_source_pipeline_v1/results/20_three_source_sequence_master_FINAL.tsv"
    assoc_path = root / "final_three_source_pipeline_v1/results/20_three_source_association_master_FINAL.tsv"

    outdir = kroot / "05_Secoviridae_corrected"
    outdir.mkdir(parents=True, exist_ok=True)

    seqs = read_fasta(faa_path)
    meta = {r["orf_id"]: r for r in read_tsv(meta_path)}
    hits = parse_domtbl(dom_path)
    seqmaster = {r["sequence_id"]: r for r in read_tsv(seqmaster_path)}

    host_fam = defaultdict(set)
    host_taxid = defaultdict(set)
    with open(assoc_path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        if "host_family" not in rd.fieldnames:
            raise RuntimeError("association master lacks host_family")
        for r in rd:
            sid = r["sequence_id"]
            hf = (r.get("host_family") or "").strip()
            ht = (r.get("host_taxid") or "").strip()
            if hf and hf.lower() not in {"na","nan","none","unclassified"}:
                host_fam[sid].add(hf)
            if ht and ht.lower() not in {"na","nan","none"}:
                host_taxid[sid].add(ht)

    # ---------- Pass 1: identify RdRP ORFs and high-confidence GxCG -> GDD pairs ----------
    per_orf = []
    by_sid = defaultdict(list)

    for oid, m in meta.items():
        seq = seqs.get(oid)
        if not seq:
            continue
        hs = hits.get(oid, [])
        pol = best_rdrp(hs)
        if not pol:
            continue
        gdd = choose_gdd(seq, pol)
        if not gdd:
            continue

        sid = m["sequence_id"]
        sm = seqmaster.get(sid, {})
        genus = (sm.get("virus_genus") or "").strip() or "NA"
        species = (sm.get("virus_species") or "").strip() or "NA"

        proteases = significant_upstream_proteases(hs, pol)

        gx = [c for c in gx_cg_starts(seq) if c < gdd and plausible(c, gdd)]
        gx_pro = [c for c in gx if near_protease(c, proteases)]

        # High-confidence choice:
        # 1) GxCG near a significant protease-like HMM hit.
        # 2) Otherwise the closest upstream GxCG compatible with Pro-Pol length.
        high_c = None
        high_method = None
        if gx_pro:
            high_c = max(gx_pro)
            high_method = "GxCG_near_protease_HMM"
        elif gx:
            high_c = max(gx)
            high_method = "GxCG_motif"

        rec = {
            "sequence_id": sid,
            "virus_genus": genus,
            "virus_species": species,
            "orf_id": oid,
            "translation_table": m["translation_table"],
            "orf_aa_length": len(seq),
            "RdRP_model": pol["model"],
            "RdRP_iE": pol["i_e"],
            "GDD_pos": gdd,
            "protease_models": ";".join(sorted({h["model"] for h in proteases})),
            "high_CG_pos": high_c,
            "high_method": high_method or "",
            "high_length": (gdd + 2 - high_c + 1) if high_c else None,
        }
        per_orf.append(rec)
        by_sid[sid].append(rec)

    high_lengths = [r["high_length"] for r in per_orf if r["high_length"]]
    if not high_lengths:
        raise RuntimeError("No high-confidence GxCG->GDD markers found")
    family_median = statistics.median(high_lengths)

    genus_high = defaultdict(list)
    for r in per_orf:
        if r["high_length"]:
            genus_high[r["virus_genus"]].append(r["high_length"])

    genus_target = {}
    for genus, vals in genus_high.items():
        if len(vals) >= 5:
            genus_target[genus] = statistics.median(vals)

    # ---------- Pass 2: choose one corrected marker per sequence ----------
    selected = []

    for sid, rows in by_sid.items():
        candidates = []

        for r in rows:
            oid = r["orf_id"]
            seq = seqs[oid]
            gdd = r["GDD_pos"]
            hs = hits.get(oid, [])
            pol = best_rdrp(hs)
            proteases = significant_upstream_proteases(hs, pol)

            if r["high_CG_pos"]:
                c = r["high_CG_pos"]
                method = r["high_method"]
                L = r["high_length"]
                # Strong preference for genuine GxCG candidates.
                score_class = 3 if method == "GxCG_near_protease_HMM" else 2
            else:
                # Fallback is the literal ICTV CG motif.
                cg = [c for c in plain_cg_starts(seq) if c < gdd and plausible(c, gdd)]
                if not cg:
                    continue

                target = genus_target.get(r["virus_genus"], family_median)

                # Protease-supported plain CG gets priority.
                cg_pro = [c for c in cg if near_protease(c, proteases)]
                pool = cg_pro if cg_pro else cg
                c = min(
                    pool,
                    key=lambda x: abs((gdd + 2 - x + 1) - target)
                )
                L = gdd + 2 - c + 1
                method = "CG_near_protease_HMM" if cg_pro else "CG_motif_length_guided"
                score_class = 1 if cg_pro else 0

            # Within one sequence, prefer stronger motif class, better RdRP E-value,
            # and length closer to the genus/family high-confidence median.
            target = genus_target.get(r["virus_genus"], family_median)
            e_score = -math.log10(max(float(r["RdRP_iE"]), 1e-300))
            len_penalty = abs(L - target)
            candidates.append((
                score_class,
                e_score,
                -len_penalty,
                oid,
                c,
                gdd,
                L,
                method,
                r
            ))

        if not candidates:
            continue

        candidates.sort(reverse=True)
        _, _, _, oid, c, gdd, L, method, r = candidates[0]
        seq = seqs[oid]
        marker = seq[c-1:gdd+2]  # include CG start through complete GDD

        sm = seqmaster.get(sid, {})
        selected.append({
            "sequence_id": sid,
            "virus_taxid": (sm.get("virus_taxid") or "").strip(),
            "virus_genus": r["virus_genus"],
            "virus_species": r["virus_species"],
            "orf_id": oid,
            "translation_table": r["translation_table"],
            "marker_start_CG_aa": c,
            "marker_end_GDD_aa": gdd + 2,
            "marker_length_aa": len(marker),
            "extraction_method": method,
            "RdRP_model": r["RdRP_model"],
            "RdRP_iE": f"{float(r['RdRP_iE']):.6g}",
            "protease_models": r["protease_models"],
            "host_annotated": "YES" if host_fam.get(sid) else "NO",
            "host_family_n": len(host_fam.get(sid, set())),
            "_marker": marker,
        })

    # ---------- outputs ----------
    fasta_path = outdir / "Secoviridae.ProPol_CG_to_GDD.corrected.faa"
    with open(fasta_path, "w") as f:
        for r in selected:
            f.write(
                f">{r['sequence_id']}|marker=ProPol_CG_to_GDD|method={r['extraction_method']}\n"
            )
            s = r["_marker"]
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")

    map_path = outdir / "Secoviridae.ProPol_CG_to_GDD.corrected.tsv"
    fields = [
        "sequence_id","virus_taxid","virus_genus","virus_species","orf_id",
        "translation_table","marker_start_CG_aa","marker_end_GDD_aa",
        "marker_length_aa","extraction_method","RdRP_model","RdRP_iE",
        "protease_models","host_annotated","host_family_n"
    ]
    with open(map_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in selected:
            w.writerow({k:r[k] for k in fields})

    lengths = [r["marker_length_aa"] for r in selected]
    species = {r["virus_species"] for r in selected if r["virus_species"] != "NA"}
    genera = {r["virus_genus"] for r in selected if r["virus_genus"] != "NA"}
    taxids = {r["virus_taxid"] for r in selected if r["virus_taxid"]}
    host_sids = {r["sequence_id"] for r in selected if r["host_annotated"] == "YES"}

    hf_union = set()
    ht_union = set()
    for r in selected:
        hf_union |= host_fam.get(r["sequence_id"], set())
        ht_union |= host_taxid.get(r["sequence_id"], set())

    methods = defaultdict(int)
    by_genus_selected = defaultdict(int)
    by_genus_species = defaultdict(set)
    for r in selected:
        methods[r["extraction_method"]] += 1
        by_genus_selected[r["virus_genus"]] += 1
        if r["virus_species"] != "NA":
            by_genus_species[r["virus_genus"]].add(r["virus_species"])

    summary_path = outdir / "Secoviridae.corrected_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric","value"])
        w.writerow(["input_family_sequences", 3431])
        w.writerow(["significant_RdRP_sequences", len(by_sid)])
        w.writerow(["corrected_ProPol_marker_sequences", len(selected)])
        w.writerow(["corrected_marker_recovery_pct_all_sequences", f"{100*len(selected)/3431:.4f}"])
        w.writerow(["corrected_marker_species_label_n", len(species)])
        w.writerow(["corrected_marker_genus_n", len(genera)])
        w.writerow(["corrected_marker_virus_taxid_n", len(taxids)])
        w.writerow(["host_annotated_marker_sequence_n", len(host_sids)])
        w.writerow(["host_annotation_pct", f"{100*len(host_sids)/len(selected):.4f}" if selected else "NA"])
        w.writerow(["host_taxid_n", len(ht_union)])
        w.writerow(["host_family_n", len(hf_union)])
        if lengths:
            w.writerow(["marker_length_min", min(lengths)])
            w.writerow(["marker_length_q05", f"{q(lengths,0.05):.1f}"])
            w.writerow(["marker_length_median", f"{statistics.median(lengths):.1f}"])
            w.writerow(["marker_length_q95", f"{q(lengths,0.95):.1f}"])
            w.writerow(["marker_length_max", max(lengths)])
        w.writerow(["high_confidence_family_median_length", f"{family_median:.1f}"])

    method_path = outdir / "Secoviridae.corrected_method_counts.tsv"
    with open(method_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["extraction_method","sequence_n"])
        for k,v in sorted(methods.items()):
            w.writerow([k,v])

    genus_path = outdir / "Secoviridae.corrected_genus_coverage.tsv"
    with open(genus_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["virus_genus","marker_sequence_n","marker_species_label_n"])
        for genus in sorted(by_genus_selected):
            w.writerow([genus, by_genus_selected[genus], len(by_genus_species[genus])])

    print("03K-B3 corrected Secoviridae extraction PASS")
    print(f"Output: {outdir}")
    print(f"significant_RdRP_sequences={len(by_sid)}")
    print(f"corrected_ProPol_marker_sequences={len(selected)}")
    print(f"species_labels={len(species)}")
    print(f"named_genera={len(genera)}")
    print(f"host_families={len(hf_union)}")
    if lengths:
        print(
            f"marker_length=min:{min(lengths)} "
            f"q05:{q(lengths,0.05):.1f} "
            f"median:{statistics.median(lengths):.1f} "
            f"q95:{q(lengths,0.95):.1f} "
            f"max:{max(lengths)}"
        )
    print("method_counts:")
    for k,v in sorted(methods.items()):
        print(k, v)
    print("genus_coverage:")
    for genus in sorted(by_genus_selected):
        print(genus, by_genus_selected[genus], len(by_genus_species[genus]))

if __name__ == "__main__":
    main()
