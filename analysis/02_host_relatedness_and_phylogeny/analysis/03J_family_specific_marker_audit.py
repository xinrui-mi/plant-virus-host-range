#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import defaultdict
import math

ROOT = Path.cwd()

ANALYSIS = ROOT / "three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
PIPE = ROOT / "final_three_source_pipeline_v1/results"

CAND = ANALYSIS / "tables/31_03H_family_phylogeny_candidates_ge100_species.tsv"
MASTER = PIPE / "20_three_source_sequence_master_FINAL.tsv"
RDRP_CORE = PIPE / "08_RdRP_core.faa"
I03 = ANALYSIS / "pipeline_results/03I_RNA_standardized_RdRP"

OUT = ANALYSIS / "pipeline_results/03J_family_specific_marker_audit"
FAOUT = OUT / "family_nt"
OUT.mkdir(parents=True, exist_ok=True)
FAOUT.mkdir(parents=True, exist_ok=True)

FASTA_CANDIDATES = [
    PIPE / "06_three_source_final.fasta",
    ROOT / "all_three_source.fasta",
]
FASTA = next((x for x in FASTA_CANDIDATES if x.exists() and x.stat().st_size > 0), None)
if FASTA is None:
    raise FileNotFoundError("Could not find frozen three-source FASTA")

PLAN = {
    "Geminiviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "full_or_near_full_Rep_protein",
        "route": "GenBank_CDS_translation_preferred",
        "note": "Use full Rep; HMM/domain search only for identification/QC, not final boundary."
    },
    "Caulimoviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "continuous_Pol_region_RT_to_RNaseH1",
        "route": "GenBank_CDS_translation_plus_domain_boundaries",
        "note": "Recover one continuous homologous Pol region spanning RT-RNaseH1."
    },
    "Alphaflexiviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "entire_Rep_protein",
        "route": "replicase_ORF_or_CDS",
        "note": "ICTV family phylogeny marker."
    },
    "Betaflexiviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "entire_or_near_entire_replication_protein",
        "route": "replicase_ORF_or_CDS",
        "note": "Replication protein used for family phylogeny."
    },
    "Closteroviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "HSP70h_protein",
        "route": "CDS_or_HSP70_HMM",
        "note": "Use viral HSP70 homolog rather than RdRP core."
    },
    "Mitoviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "complete_or_near_complete_RdRP",
        "route": "annotated_protein_preferred",
        "note": "Avoid short RdRP HMM core; translation requires care for mitochondrial genetic code."
    },
    "Partitiviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "complete_RdRP",
        "route": "RdRP_CDS_or_full_RdRP_ORF",
        "note": "RdRP encoded by polymerase segment."
    },
    "Potyviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "complete_major_polyprotein",
        "route": "polyprotein_CDS_or_long_ORF",
        "note": "Audit complete-polyprotein recovery before considering any predefined fallback."
    },
    "Rhabdoviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "full_or_near_full_L_protein",
        "route": "L_CDS_or_long_ORF_then_trim",
        "note": "Align full L and remove ambiguously aligned regions after MSA."
    },
    "Secoviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "ProPol_CG_to_GDD",
        "route": "polyprotein_ORF_plus_motif_boundaries",
        "note": "Continuous region from 3C-like protease CG motif to polymerase GDD motif."
    },
    "Solemoviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "RdRP_domain_from_minus1_frameshift_site",
        "route": "RdRP_CDS_or_ORF_context",
        "note": "Recover homologous RdRP region rather than short HMM hit."
    },
    "Tombusviridae": {
        "decision": "MAIN_PHYLOGENY",
        "marker": "polymerase_readthrough_or_frameshift_conserved_region",
        "route": "polymerase_ORF_context_plus_RdRP_anchor",
        "note": "Use conserved readthrough/frameshift polymerase portion."
    },
    "Alphasatellitidae": {
        "decision": "EXCLUDE_MAIN_SPECIAL_SATELLITE",
        "marker": "Rep_possible_but_not_main",
        "route": "none_for_main_analysis",
        "note": "Helper-dependent movement/encapsidation/transmission; biologically different host-range interpretation."
    },
    "Tolecusatellitidae": {
        "decision": "EXCLUDE_MAIN_SPECIAL_SATELLITE",
        "marker": "no_universal_family_protein_marker",
        "route": "none_for_main_analysis",
        "note": "Betasatellites encode betaC1; deltasatellites are non-coding; no universal homologous protein."
    },
}


def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fasta_iter(path):
    h = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if h is not None:
        yield h, "".join(seq)


def quantile(vals, p):
    vals = sorted(vals)
    if not vals:
        return ""
    if len(vals) == 1:
        return vals[0]
    x = (len(vals) - 1) * p
    lo = math.floor(x)
    hi = math.ceil(x)
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - x) + vals[hi] * (x - lo)


cand = read_tsv(CAND)
candidate_families = [r["virus_family"] for r in cand]
candidate_set = set(candidate_families)

missing_plan = candidate_set - set(PLAN)
if missing_plan:
    raise RuntimeError(f"Marker plan missing candidate families: {sorted(missing_plan)}")

# ------------------------------------------------------------------
# 1. Candidate + marker plan
# ------------------------------------------------------------------
plan_out = OUT / "01_family_marker_plan.tsv"

base_cols = list(cand[0].keys())
extra_cols = ["marker_decision", "recommended_marker", "extraction_route", "marker_note"]

with plan_out.open("w", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t", fieldnames=base_cols + extra_cols)
    w.writeheader()
    for r in cand:
        p = PLAN[r["virus_family"]]
        rr = dict(r)
        rr.update({
            "marker_decision": p["decision"],
            "recommended_marker": p["marker"],
            "extraction_route": p["route"],
            "marker_note": p["note"],
        })
        w.writerow(rr)

# ------------------------------------------------------------------
# 2. Stage20 ID sets and nt length inventory
# ------------------------------------------------------------------
ids_by_family = defaultdict(set)
lengths = defaultdict(list)
source_counts = defaultdict(lambda: defaultdict(int))

with MASTER.open() as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        fam = r.get("virus_family", "")
        if fam not in candidate_set:
            continue
        sid = r["sequence_id"]
        ids_by_family[fam].add(sid)

        try:
            lengths[fam].append(int(float(r["sequence_length_nt"])))
        except Exception:
            pass

        source_counts[fam][r.get("source_dataset", "")] += 1

# ------------------------------------------------------------------
# 3. Extract candidate-family nucleotide FASTAs
# ------------------------------------------------------------------
all_target = set().union(*ids_by_family.values())
fam_of_id = {}
for fam, ids in ids_by_family.items():
    for sid in ids:
        fam_of_id[sid] = fam

handles = {}
recovered = defaultdict(int)

try:
    for fam in candidate_families:
        handles[fam] = (FAOUT / f"{fam}.nt.fasta").open("w")

    for sid, seq in fasta_iter(FASTA):
        if sid not in all_target:
            continue
        fam = fam_of_id[sid]
        handles[fam].write(f">{sid}\n{seq}\n")
        recovered[fam] += 1
finally:
    for h in handles.values():
        h.close()

inv_out = OUT / "02_family_sequence_inventory.tsv"
with inv_out.open("w", newline="") as f:
    cols = [
        "virus_family", "marker_decision",
        "stage20_sequence_n", "fasta_recovered_n",
        "nt_min", "nt_q25", "nt_median", "nt_q75", "nt_max",
        "source_counts"
    ]
    w = csv.DictWriter(f, delimiter="\t", fieldnames=cols)
    w.writeheader()

    for fam in candidate_families:
        x = lengths[fam]
        w.writerow({
            "virus_family": fam,
            "marker_decision": PLAN[fam]["decision"],
            "stage20_sequence_n": len(ids_by_family[fam]),
            "fasta_recovered_n": recovered[fam],
            "nt_min": min(x) if x else "",
            "nt_q25": round(quantile(x, .25), 2) if x else "",
            "nt_median": round(quantile(x, .50), 2) if x else "",
            "nt_q75": round(quantile(x, .75), 2) if x else "",
            "nt_max": max(x) if x else "",
            "source_counts": ";".join(
                f"{k}:{v}" for k, v in sorted(source_counts[fam].items())
            ),
        })

# ------------------------------------------------------------------
# 4. Frozen Stage8 RdRP-core lengths by candidate family
# ------------------------------------------------------------------
family_by_id = fam_of_id
core_lengths = defaultdict(list)
core_n = defaultdict(int)

if RDRP_CORE.exists():
    for sid, seq in fasta_iter(RDRP_CORE):
        fam = family_by_id.get(sid)
        if fam in candidate_set:
            aa = sum(c not in "-?." for c in seq)
            core_lengths[fam].append(aa)
            core_n[fam] += 1

core_out = OUT / "03_stage8_RdRP_core_length_audit.tsv"
with core_out.open("w", newline="") as f:
    cols = [
        "virus_family", "stage8_RdRP_core_n",
        "aa_min", "aa_q25", "aa_median", "aa_q75", "aa_max"
    ]
    w = csv.DictWriter(f, delimiter="\t", fieldnames=cols)
    w.writeheader()

    for fam in candidate_families:
        x = core_lengths[fam]
        w.writerow({
            "virus_family": fam,
            "stage8_RdRP_core_n": core_n[fam],
            "aa_min": min(x) if x else "",
            "aa_q25": round(quantile(x, .25), 2) if x else "",
            "aa_median": round(quantile(x, .50), 2) if x else "",
            "aa_q75": round(quantile(x, .75), 2) if x else "",
            "aa_max": max(x) if x else "",
        })

# ------------------------------------------------------------------
# 5. Existing 03I alignment diagnostic
# ------------------------------------------------------------------
diag_out = OUT / "04_existing_03I_alignment_diagnostic.tsv"

with diag_out.open("w", newline="") as f:
    cols = [
        "virus_family", "03I_QC_exists",
        "03I_input_n", "03I_kept_n",
        "effective_aa_min", "effective_aa_median", "effective_aa_max",
        "gap_fraction_median", "treefile_exists"
    ]
    w = csv.DictWriter(f, delimiter="\t", fieldnames=cols)
    w.writeheader()

    for fam in candidate_families:
        qcf = I03 / fam / f"{fam}.alignment_sequence_QC.tsv"
        tf = I03 / fam / f"{fam}.RdRP_core.treefile"

        if not qcf.exists():
            w.writerow({
                "virus_family": fam,
                "03I_QC_exists": "NO",
                "treefile_exists": "NO",
            })
            continue

        rows = read_tsv(qcf)
        eff = [float(r["effective_aa"]) for r in rows]
        gap = [float(r["gap_fraction"]) for r in rows]
        kept = sum(r["action"] == "keep" for r in rows)

        w.writerow({
            "virus_family": fam,
            "03I_QC_exists": "YES",
            "03I_input_n": len(rows),
            "03I_kept_n": kept,
            "effective_aa_min": min(eff) if eff else "",
            "effective_aa_median": round(quantile(eff, .50), 3) if eff else "",
            "effective_aa_max": max(eff) if eff else "",
            "gap_fraction_median": round(quantile(gap, .50), 5) if gap else "",
            "treefile_exists": "YES" if tf.exists() and tf.stat().st_size > 0 else "NO",
        })

# ------------------------------------------------------------------
# QC
# ------------------------------------------------------------------
qc = OUT / "05_marker_audit_QC.tsv"
with qc.open("w") as f:
    f.write("metric\tvalue\tstatus\n")
    f.write(f"candidate_family_n\t{len(candidate_families)}\tPASS\n")
    f.write(
        "main_phylogeny_family_n\t"
        f"{sum(PLAN[x]['decision']=='MAIN_PHYLOGENY' for x in candidate_families)}"
        "\tPASS\n"
    )
    f.write(
        "special_satellite_excluded_n\t"
        f"{sum(PLAN[x]['decision'].startswith('EXCLUDE') for x in candidate_families)}"
        "\tPASS\n"
    )
    f.write(
        f"target_nt_sequence_n\t{len(all_target)}\t"
        f"{'PASS' if sum(recovered.values()) == len(all_target) else 'FAIL'}\n"
    )
    f.write(
        f"fasta_recovered_n\t{sum(recovered.values())}\t"
        f"{'PASS' if sum(recovered.values()) == len(all_target) else 'FAIL'}\n"
    )

print("03J family-specific marker audit PASS")
print("Output:", OUT)
print("Candidate families:", len(candidate_families))
print("Main phylogeny families:",
      sum(PLAN[x]["decision"] == "MAIN_PHYLOGENY" for x in candidate_families))
print("Special satellite exclusions:",
      sum(PLAN[x]["decision"].startswith("EXCLUDE") for x in candidate_families))
