# Public release sign-off

Pre-release audit of `plant-virus-host-range`, carried out as the last fixing pass
before publication. Every item below was checked against the frozen cluster
code, the final Data S1, and the manuscript's own figure and section numbering.
Nothing was recomputed to make a number match a document.

**Overall verdict: READY FOR PUBLIC RELEASE — with two narrow, documented
caveats.** All six figure wrappers (Fig. 1–5 and Fig. S1) and the twelve-family
aggregate now run end-to-end from a fresh clone. The two remaining caveats are
documented rather than hidden: (1) figure comparisons are runtime-level for
Fig. 1–5 and S2–S5 (exit codes, logged values, panel counts); Fig. S1 is the
one figure verified directly against the submitted PDF (item 9).

| # | Item | Before | After | Verdict |
|---|---|---|---|---|
| 1 | Twelve-family significance counts | `9 of 12` used in three places as if it were the main result | Three counts separated and labelled: main **11/12** (BH q), support-stratified **9/12**, ≥2-sequence-tip **12/12** (10/12 at ≥3) | PASS |
| 2 | Network modules 10 vs 12 | Validation report said `10 modules; Q = 0.3446`; release tables had the 10-module partition | Root cause identified as Louvain non-determinism, not a logic or data difference; frozen 12-module tables installed (MD5-matched); Fig. 4 now draws 12 / Q = 0.3532 | PASS, with a documented caveat |
| 3 | Figure numbering | `Fig. 2D–F` in 7 files; `Fig. 2F`; `Fig. 6` in two scripts; `Results 3.7`; TSWV `29 species` | All corrected to the final numbering (Fig. 2 has A–D only); `29` is now `VC80 clusters` | PASS |
| 4 | Reproducibility commands | Six of six plot paths wrong (case, invented file names, Fig. 3/4/5 pointed at the wrong module) | Paths and entry scripts rewritten from the actual tree; every command executed in a fresh clone | PASS for Fig. 1–5 and Fig. S1 (verified against the submitted PDF, item 9) |
| 5 | `MANUSCRIPT_OUTPUT_MAP.tsv` | S2E2/S2E3 attributed to `03G`; S2E4 was a stale pre-`03L` export (4 families) | Reattributed to the `03L` chain; S2E4 re-exported byte-for-byte from the frozen `04_FINAL_min_support_sensitivity_summary.tsv` (12 families × 2 thresholds, 24 rows; matches `Data_S1.xlsx`) — see `UNRESOLVED_ITEMS.md` U24 | PASS |
| 6 | Git-tracked content | Per-family `03L` tables were ignored although `03L_FINAL_aggregate.R` reads them — a clone could not rebuild the twelve-family block | Ignore rule narrowed to scratch only; 96 files added; aggregate verified to run in a clone with zero diff | PASS |
| 7 | Data directory layout | — | Unchanged: `data/curated/`, `data/analysis_data/`, `data/external_evidence/`. No `datasets/` added | PASS (no change made) |
| 8 | Internal file names | `prioritized_genera_external_evidence_FINAL.tsv`, `priority_genus_external_evidence_provenance_FINAL.tsv` | `prioritized_genera_evidence_summary.tsv`, `prioritized_genera_evidence_sources.tsv`; analysis scripts do not read either, so nothing broke | PASS |
| 9 | Fig. S1 inputs | Two 03G summary tables and the `03E` master were never shipped | The submitted S1 is the single-panel RdRP65 decay figure; standalone entry built (reuses frozen `make_distance_decay_plot()` + `tables/04_RdRP65_host_sharing_by_distance.tsv`), labels matched to the submitted wording, verified against the submitted `Supplementary Fig1.pdf`. Reads no tree files; former U23 withdrawn | PASS |
| 10 | Licensing and metadata | `LICENSE_DATA` scope paragraph could be read as licensing third-party records | Scope split into "what the authors produced" vs "what remains under upstream terms"; curated evidence described as traceable, not experimentally validated | PASS |

---

## 1. Twelve-family significance: three counts, not one

`S2E1_twelve_family_phylogenetic_summary.tsv` was read column by column
(`empirical_p_greater`, `BH_q_greater`), not taken from any summary sentence:

| Family | p (greater) | BH q (greater) | Significant |
|---|---|---|---|
| Alphaflexiviridae … Rhabdoviridae (8) | 0.001–0.004 | 0.0015–0.0053 | yes |
| Tombusviridae | 0.015 | 0.0180 | yes |
| Partitiviridae | 0.041 | 0.0447 | yes |
| Mitoviridae | 0.095 | 0.0950 | **no** |

That is **11 of 12**. The other two counts come from different analyses and are
now labelled as such:

- support-stratified sensitivity, `S2E2`, BH q < 0.05 → **9 of 12**;
- minimum support ≥2 sequences per tree tip,
  `03L/04_FINAL_min_support_sensitivity_summary.tsv` → **12 of 12** (≥3 → 10 of 12).

`S2E4_high_support_sensitivity.tsv` was found to be a stale pre-`03L` export
carrying the old four-family (old-tree) numbers; it has been re-exported
byte-for-byte from the frozen `04_FINAL_min_support_sensitivity_summary.tsv`
(12 families × 2 thresholds = 24 rows) and now matches `Data_S1.xlsx` and the
manuscript numbers exactly. See `UNRESOLVED_ITEMS.md` U24. The frozen tables
were not touched.

## 2. Network modules: why 10 ≠ 12, and what ships

`05A_network_main.R` calls `cluster_louvain()`. Louvain is a randomised
multi-start heuristic, so it returns *an* optimum, not *the* optimum.

| Check | Result |
|---|---|
| `02_network_edge_weights.tsv` MD5, cluster vs release | identical (`13abce55…`) |
| `01_network_input…VC80.tsv` MD5, cluster vs release | identical (`f00b540a…`) |
| Script diff `05A` (+ helper) | path resolution only |
| Master seed, both sides | `20260813L` |
| 500 draws on the shipped edges | 500 distinct partitions; module counts 9–13 (10: 90, 11: 252, 12: 128) |
| Frozen Q = 0.353235901330736 reachable in those 500 draws? | **no** |
| Betweenness / strength / degree, 269 genera | max abs diff **3.7e-16 / 0 / 0** |
| Prioritised candidate set | **identical, 12 of 12** |

The frozen Q lies inside the spread but does not reappear locally, which points
at the cluster's igraph build producing a different random stream — not at a
difference in data, parameters or code. Everything the conclusion rests on
(genera, edges, all partition-independent metrics, the 12 candidates) is
reproduced exactly; only the partition differs.

Because a re-run cannot land on the published partition, the four tables the
figure reads are shipped as the frozen cluster outputs, MD5-matched:

| File | MD5 (release = cluster) |
|---|---|
| `03_network_module_membership.tsv` | `171e8d46ce9604ea7726958c78185dc6` |
| `04_current_ICTV_genus_network_metrics.tsv` | `ea5fe9762bd0e6e5038e4af7bb21ca04` |
| `05_candidate_thresholds.tsv` | `3592eccb354d1ac7f1946c7741c94542` |
| `06_network_candidates.tsv` | `4b8b6f05d53df81fee4710bdd150c064` |

Fig. 4 re-rendered against them: exit 0, log reports **12**. The 10-module
re-run is documented in `FINAL_RELEASE_VALIDATION.md` §B4 rather than hidden.

## 3. Figure numbering and units

| Stale reference | Files affected | Corrected to |
|---|---|---|
| `Fig. 2D–F`, `Fig. 2D-F`, `Fig. 2D-E`, `Fig. 2F` | README, `supplementary/README.md`, output map, dependency closure, diff audit, validation, unresolved items | `Fig. 2C–D` (Fig. 2 has A–D only) |
| `Fig.6 accession-level`, `Fig.6` in `06B` comments | evidence matrix, `06B_…R` | `Fig. 4` |
| `Results 3.7` | diff audit row I added earlier | `Results 3.7` was correct for that row's manuscript mapping; the release itself uses §3.1–3.6 |
| `TSWV sequences / species / countries  1,260 / 29 / 30` | README | `VC80 clusters`, with a note that S5 counts isolates and gene datasets count sequences |

Fig. S1 vs S2 was also disentangled: **Fig. S1 = RdRp65 distance decay**,
**Fig. S2 = twelve-family phylogeny**. `analysis/PLOT_MAP.md` now carries an
explicit table for this, because the two are easy to swap.

## 4. Reproducibility commands

All six plot paths in `docs/reproducibility.md` were wrong before this pass:
upper-case directory names that do not exist (`plot/Fig1/` → real name
`plot/fig1/`), invented simplified script names (`01_Fig1_publication_plot.R` →
real `01_Fig1_publication_plot_FINAL_v14.R`), and Fig. 3/4/5 pointed at the
module that *reads* them rather than the module that *draws* them (Fig. 4 lives
under `04_genus_host_network`, not `05_vector_contextualization`).

Every row was rewritten from the actual tree and then executed in a **fresh
clone** (`git clone`, not the working copy):

| Command | Result |
|---|---|
| `Rscript analysis/01_dataset_summary/analysis/01_Fig1_data_framework_QC.R` (after `source 00_paths.sh`) | exit 0 — 80,028 / 83,034 / 50,807 / 72,494 |
| `03L_FINAL_aggregate.R analysis` | exit 0, **zero diff** against the tracked summaries |
| `run_fig2_publication_FINAL_v12.sh` (no env) | exit 0 — main + Fig. S2 3×4 and 2×6 |
| `run_fig3_publication_FINAL_v22.sh` (no env) | exit 0 |
| `run_Fig4_publication_FINAL_v1_1.sh` (no env) | exit 0 — log reports 12 modules |
| `run_fig7_publication_FINAL_v1_3_6.sh` (no env) | exit 0 — Fig. 5 complete |
| `run_FigS1_standalone.sh` (Fig. S1) | exit 0 — `FigS1_RdRP65_host_sharing_decay.{pdf,png,tiff}`; verified against the submitted `Supplementary Fig1.pdf` (item 9) |

The wrappers resolve the repository from their own location, so
`ANALYSIS_ROOT` does not have to be exported. Fig. 1 has no wrapper and needs
`source analysis/00_config/00_paths.sh` first — that is stated in the doc.

These are **runtime checks**, not full numeric re-derivations. No figure was
compared against the submitted artwork pixel by pixel.

## 5. Output map corrections

`S2E2` and `S2E3` were attributed to `03G_phylogeny_host_spectrum_sensitivity.R`.
They are not produced there. `03L_FINAL_aggregate.R` reads
`05_support_stratified_summary.tsv` and `04_pair_sampling_stability_summary.tsv`
from each of the twelve family directories and writes
`03_FINAL_support_stratified_summary.tsv` / `02_FINAL_pair_sampling_stability_summary.tsv`,
which are the sources of S2E2 and S2E3. `03G` is an independent earlier
implementation and is kept as such. `03F`/`03G` were not deleted and no
sensitivity analysis was silently dropped.

## 6. What Git actually tracks

Checked in a throwaway clone, not by looking at the tarball — a file being
present in the archive does not mean Git will publish it.

This found a genuine defect the previous passes had missed: the per-family
`03L_FINAL_12family_host_spectrum/<Family>/` directories were excluded by a
`*/*/**` rule, yet `03L_FINAL_aggregate.R` reads `01_main_summary.tsv`,
`04_…`, `05_…`, `07_…` from inside them and stops with a named-file error if
they are absent. In a clone, the twelve-family block was unrebuildable.

Fix: the rule now excludes only the scratch it was meant to exclude
(`*iteration*.tsv`, `*.fasta`, `*.faa`), and 96 per-family tables (3.0 MB) are
tracked. Verified after the fix: aggregate runs in a clone, exit 0, zero diff.

Tracked totals after the pass:

| Asset class | Tracked in clone |
|---|---|
| Files total | **629** |
| Twelve-family treefiles | 12 |
| `.nr95.membership.tsv` / `.alignment_sequence_QC.tsv` | 4 / 4 (this is the real extent of the shipped set) |
| TSWV five-gene trees / segment trees | 5 / 3 |
| `03L` per-family + pooled summaries | 96 + 7 |
| `data/**/*.gz` | 26 |
| `data/external_evidence/*.tsv` | 2 |
| `supplementary/Data_S1.xlsx` | 1 |

## 7–8. Data layout and public names

No `datasets/` directory was created; `data/curated/`, `data/analysis_data/` and
`data/external_evidence/` are as they were. The two `.gz` files under
`data/analysis_data/` that still carry the `20_…_FINAL` prefix are named that way
because the analysis scripts reference those exact paths; they were left alone
rather than renamed into something that would break the chain. The two
reader-facing evidence files were renamed because nothing reads them by name:

| Before | After |
|---|---|
| `prioritized_genera_external_evidence_FINAL.tsv` | `prioritized_genera_evidence_summary.tsv` |
| `priority_genus_external_evidence_provenance_FINAL.tsv` | `prioritized_genera_evidence_sources.tsv` |

Checked before renaming: `06A_candidate_external_context.R` writes its own
matrix and reads neither file, so no script depends on the old names. All
references in `docs/MANUSCRIPT_OUTPUT_MAP.tsv` and `data/README.md` were updated.

## 9. Fig. S1: resolved — standalone entry, verified against the submitted PDF

The submitted Supplementary Fig. S1 is the **single-panel RdRP65 host-sharing
decay figure** (98 × 75 mm), not the old four-panel composite. This pass:

1. **Built a standalone entry**,
   `plot/figS1/01_FigS1_publication_plot_v2_11_1_standalone.R` (wrapper
   `run_FigS1_standalone.sh`), which reuses the frozen helper
   `make_distance_decay_plot()` from the pre-final code
   (`03_Fig3_publication_plot_v2_11_1.R`, kept untouched for provenance) and
   the frozen input `tables/04_RdRP65_host_sharing_by_distance.tsv` (shipped
   since the first release). The frozen code's four-panel main figure and its
   phylogram atlases are NOT part of the submitted S1 and are not drawn here.
2. **Matched the submitted labels**: y-axis "RdRP65 sharing probability (%)"
   and size-legend title "Number of host pairs" (the frozen code's internal
   copy of the same panel — `FigS3A_RdRP65_host_taxonomic_decay` — says
   "Probability of sharing RdRP65 units (%)" and "Host-pair count").
3. **Verified against the submitted `Supplementary Fig1.pdf`**: page size
   matches (277 × 212 pt = 98 × 75 mm), and the rendered comparison shows the
   same data points, decay line, colours, bubble scaling, axis categories and
   label wording. **Fig. S1 reads no tree files** — the earlier concern about
   the four-family trees (former U23) does not apply to the submitted S1 and
   was withdrawn; the provisional `S1_frozen_trees/` directory was removed.

Runs end-to-end from a clone (exit 0, three outputs: PDF/PNG/TIFF).

## 10. Licensing and public safety

- `LICENSE_DATA` now separates what the authors produced (CC BY 4.0) from
  upstream records (NCBI Virus, Virus-Host DB, Land Flora, VIP-DB, ICTV, EPPO,
  FAO/IPPC), which stay under their own terms. It no longer reads as a blanket
  CC BY 4.0 claim over third-party content.
- Curated evidence is described as **traceable**, with each row pointing at a
  public source, and explicitly not as independently validated.
- The internal archive document name and its page ranges were removed from the
  evidence file. All 28 rows retain a public URL and a verification date, so no
  row lost its origin.
- Credential scan over the whole tree: no password, token, API key or private
  key. `token` appears only as a string-parsing variable in three scripts.
- Personal path scan: no `/Users/`, `~/Desktop/`, `mxrconda` in executable code.
  The only remaining `/share/home/...` strings are in three LSF submission
  scripts under `cluster_job_scripts/`, each marked `ARCHIVAL ONLY`.

---

## Residual risks

| Risk | Severity | Note |
|---|---|---|
| Louvain module count differs on re-run | mitigated | Documented in README, §B4 and `UNRESOLVED_ITEMS.md` U20; frozen tables are also now protected: a re-run of `05A` writes to `tables/rerun/` and cannot overwrite the published partition, and the Fig. 4 script verifies the frozen MD5 before plotting (U26) |
| Upstream licences verified 2026-10-27 | resolved | `docs/EXTERNAL_INPUTS.md` §5 records each source's terms and verification URL; Virus-Host DB content stays CC BY-NC-SA 4.0 and is excluded from the CC BY 4.0 grant; no raw third-party dump is redistributed |
| Figure comparisons are exit-code level for Fig. 1–5, S2–S5 | moderate | Only Fig. S1 was compared directly against the submitted PDF; no pixel comparison against the submitted artwork for the others; the journal's production step may differ typographically |
| `data/analysis_data/20_*_FINAL.tsv.gz` names | cosmetic | Kept because scripts reference those exact paths |

Not pushed anywhere. The archive in `/workspace` is the deliverable for review.
