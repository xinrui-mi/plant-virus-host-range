# Reproducibility guide

This repository supports three *distinct* levels of reproduction. They are not
interchangeable; state which one you ran when you cite a number.

| Level | Question it answers | Ships with | Needs external data | Deterministic |
|---|---|---|---|---|
| **L1 — rebuild from raw sources** | Can the whole database be reconstructed from NCBI Virus, Virus-Host DB, Land Flora RNA virome and VIPDB? | pipeline scripts only | **yes**, all of them | interactive, **not** fully deterministic |
| **L2 — re-run the published analyses** | Do the released analysis scripts reproduce the released result tables from the released analysis data? | `data/` + `analysis/` | **no** | **mostly yes (seeded)** — with one known exception, see the Louvain note below |
| **L3 — results to figures** | Do the released result tables reproduce the submitted figures? | `analysis/` + plot scripts | **no** | **yes** — figures are drawn from the frozen result tables |

Most readers want **L2**, then **L3**. L1 is documented (Methods 2.1-2.2 and
`docs/EXTERNAL_INPUTS.md`) but is deliberately **not** claimed as a one-command
reproduction: the upstream harmonisation was interactive, spanned a local R
session and an HPC cluster, and involved manual curation decisions. Rebuilding
it after the fact would produce a *different* database, not a check on the
published one.

---

## L2 — re-run the published analyses (recommended)

Runs on a laptop with ~16 GB RAM. Nothing outside the repository is needed.

#### Determinism boundary within L2

Most L2 steps are seeded and re-run to the same values on the same inputs. The
one known exception is the **Louvain community detection** in
`analysis/04_genus_host_network` (`05A_network_main.R`,
`igraph::cluster_louvain()`): it is a randomised multi-start heuristic, so a
re-run can return a different (equally valid) modularity partition — in local
testing, 500 seeded draws produced module counts from 9 to 13. What *is*
reproducible: the edge weights, all partition-independent node metrics
(betweenness, strength, degree — identical to < 1e-15) and the resulting set of
12 prioritised genera. What is **not** guaranteed to reproduce value-for-value:
the exact 12-module partition and Q = 0.3532 reported in the paper. For that
reason the four network tables the figure reads are shipped as the **frozen
outputs the manuscript actually used** (MD5-matched to the analysis run), and
Fig. 4 always draws 12 modules from them. This is documented in
`FINAL_RELEASE_VALIDATION.md` §B4 and `UNRESOLVED_ITEMS.md` U20.

Also note what L2 does **not** cover: it re-runs analyses from the released
`data/` tables. Several of those tables are themselves derived from external
raw sources (NCBI Virus / GenBank records, Virus-Host DB, the Land Flora RNA
virome study, the VIP-DB workbook — see `docs/EXTERNAL_INPUTS.md`); rebuilding
*those* from raw data is L1 territory and is not part of L2.

```bash
git clone https://github.com/xinrui-mi/plant-virus-host-range.git
cd plant-virus-host-range

source analysis/00_config/00_paths.sh      # exports ANALYSIS_ROOT, REPO_ROOT, DATA_DIR, ...
bash   analysis/00_config/00_preflight.sh  # checks every released input and required column

Rscript analysis/01_dataset_summary/analysis/01_Fig1_data_framework_QC.R
Rscript analysis/01_dataset_summary/analysis/02_Fig2_global_landscape.R

Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03A_host_constraint_main.R
Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03B_degree_preserving_null.R
Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03C_sensitivity.R
Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03F_phylogeny_host_spectrum_structure.R
Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03G_phylogeny_host_spectrum_sensitivity.R

Rscript analysis/03_arthropod_transmission_ecology/analysis/04A_VIPDB_normalize.R
Rscript analysis/03_arthropod_transmission_ecology/analysis/04B_arthropod_transmission_analysis.R
Rscript analysis/03_arthropod_transmission_ecology/analysis/04C_host_distance_interaction_sensitivity.R
Rscript analysis/03_arthropod_transmission_ecology/analysis/04D_transmission_biology_extended_v3.R

Rscript analysis/04_genus_host_network/analysis/05A_network_main.R
Rscript analysis/04_genus_host_network/analysis/05B_network_robustness.R
Rscript analysis/04_genus_host_network/analysis/05C_RdRP65_network_sensitivity.R
Rscript analysis/04_genus_host_network/analysis/05D_projection_randomization.R
Rscript analysis/04_genus_host_network/analysis/05E_cross_genus_VC80_sensitivity.R

Rscript analysis/05_vector_contextualization/analysis/06A_candidate_external_context.R
Rscript analysis/05_vector_contextualization/analysis/06B_candidate_vector_support_adjusted_sensitivity.R

Rscript analysis/06_TSWV_evolution/analysis/07A_candidate_species_decomposition.R
Rscript analysis/06_TSWV_evolution/analysis/07N_segment_concordance_reassortment.R
```

Each script writes into `<module>/tables/` and `<module>/qc/`; nothing is
written back into `data/`.

Two scripts in module `02` produce files whose numeric prefixes overlap but whose names do not
collide, so running them in either order is safe and neither overwrites the other:

| Script | Writes |
|---|---|
| `03B_degree_preserving_null.R` | `06_formal_taxid_degree_preserving_null_{beta,distance}_{iterations,summary}.tsv`, `10_RdRP65_degree_preserving_null_*.tsv` |
| `03C_sensitivity.R` | `10_global_decay_sensitivity.tsv`, `11_global_decay_sensitivity_summary.tsv`, `12_categorical_distance_sensitivity_model.tsv` |

If a long run is interrupted, the output files it had already written are **stale and partial**.
Delete the affected table group and re-run the script rather than mixing outputs from two runs —
`03B` writes all eight of its files only at the end, so a partially populated set indicates an
interrupted run.

### Ordering constraints

1. `01` before `02` — `02` reads `01_dataset_summary/tables/01_VC80_internal_consistency_QC.tsv`.
2. `03A` before `03B`/`03C` — `03A` writes `02_host_relatedness_and_phylogeny/work/03_main_objects.rds`,
   which `03B` (line 5) and `03C` (line 3) read. That `work/` directory is git-ignored, so the file
   does **not** ship with the repository; it is a local cache created by `03A`.
3. `04A` before `04B`/`04C`/`04D` — the VIPDB normalization feeds all three.
4. `05A` before `05B`-`05E` and before `06A`/`06B` — the candidate thresholds come from `05A`.
5. All twelve `03L_FINAL_run_one_family.R` runs before `03L_FINAL_aggregate.R`.
6. `03F`/`03G` read the frozen trees directly, so they do **not** need `03K`.

### Memory and worker count

`03B_degree_preserving_null.R` runs iterations in bounded parallel chunks of `NCORES` workers. Peak
resident memory is about 480 MB per worker, so a machine with limited free RAM should lower the
worker count:

```bash
NCORES=4 Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03B_degree_preserving_null.R
```

`NCORES` affects only process management; each iteration derives its own seed from `SEED`, so the
null distribution does not depend on the worker count.

**Why this script validates every iteration.** An earlier version dispatched all iterations in a
single `mclapply` call. When a fork worker is reaped without returning a value, `mclapply`
substitutes `NULL` and does *not* raise an error — so a run can print `PASS` while using far fewer
resamples than requested. That happened once during release preparation: 999 iterations were
requested and only 251 were actually fitted, with the p-value shifting from 0.001 to 0.0040. The
current script therefore checks each returned element, recomputes any missing iteration serially, and
**aborts** if the number of distinct iterations still falls short. The QC table reports
`iterations_requested` and `iterations_completed` as separate columns, so a short run is visible
rather than silent. If you ever see those two columns differ, the run did not complete and its
p-value must not be quoted.

The full manuscript run uses `FIG3_NULL_ITER = 999` and `FIG3_RDRP_NULL_ITER = 999`, which takes
roughly 30–40 minutes at `NCORES=12` on a 32-core machine. A quick smoke test sets both to a small
number (see *Fast pass* below), which is enough to confirm that the pipeline runs but produces a
*p-value that must not be quoted*.

### Two scripts stop by design (a GATE, exit status != 0)

| script | GATE condition | what to read |
|---|---|---|
| `04A_VIPDB_normalize.R` | a VIPDB transmission mode cannot be mapped to one of the six audited raw forms | the gate table it writes into `tables/` |
| `06A_candidate_external_context.R` | a newly prioritised genus has no manually curated external record | `GATE_*.tsv` in `tables/` |

These are deliberate stop conditions, not bugs.

### Fast pass

Every randomised step is seeded (`SEED = 20260813` in
`analysis/00_config/00_config.R`). To check that the code runs end to end
without waiting for the full permutation counts, override them:

```bash
export FIG3_NULL_ITER=50 FIG3_RDRP_NULL_ITER=50 FIG3_SENSITIVITY_ITER=10 \
       FIG3_TREE_PERM_ITER=99 FIG5_ROBUSTNESS_ITER=20 FIG5_PROJECTION_NULL_ITER=99
```

Reduced-iteration output is for smoke-testing only. **Do not compare it with
the published numbers** — the published numbers use the counts frozen in
`00_config.R`.

---

## L3 — results to figures

The plot code is released separately from the analysis code because the figures
were type-set as journal artwork. Each figure has one entry script and one
`00_paths.R` that resolves its inputs relative to `$ANALYSIS_ROOT`.

| manuscript item | plot directory | entry script | one-command wrapper |
|---|---|---|---|
| Fig. 1 | `analysis/01_dataset_summary/plot/fig1/` | `01_Fig1_publication_plot_FINAL_v14.R` | (none — run with `Rscript`, see below) |
| Fig. 2 | `analysis/02_host_relatedness_and_phylogeny/plot/fig2/` | `03_Fig2_publication_plot_FINAL_v3.R` | `run_fig2_publication_FINAL_v12.sh` |
| Fig. 3 | `analysis/03_arthropod_transmission_ecology/plot/fig3/` | `03_Fig3_publication_plot_FINAL_v22.R` | `run_fig3_publication_FINAL_v22.sh` |
| Fig. 4 | `analysis/04_genus_host_network/plot/fig4/` | `04_Fig4_publication_FINAL_v1_1.R` | `run_Fig4_publication_FINAL_v1_1.sh` |
| Fig. 5 | `analysis/06_TSWV_evolution/plot/fig5/` | `07Q_Fig7_publication_plot_FINAL_v1_3_6.R` | `run_fig7_publication_FINAL_v1_3_6.sh` |
| Fig. S1 (RdRP65 host-sharing decay, single panel) | `analysis/02_host_relatedness_and_phylogeny/plot/figS1/` | `01_FigS1_publication_plot_v2_11_1_standalone.R` | `run_FigS1_standalone.sh` |

The supplementary figures S2–S5 are **not** separate scripts: they are rendered
by the same entry script as their main figure and land in that directory's
`supplementary/` output folder. Fig. S1 is the exception — it is the submitted
single-panel RdRP65 host-sharing decay figure (98 × 75 mm), drawn by a
**standalone entry** (`01_FigS1_…standalone.R`) that reuses the frozen plotting
logic (`make_distance_decay_plot()` from the pre-final code
`03_Fig3_publication_plot_v2_11_1.R`, kept alongside for provenance) and the
frozen input `tables/04_RdRP65_host_sharing_by_distance.tsv`. Its y-axis and
legend titles carry the submitted figure's wording ("RdRP65 sharing
probability (%)", "Number of host pairs"); the frozen code's internal copy of
the same panel (`FigS3A_RdRP65_host_taxonomic_decay`) differs only in those two
label strings. The submitted S1 reads no tree files. Output:
`main/FigS1_RdRP65_host_sharing_decay.{pdf,png,tiff}`, verified against the
submitted `Supplementary Fig1.pdf` (same data, geometry, colours and labels).
Fig. S2 (the twelve-family phylogeny panels) comes out of the Fig. 2 run;
Fig. S3 out of Fig. 3; Fig. S4 out of Fig. 4; Fig. S5 out of Fig. 5.

Directory names are lower-case (`fig1`, `fig2`, `fig3`, `fig4`, `fig5`,
`figS1`). Script file names keep the version suffixes of the frozen artwork runs;
they are deliberately **not** renamed to match the manuscript figure numbers.
Only `fig5/` is inconsistent — its script still carries the old internal name
`07Q_Fig7_…` because that is the frozen artefact's name.

```bash
source analysis/00_config/00_paths.sh
# Easiest: run the wrapper. It resolves its own directory and needs no arguments.
bash analysis/02_host_relatedness_and_phylogeny/plot/fig2/run_fig2_publication_FINAL_v12.sh
# Fig. 1 has no wrapper; point Rscript at the script, which finds the repo itself.
Rscript analysis/01_dataset_summary/plot/fig1/01_Fig1_publication_plot_FINAL_v14.R
```

Run these from the repository root, or from anywhere at all — each wrapper and
each `00_paths.R` resolves the repository layout from its own file location, so
`$ANALYSIS_ROOT` does not have to be exported first.

The plot scripts write PDF/PNG/TIFF into `<plot dir>/main/`,
`<plot dir>/panels/` and `<plot dir>/supplementary/`. Figure assembly (panel
composition, journal sizing) is done by the plotting code, but the exact
submitted artwork also passed through the journal's production step; minor
typographic differences are possible and are not a reproducibility failure.

If a level-3 run is all you need, you can skip the level-2 analysis entirely:
the frozen result tables are shipped and the plot scripts read them.

---

## What is *not* reproduced here

| item | why | where it is documented |
|---|---|---|
| Upstream harmonisation of the three sources + VIPDB | interactive; spans a local session and an HPC cluster | Methods 2.1; `docs/EXTERNAL_INPUTS.md` |
| Rebuilding the 12 family trees from raw sequences | needs the full `all_three_source.fasta`, the DNA GenBank archive and Pfam-A HMM, none of which can be redistributed | `analysis/02_*/pipeline/03K_*`; `docs/EXTERNAL_INPUTS.md` |
| Rebuilding the TSWV L/M/S trees and HyPhy runs | needs the TSWV source GenBank archive and HyPhy/GARD | `analysis/06_TSWV_evolution/pipeline/07C`-`07M` |
| `supplementary/Data_S1.xlsx` and Tables S1-S5 | frozen submission products, published exactly as submitted | `supplementary/README.md` |

The *outputs* of every one of these steps are shipped, so the L2 and L3 chains
never require them.

---

## Environment

See `environment/README.md`. The frozen versions are:

- R 4.4.3 (`data.table`, `ape`, `igraph`, `MASS`, `Matrix`, `readxl`)
- Python 3.12 (`pandas`, `numpy`, `biopython`)
- CD-HIT 4.8.1, MAFFT 7.526, trimAl 1.5.rev1, IQ-TREE 3.0.1, HMMER 3.4, HyPhy 2.5.101
  (phylogeny rebuilding only)

```bash
conda create -n pvh -c conda-forge -c bioconda \
  r-base=4.4.3 r-data.table r-ape r-igraph r-mass r-matrix r-readxl \
  python=3.12 pandas numpy biopython
conda activate pvh
```

---

## Verifying a run

`checksums/SHA256SUMS.txt` lists the digest of every released file, including
the frozen result tables that the L2 chain should reproduce.

```bash
sha256sum -c checksums/SHA256SUMS.txt            # release integrity
```

`FINAL_RELEASE_VALIDATION.md` records, module by module, which numbers were
re-derived from the frozen cluster outputs and which were not, with the
evidence for each verdict.
