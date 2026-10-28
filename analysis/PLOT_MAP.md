# Plotting code

This directory holds the scripts that turn the released analysis tables into the manuscript figures.
It is organised by **manuscript figure number**, because that is how a reader looks for a figure. The
cluster working directories are organised by *analysis module*, so the two numbering schemes are
offset — this file records the mapping so no one has to guess.

## Numbering: module directory vs manuscript figure

| Manuscript | Cluster plotting directory | Analysis module it reads from | Released here |
|---|---|---|---|
| Fig. 1 | `Fig2_global_landscape/Fig1_publication_FINAL_v14/` | `01_dataset_summary` | `01_dataset_summary/plot/fig1/` |
| Fig. 2 | `Fig3_host_phylogenetic_constraint/Fig2_publication_FINAL_v17/` | `02_host_relatedness_and_phylogeny` | `02_host_relatedness_and_phylogeny/plot/fig2/` |
| Fig. 3 | `Fig4_arthropod_transmission_ecology/Fig3_transmission_publication_FINAL_v22/` | `03_arthropod_transmission_ecology` | `03_arthropod_transmission_ecology/plot/fig3/` |
| Fig. 4 | `Fig5_network_prioritization/Fig4_publication_FINAL_v1_1/` | `04_genus_host_network` | `04_genus_host_network/plot/fig4/` |
| Fig. 5 | `Fig7_TSWV_case_study/Fig7_publication_FINAL_v1_3_6/` | `06_TSWV_evolution` + `05_vector_contextualization` | `06_TSWV_evolution/plot/fig5/` |
| Fig. S1 | `Fig3_host_phylogenetic_constraint/plot_no/Fig3_publication_plot_v2_11_1/` | `02_host_relatedness_and_phylogeny` | `02_host_relatedness_and_phylogeny/plot/figS1/` |
| Fig. S2 | `Fig3_host_phylogenetic_constraint/Fig2_publication_FINAL_v17/` | `02_host_relatedness_and_phylogeny` | same directory as Fig. 2 |
| Fig. S3 | `Fig4_arthropod_transmission_ecology/Fig3_transmission_publication_FINAL_v22/` | `03_arthropod_transmission_ecology` | same directory as Fig. 3 |
| Fig. S4 | `Fig5_network_prioritization/Fig4_publication_FINAL_v1_1/` | `04_genus_host_network` | same directory as Fig. 4 |
| Fig. S5 | `Fig7_TSWV_case_study/Fig7_publication_FINAL_v1_3_6/` | `06_TSWV_evolution` | same directory as Fig. 5 |

**What each supplementary figure actually is** — these two are easy to confuse
and are *not* interchangeable:

| Supplementary | Content | Produced by |
|---|---|---|
| **Fig. S1** | **Single panel**: RdRP65 host-sharing probability (%) vs host taxonomic distance (same species → same kingdom), bubble size = number of host pairs, 98 × 75 mm. Drawn by the standalone entry `01_FigS1_publication_plot_v2_11_1_standalone.R` (wrapper `run_FigS1_standalone.sh`), which reuses the frozen helper `make_distance_decay_plot()` and the frozen input `tables/04_RdRP65_host_sharing_by_distance.tsv`. Verified against the submitted `Supplementary Fig1.pdf`: same data, geometry, colours, y-axis title ("RdRP65 sharing probability (%)") and legend title ("Number of host pairs"); the frozen code's internal copy of the same panel (saved as `FigS3A_RdRP65_host_taxonomic_decay`) differs only in those two label strings. **Fig. S1 reads no tree files.** | its own `plot/figS1/` directory |
| **Fig. S2** | Twelve-family phylogenetic / host-spectrum atlas | the Fig. 2 run (lands in `fig2/supplementary/`) |

The phylogram panels inside the frozen `03_Fig3_publication_plot_v2_11_1.R` script
resolve trees via `find_tree_for_family()` (final 12-family set), but those panels
belong to the superseded Fig. 3 design — the submitted Fig. S1 does not execute
them and reads no tree files. Fig. S3–S5 come out of the Fig. 3, 4
and 5 runs respectively.

All plot directories are lower-case (`fig1`, `fig2`, `fig3`, `fig4`, `fig5`,
`figS1`). The script file names keep the version suffixes of the frozen artwork
runs and are deliberately not renamed to match manuscript figure numbers — in
particular `fig5/` still carries the internal name `07Q_Fig7_…`.

The offset exists because the analysis modules were numbered as the analysis was built up, while the
figures were renumbered during writing. Nothing is wrong with either scheme; they simply count
different things.

## Files per directory

The directories are not identical. The two "legacy" name pairs (`fig2` / `figS1`) come from an older
template; the other three carry one extra file. In every case the main script and `00_paths.R` are
present.

| Directory | Main script | Auxiliary files |
|---|---|---|
| `01_dataset_summary/plot/fig1/` | `01_Fig1_publication_plot_FINAL_v14.R` | `00_check_packages.R`, `00_display_rules.R`, `00_helpers.R`, `00_paths.R`, `00_plot_theme.R` |
| `02_.../plot/fig2/` | `03_Fig2_publication_plot_FINAL_v3.R` | `00_check_packages.R`, `00_helpers_Fig2_v3.R`, `00_helpers.R`, `00_paths.R`, `00_plot_theme.R`, `README.txt`, `FINAL_FREEZE_NOTE.txt` |
| `02_.../plot/figS1/` | `03_Fig3_publication_plot_v2_11_1.R` | `00_check_packages.R`, `00_helpers.R`, `00_paths.R`, `00_plot_theme.R`, `README.txt` |
| `03_.../plot/fig3/` | `03_Fig3_publication_plot_FINAL_v22.R` | `00_check_packages.R`, `00_display_rules.R`, `00_helpers.R`, `00_paths.R`, `00_plot_theme.R`, `README.txt`, `FINAL_FREEZE_NOTE.txt` |
| `04_.../plot/fig4/` | `04_Fig4_publication_FINAL_v1_1.R` | `00_check_packages.R`, `00_display_rules.R`, `00_helpers.R`, `00_paths.R`, `00_plot_theme.R`, `README.txt`, `FINAL_FREEZE_NOTE.txt`, `PATCH_NOTE_v1_1.txt` |
| `06_.../plot/fig5/` | `07Q_Fig7_publication_plot_FINAL_v1_3_6.R` | `00_check_packages.R`, `00_helpers.R`, `00_paths.R`, `00_theme.R`, `README.txt` |

Each directory except `01_.../plot/fig1/` also has an executable `run_*.sh` that resolves the
repository root from its own location and calls `Rscript` on the main script. Fig. 1 has no wrapper:
its plotting script is self-configuring (`00_paths.R` resolves the repository root from
`sys.frame(1)$ofile`) and is run directly as

```bash
Rscript analysis/01_dataset_summary/plot/fig1/01_Fig1_publication_plot_FINAL_v14.R
```

which writes to the git-ignored `figures/fig1/{main,panels,logs}/`.

`cluster_job_scripts/submit_*.sh` is the **original LSF submission wrapper** from the computing
cluster (`bsub` flags, queue names, conda interpreter path). It is archived for audit only, carries a
header saying so, and is not part of the reproduction path.

`*_publication_*` means **plotting code**. A sibling directory named `publication_figures_*` is by
contrast **rendered output** (PDF/TIFF/PNG plus logs) and is not included in this repository; those
files are the figure files submitted with the manuscript, not code.

## Do the figures read the released tables, or do they hard-code numbers?

They read the released tables. Verified for Fig. 1:

```r
raw01 <- read_tsv_locked("01_virus_family_evidence.tsv")
raw02 <- read_tsv_locked("02_host_family_evidence.tsv")
raw03 <- read_tsv_locked("03_country_evidence.tsv")
raw04 <- read_tsv_locked("04_continent_evidence.tsv")
raw07 <- read_tsv_locked("07_evidence_rarefaction_summary.tsv")
```

`read_tsv_locked()` resolves every path through `TABLE_DIR` and raises an error if the file is absent,
so a figure cannot be drawn from stale or partial numbers. This is what makes level L3
(result → figure) of `docs/reproducibility.md` meaningful: re-running the plotting script on a
re-run table set reproduces the panel.

## Cluster paths

The plotting scripts originally carried their input locations as `Sys.getenv(..., unset = "/share/...")`
defaults inside `00_paths.R`, plus a `R_LIB` pointing at the cluster's conda R library and a
`R_LIBS_USER` export in the `run_*.sh` wrappers. All of that has been removed:

| Was | Now |
|---|---|
| `R_LIB <- "/share/home/.../mxrconda/lib/R/library"` | `R_LIB <- Sys.getenv("R_LIB", unset = "")`, applied only if set and existing |
| `DEFAULT_*_ROOT <- "/share/home/.../three_source_manuscript_analysis_v1_1/FigXX_..."` | derived from the script's own location (`sys.frame(1)$ofile` / `--file=`), four levels up to the repo root |
| `export R_LIBS_USER="/share/home/.../mxrconda/lib/R/library"` | a comment; honoured only if the user exports it |
| `RSCRIPT="${RSCRIPT:-/share/home/.../mxrconda/bin/Rscript}"` | `RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"` |

Every path still honours the original environment variable (`FIG3_ROOT`, `FIG7_ROOT`,
`NETWORK_ROOT`, `CONTEXT_ROOT`, `ASSOC_MASTER`, `FINAL_THREE_SOURCE_ROOT`, `FIG3_ANALYSIS_ROOT`,
`FINAL_PIPELINE_RESULTS`, `PROJECT_ROOT`, `R_LIB`), so a user with a different layout can override
without editing code. A grep for `share/home`, `bgi_mixr`, `mxrconda`, `/Users/` or `~/Desktop`
across all released scripts returns nothing **except** inside `cluster_job_scripts/submit_*.sh`,
which are the archived LSF submission scripts described above: those necessarily retain the
cluster's absolute paths, are marked `ARCHIVAL ONLY` in their header, and are never invoked by the
reproduction path. Nothing in the reproduction path contains a machine-specific path.

## Status of the figure code

All six plotting directories are now present: Fig. 1, Fig. 2 (+ S2), Fig. 3 (+ S3), Fig. 4 (+ S4),
Fig. 5 (+ S5) and Fig. S1. See `UNRESOLVED_ITEMS.md` section U5.
