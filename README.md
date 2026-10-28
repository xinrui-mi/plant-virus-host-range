# Plant virus host range: data and analysis repository

Companion repository for the manuscript
**"Phylogenetic structure, lineage differences and transmission ecology of global
plant virus host associations"**.

The study integrates **80,028 plant-virus sequences** and **83,034 virus–host
associations** (76 virus families, 296 host families, 184 countries) to examine
the global structure of plant virus host ranges. Six result blocks are covered:

1. the composition and evidence coverage of the integrated database (§3.1, Fig. 1);
2. the decay of virus sharing with host taxonomic distance, and the 24
   family-stratified slopes (§3.2, Fig. 2A–C);
3. the relationship between host spectrum and phylogenetic position across
   **12 family-specific marker trees** (§3.3, Fig. 2D, Fig. S2);
4. arthropod-associated sequence evidence and the transmission-ecology
   contrasts, including the adjusted host-range and country-breadth models and a
   within-mode lineage-composition sensitivity (§3.4, Fig. 3, Fig. S3);
5. the current-genus × host-family weighted network and the prioritised genera,
   with resampling, threshold-perturbation and small-world controls (§3.5, Fig. 4,
   Fig. S4);
6. tomato spotted wilt virus (TSWV; current ICTV species
   *Orthotospovirus tomatomaculae*) as a gene- and segment-level case study
   (§3.6, Fig. 5, Table S5).

---

## Data at a glance

| Quantity | Value |
|---|---|
| Virus sequences | 80,028 |
| Virus–host associations | 83,034 |
| Source composition (sequences) | NCBI Virus 74,100 · Virus-Host DB 919 · shared/coordinated 2,493 · Land Flora RNA virome 2,516 |
| VC80 units (retained) | 11,708 |
| Sequences with an RdRP65 assignment | 27,549 → 1,700 RdRP65 units |
| Virus families / host families / countries | 76 / 296 / 184 |
| Arthropod-evidence union (S0 ∪ S1 ∪ S2) | 4,401 sequences · 340 virus TaxIDs |
| VIP-DB documented-vector virus TaxIDs | 410 |
| Network genera / host families / weighted edges | 269 / 197 / 1,902 |
| Louvain modules | 12 (Q = 0.3532, frozen run — see the note below) |
| Prioritised genera | 12 |
| TSWV sequences / VC80 clusters / countries | 1,260 / 29 / 30 |

> The 29 here are **VC80 clusters**, not species. The TSWV species-level subsets
> in `Table_S5` are counted in isolates, and the five gene datasets are counted
> in sequences; `data/README.md` defines every unit used in the release.

> **On the 12 modules.** The genus–host network is partitioned with Louvain
> modularity, a randomised multi-start heuristic. Four facts bound what this
> release claims:
>
> 1. Fig. 4 has been re-rendered from the frozen result tables and passes
>    (12 modules, Q = 0.3532).
> 2. Node metrics (degree / strength / betweenness) and the set of 12
>    prioritised genera are identical between the frozen run and a re-run.
> 3. The Louvain community partition itself is **not** reproduced value-for-value:
>    a re-run of `05A_network_main.R` on the same edges and seed lands on a
>    different, equally valid optimum (typically 10–12 modules). The mechanism
>    has not been pinned down and no single cause is claimed here.
> 4. The published module tables are the frozen partition the paper uses,
>    MD5-matched to the analysis run. Details in
>    `FINAL_RELEASE_VALIDATION.md` §B4.

Analysis units are defined once, in `data/README.md`; the distinction between a
**sequence**, a **VC80 unit**, a **virus TaxID**, and an **RdRP65 operational
lineage** matters throughout and is not a naming convention.

---

## Where to start

| Goal | Start here |
|---|---|
| Explore the released datasets | `data/README.md`, `data/curated/` |
| Trace one manuscript number back to its script | `docs/MANUSCRIPT_OUTPUT_MAP.tsv` |
| Re-run an analysis module | `analysis/README.md` |
| Understand the three reproducibility levels | `docs/reproducibility.md` |
| Check what must be obtained from outside | `docs/EXTERNAL_INPUTS.md` |
| See what was checked, and what was not | `FINAL_RELEASE_VALIDATION.md`, `UNRESOLVED_ITEMS.md` |
| See how this release differs from the previous one | `RELEASE_DIFF_AUDIT.tsv` |

---

## Quick start

```bash
git clone https://github.com/xinrui-mi/plant-virus-host-range.git
cd plant-virus-host-range

source analysis/00_config/00_paths.sh      # exports ANALYSIS_ROOT, REPO_ROOT, DATA_DIR, ...
bash   analysis/00_config/00_preflight.sh  # verifies every released input and required column
sha256sum -c checksums/SHA256SUMS.txt      # verifies file integrity

Rscript analysis/01_dataset_summary/analysis/01_Fig1_data_framework_QC.R
```

No additional external data are required for the released analysis chain.
`EXTERNAL_DATA_DIR` is only needed to rebuild the family phylogenies, or to
re-run the TSWV gene-extraction / GARD / selection workflow from third-party
sequence resources.

Environment: **R 4.4.3** and **Python 3.12**, with the exact package versions
recorded under `environment/`.

---

## The six analysis modules

| Module | Manuscript location | What it produces |
|---|---|---|
| `01_dataset_summary` | §3.1 · Fig. 1 · Table S1 | dataset integration, taxonomy/metadata coverage, analysis-unit coverage, VC80 internal consistency, RNA evidence space, 100-permutation evidence accumulation |
| `02_host_relatedness_and_phylogeny` | §3.2–3.3 · Fig. 2A–C, 2D · Figs. S1–S2 · Table S2 | host distance decay 0–6, degree-preserving null, sampling and categorical-distance sensitivity, 24 family-stratified slopes, **12** family-specific marker trees and host-spectrum correlations |
| `03_arthropod_transmission_ecology` | §3.4 · Fig. 3 · Fig. S3 · Table S3 | S0/S1/S2 evidence tiers, adjusted host-range models, four transmission modes, five-group contrasts, country-breadth models, within-mode lineage-composition sensitivity |
| `04_genus_host_network` | §3.5 · Fig. 4 · Fig. S4 · Table S4 | genus–host-family network, module detection, bridge and candidate scores, 200× resampling, threshold perturbation, small-world randomisation, cross-genus VC80 sensitivity |
| `05_vector_contextualization` | §3.5 · Fig. 4 · Table S4 | curated external evidence and support-adjusted VIP-DB contextualisation of the prioritised genera |
| `06_TSWV_evolution` | §3.6 · Fig. 5 · Table S5 | TSWV membership gate, five-gene phylogenies, study-aware segment pairing, segment concordance, GARD recombination and selection |

Shared configuration, helper functions and the preflight check live in
`analysis/00_config/`.

---

## Repository layout

```
plant-virus-host-range/
├── README.md
├── FINAL_RELEASE_VALIDATION.md     what was verified, module by module
├── UNRESOLVED_ITEMS.md             what is still open, and why
├── RELEASE_DIFF_AUDIT.tsv          old path -> final path, decision, change type
├── LICENSE                         code: MIT
├── LICENSE_DATA                    project-produced data: CC BY 4.0
├── CITATION.cff
│
├── data/
│   ├── README.md                   dataset descriptions and key fields
│   ├── curated/                    released curated datasets
│   ├── analysis_input/             analysis-ready inputs read by the code
│   ├── analysis_data/
│   │   └── phylogeny/
│   │       ├── twelve_family/      12 marker trees + AA95 membership + tip-host tables
│   │       └── TSWV/               L / M / S segment trees + matched metadata
│   └── external_evidence/          curated external evidence + provenance
│
├── analysis/
│   ├── README.md                   module descriptions and execution order
│   ├── 00_config/                  configuration, helpers, preflight
│   ├── 01_dataset_summary/
│   ├── 02_host_relatedness_and_phylogeny/
│   ├── 03_arthropod_transmission_ecology/
│   ├── 04_genus_host_network/
│   ├── 05_vector_contextualization/
│   └── 06_TSWV_evolution/
│
├── additional_validation/          NOT_IN_PAPER — see its README
├── supplementary/                  Data S1 + machine-readable Tables S1–S5
├── environment/                    software versions and package information
├── docs/
│   ├── MANUSCRIPT_OUTPUT_MAP.tsv
│   ├── FROZEN_DEPENDENCY_CLOSURE.tsv
│   ├── EXTERNAL_INPUTS.md
│   └── reproducibility.md
└── checksums/
    └── SHA256SUMS.txt
```

---

## Reproducibility, in one paragraph

There are three levels, and they are not interchangeable. **L2** (recommended)
re-runs the published analyses from the released analysis-ready data and needs
nothing outside this repository. **L3** re-renders the figures from the frozen
result tables. **L1** rebuilds the database from raw NCBI Virus / Virus-Host DB /
Land Flora / VIP-DB downloads; it is documented but explicitly *not* claimed as a
one-command reproduction, because the upstream harmonisation was interactive and
embeds manual curation decisions. Full detail, including how to run a fast
reduced-iteration pass, is in `docs/reproducibility.md`.

---

## What is deliberately *not* in this repository

- The raw third-party database dumps (NCBI Virus, Virus-Host DB, Land Flora RNA
  virome, VIP-DB workbook, Pfam, GenBank flat files). See `docs/EXTERNAL_INPUTS.md`.
- **Rendered figure files.** This release contains the datasets, the analysis
  code and the plotting code — no finished figure images (PDF/PNG/TIFF). Every
  figure can be re-rendered from the frozen result tables: Figs. 2–5 with the
  shipped `run_*.sh` wrappers, Fig. 1 by calling its plotting script directly
  (it has no wrapper). See `analysis/PLOT_MAP.md`.
- The 884 intermediate scripts from the project's history. Only the closed
  dependency set that produces a released result is shipped; the selection rule
  and the per-script trace are in `docs/FROZEN_DEPENDENCY_CLOSURE.tsv`.
- Half-finished plot attempts. One plotting directory per manuscript figure is
  shipped.
- IQ-TREE intermediate files (`.mldist`, `.ckp.gz`, `.splits.nex`, …). The
  trees, the model reports, the filtered alignments and the tip-to-host tables
  are shipped; the search scratch files are not.

---

## Licence

Code is released under the **MIT License** (`LICENSE`). Project-generated derived
tables, curated evidence, supplementary tables and the compilation layer of the
released datasets are released under **CC BY 4.0** (`LICENSE_DATA`).

`LICENSE_DATA` covers what this project produced. It does **not** re-licence any
third-party source data, and nothing here should be read as claiming that NCBI
Virus, Virus-Host DB, Land Flora RNA virome or VIP-DB content is available under
CC BY 4.0. See `docs/EXTERNAL_INPUTS.md` §5.

## How to cite

See `CITATION.cff`, and cite the manuscript itself. Please also cite the upstream
resources listed in `docs/EXTERNAL_INPUTS.md` §6.

## Contact

For reproducibility questions or discrepancies, please open a GitHub issue. For
manuscript correspondence, see the paper.
