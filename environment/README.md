# Environment

Computational environment in which the manuscript analyses were run.

| File | Contents |
|---|---|
| `software_versions.txt` | Tool versions captured at freeze time (2026-08-14). Personal conda prefixes replaced by `<CONDA_PREFIX>`. |
| `R_sessionInfo.txt` | Full `sessionInfo()` from the analysis environment (R 4.4.3, CentOS 7). Same placeholder substitution. |
| `R_packages_used.txt` | Every R package referenced by `analysis/`, with the version actually installed — derived by static scanning of `library()` / `require()` / `::` calls. |
| `python_version.txt` | Python interpreter version record. |

The executables below (CD-HIT, MAFFT, trimAl, IQ-TREE, HMMER, HyPhy) are needed only by the optional
`pipeline/` scripts that rebuild trees; the released L2 chain runs without them.

## Headline versions

| Tool | Version |
|---|---|
| R | 4.4.3 (2025-02-28) |
| Python | 3.12.12 |
| CD-HIT | 4.8.1 |
| MAFFT | 7.526 |
| trimAl | 1.5.rev1 |
| IQ-TREE | 3.0.1 |
| HMMER | 3.4 |
| HyPhy | 2.5.101 (MP) |

Key R packages: `ape` 5.8-1, `data.table` 1.18.4, `igraph` 2.1.4, `MASS` 7.3-66, `Matrix` 1.7-6,
`readxl` 1.5.0. Plot code additionally uses `ggplot2` 4.0.3, `ggtree` 3.14.0, `patchwork` 1.3.2,
`sf` 1.0-20, `rnaturalearth` 1.2.0, `ragg` 1.5.2. Full list in `R_packages_used.txt`.

## Reproducing the environment

The L2 analysis chain needs only the analysis packages:

```bash
conda create -n pvh -c conda-forge \
  r-base=4.4.3 r-data.table r-ape r-igraph r-mass r-matrix r-readxl \
  python=3.12 pandas numpy biopython
```

The L3 figure chain additionally needs the plotting packages:

```bash
conda install -c conda-forge \
  r-ggplot2 r-ggtree r-patchwork r-ggrepel r-sf r-rnaturalearth r-scales r-ragg \
  r-stringr r-dplyr r-tidyr r-tibble r-readr r-systemfonts r-grdevices
```

The optional tree-rebuilding pipelines additionally need:

```bash
conda install -c bioconda cd-hit mafft trimal iqtree hmmer hyphy
```

Pin versions before publishing results — several packages (`dplyr`, `ggplot2`, `igraph`) have active
release cycles and defaults do change.
