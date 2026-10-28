# `analysis/` — manuscript analysis modules

Six modules, one per manuscript result block. Each has the same internal layout:

```
analysis/<module>/
├── analysis/      # R / Python scripts that produce the result tables
├── pipeline/      # (only where needed) scripts that rebuild upstream artefacts
└── plot/          # (level-3 figures) scripts that render the submitted artwork
```

Scripts write into `<module>/tables/`, `<module>/qc/`, `<module>/logs/` and
`<module>/work/`; those directories are created on demand and are ignored by git.

| module | manuscript | results produced |
|---|---|---|
| `00_config/` | — | shared configuration, path resolution, helpers, preflight |
| `01_dataset_summary/` | Results 3.1 · Fig. 1 | source composition, taxonomy/metadata coverage, analysis-unit coverage, VC80 internal consistency, RNA evidence space, 100-permutation evidence accumulation |
| `02_host_relatedness_and_phylogeny/` | Results 3.2–3.3 · Fig. 2 · Figs. S1–S2 | host taxonomic distance decay (0–6), degree-preserving null, sensitivity, 24 family-stratified slopes, 12 family-specific marker trees and host-spectrum correlations |
| `03_arthropod_transmission_ecology/` | Result 3.4 · Fig. 3 · Fig. S3 | S0/S1/S2 arthropod evidence tiers, support- and family-adjusted host-range models, four transmission modes, five-group contrasts, country-breadth models, within-mode lineage-composition sensitivity |
| `04_genus_host_network/` | Results 3.5 · Fig. 4 · Fig. S4 | 269-genus × 197-host-family weighted network, module detection, bridge/candidate scores, 200× resampling robustness, threshold perturbation, small-world randomisation, cross-genus VC80 sensitivity |
| `05_vector_contextualization/` | Results 3.5 · Fig. 4 | curated external evidence and support-adjusted VIP-DB contextualization of the prioritised genera |
| `06_TSWV_evolution/` | Results 3.6 · Fig. 5 | TSWV membership gate, five-gene phylogenies, study-aware segment pairing, segment concordance, GARD recombination and selection scans |

## Running order

```bash
source analysis/00_config/00_paths.sh
bash   analysis/00_config/00_preflight.sh
```

Then run the modules in numeric order. The only hard cross-module dependencies are:

1. `01` → `02` (`02` reads `01_dataset_summary/tables/01_VC80_internal_consistency_QC.tsv`)
2. `03/04A` → `03/04B,04C,04D` (normalized VIPDB)
3. `04/05A` → `04/05B-05E` and → `05/06A,06B` (candidate thresholds)
4. `03L_FINAL_run_one_family.R` (×12) → `03L_FINAL_aggregate.R`

## Where the numbers come from

`../docs/FROZEN_DEPENDENCY_CLOSURE.tsv` lists, script by script, the frozen
inputs, the frozen outputs and the manuscript anchor, together with the evidence
used to establish that trace. `../docs/MANUSCRIPT_OUTPUT_MAP.tsv` gives the
reverse mapping: manuscript item → file → script → input.

## A note on the numbering

The script prefixes (`03A`, `03K`, `04D`, `07Q`, …) are the original
project-internal step numbers and are **not** manuscript figure numbers. The
folder names are what map to the manuscript; the prefixes are kept so that the
released code stays traceable to the frozen cluster run.
