# Transmission-lineage specificity control — `NOT_IN_PAPER`

This is the specificity control run after the first transmission-lineage model.
It exists because that model showed near-perfect prediction for vector x
virus-lineage combinations; before claiming that "the same insect" carries
predictive information we must test whether virus genus or family **alone**
already predicts the transmission mode equally well.

## Design

Same virus-vector pair targets, eight predictor schemes:

| # | scheme | columns used |
|---|---|---|
| 1 | virus genus only | `virus_genus` |
| 2 | virus family only | `virus_family` |
| 3 | vector species only | `vector_species_id` |
| 4 | vector genus only | `vector_genus` |
| 5 | vector species x virus genus | `vector_species_id`, `virus_genus` |
| 6 | vector genus x virus genus | `vector_genus`, `virus_genus` |
| 7 | vector species x virus family | `vector_species_id`, `virus_family` |
| 8 | vector genus x virus family | `vector_genus`, `virus_family` |

Each scheme is evaluated at three settings: `min_train` = 2 training viruses,
with dominant-mode purity thresholds 0.0 / 0.8 / 1.0.

**Primary strict rule** — at least 2 other virus TaxIDs in the training group and
100% dominant-mode purity. A target is only scored when its predictor columns are
all non-missing *and* the training group resolves to a single winning mode.

Pairs whose own records disagree on the mode are dropped first and listed in
`00_ambiguous_pairs.tsv`.

## Outputs

| file | rows | content |
|---|---|---|
| `00_ambiguous_pairs.tsv` | 1 (header only) | pairs excluded because they carry >1 mode |
| `01_predictions_all_schemes.tsv` | 13,104 | one row per target x scheme x setting |
| `02_specificity_control_summary.tsv` | 24 | coverage / accuracy / balanced accuracy |
| `03_recall_by_mode.tsv` | 96 | recall per true mode |
| `04_incremental_vector_value.tsv` | 4 | accuracy on the shared target set, baseline vs conditioned |
| `05_genus_global_mode_consistency.tsv` | 46 | per-genus dominant mode and consistency |
| `06_family_global_mode_consistency.tsv` | 18 | per-family dominant mode and consistency |

## Result

Strict rule, from `04_incremental_vector_value.tsv`:

| baseline | conditioned | common n | baseline acc. | conditioned acc. |
|---|---|---|---|---|
| virus genus only | vector species x virus genus | 365 | 1.000 | 1.000 |
| virus genus only | vector genus x virus genus | 375 | 0.997 | 0.997 |
| virus family only | vector species x virus family | 347 | 1.000 | 1.000 |
| virus family only | vector genus x virus family | 356 | 0.997 | 0.997 |

Virus lineage alone already predicts the four-class mode almost perfectly.
Conditioning on vector identity reduces coverage (0.84-0.87 -> 0.67-0.73) without
raising accuracy. This supports presenting the Results 3.4 aggregation as a
**lineage-level** property rather than as evidence that vector identity carries
extra information.

## Reproduction status

These seven tables were regenerated in this repository from
`data/analysis_data/VIPDB_relation_normalized.tsv.gz` and match the frozen
cluster run **row for row (13,104 + 24 + 96 + 4 + 46 + 18 + 0 rows, identical
values)**. Both runs use 546 eligible pairs, 414 virus TaxIDs.

Differences from the cluster run are limited to serialization: `fwrite` field
quoting and, in an earlier pass, label case. `fwrite(..., quote = FALSE)` and
`tolower()` are now pinned in the script so the output is byte-stable across
`data.table` versions.
