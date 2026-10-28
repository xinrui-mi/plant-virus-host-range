# Supplementary tables

This directory holds the supplementary tables of the manuscript, plus the
single-workbook form `Data_S1.xlsx`.

| Table | Files | Manuscript |
|---|---|---|
| `Table_S1_Data_and_units/` | 4 | Results 3.1 · Fig. 1 — source composition, coverage, VC80 consistency, RdRP threshold calibration |
| `Table_S2_Host_relatedness/` | 10 | Results 3.2 · Fig. 2A–C (`S2A1`–`S2D`); Results 3.3 · Fig. 2D · Fig. S2 (`S2E1`–`S2E5`); see the significance-count note below for 11/12 vs 9/12 vs 12/12. **Fig. S1** is the single-panel RdRP65 host-sharing decay plot (its source data is the frozen `04_RdRP65_host_sharing_by_distance` table in the analysis module, same evidence layer as `S2A1`–`S2A2`) |
| `Table_S3_Arthropod_vector/` | 15 TSV + 1 README + 1 rebuild script | Results 3.4 · Fig. 3 · Fig. S3 — sections A through L |
| `Table_S4_Prioritized_genera/` | 5 | Results 3.5 · Fig. 4 · Fig. S4 |
| `Table_S5_TSWV/` | 4 | Results 3.6 · Fig. 5 — coverage, five-gene datasets, segment concordance, recombination and selection |


## Reading Table S2E: three significance counts that must not be conflated

The twelve-family host-spectrum analysis reports **three different counts**, and
they answer three different questions. Do not quote one in place of another:

| Analysis | File | Column | Families significant |
|---|---|---|---|
| Main test (one-sided permutation, BH-corrected) | `S2E1_twelve_family_phylogenetic_summary.tsv` | `BH_q_greater` | **11 of 12** (`Mitoviridae` q = 0.095) |
| Sequence-support-stratified sensitivity | `S2E2_support_stratified_permutation.tsv` | `BH_q_greater` | **9 of 12** |
| Minimum-support sensitivity, each tip needs ≥2 sequences | `../../analysis/02_host_relatedness_and_phylogeny/pipeline_results/03L_FINAL_12family_host_spectrum/04_FINAL_min_support_sensitivity_summary.tsv` | `BH_q_greater` | **12 of 12** (10 of 12 at ≥3) |

`S2E4_high_support_sensitivity.tsv` is the minimum-support sensitivity
analysis, exported byte-for-byte from the frozen
`analysis/02_host_relatedness_and_phylogeny/pipeline_results/03L_FINAL_12family_host_spectrum/04_FINAL_min_support_sensitivity_summary.tsv`
(12 families × 2 thresholds). It matches `Data_S1.xlsx` (Table S2) and the
manuscript numbers; a former 4-family export from the pre-`03L` run was
replaced (see `UNRESOLVED_ITEMS.md` U24).

S2E1–E4 are produced by the `03L` chain (`03L_FINAL_run_one_family.R` ×12 →
`03L_FINAL_aggregate.R`). `03G_phylogeny_host_spectrum_sensitivity.R` is an
earlier independent implementation kept for inspection; it is **not** the source
of any published S2E number.


## State of this directory

The tables that the manuscript actually reports are the frozen submission
outputs. They were **not** recomputed during release preparation. Point-wise
provenance (which script, which input, which reported number) is in
`../docs/MANUSCRIPT_OUTPUT_MAP.tsv`.

Five groups of files differ from the previously published GitHub snapshot, and
their status is stated explicitly:

| File | Status | Note |
|---|---|---|
| `Table_S2_Host_relatedness/S2E1_twelve_family_phylogenetic_summary.tsv` | `FINAL` | Replaces `S2E1_four_family_phylogenetic_summary.tsv`. The final analysis uses **12** family-specific marker trees, not 4. The old four-family file was **removed** rather than left beside the new one, so it cannot be picked up by mistake. Rebuilt by `03L_FINAL_run_one_family.R` (×12) + `03L_FINAL_aggregate.R`; 12 families, 6,010 tips. |
| `Table_S2_Host_relatedness/S2E2_support_stratified_permutation.tsv` | `FINAL` | The previous release had E2 and E3 in the reverse order of the manuscript. The final numbering is **E2 = support-stratified permutation, E3 = pair-sampling stability**. Both files were regenerated in the corrected order. |
| `Table_S2_Host_relatedness/S2E3_pair_sampling_stability.tsv` | `FINAL` | See above. |
| `Table_S3_Arthropod_vector/S3I..S3L_*.tsv` | `FINAL` | The previous release stopped at section H. Sections I–L are stacked views of the 19 tables written by `04D_transmission_biology_extended_v3.R`; see `S3I_L_README.md` in that directory for the block map and the reading caveat. All 17 source blocks were compared cell by cell against the `04D` output: 0 differences. |
| `Data_S1.xlsx` | `FINAL` | One corrected cell only: `Contents & Notes!A18` now reads `4,401 plant-virus sequences` instead of `4,401 TaxIDs`, and the S0/S1/S2 definitions carry their counting units. A full cell-by-cell comparison against the submitted workbook reports **exactly 1 changed cell**; `Table S1`–`Table S5` are unchanged (0 changed cells across 46×5, 151×20, 6081×41, 348×15 and 36×9). |


`Data_S1.xlsx` is the single-workbook form of Tables S1–S5, included for
reviewers who prefer to page through a workbook.

## Reading Table S3, sections I–L

`S3I`, `S3J`, `S3K` and `S3L` are **not** ordinary rectangular tables. Each file
stacks several source tables one below the other, and **every block carries its
own header row**:

| File | Blocks | Data rows | Header widths |
|---|---:|---:|---|
| `S3I_vector_family_by_virus_genus.tsv` | 1 | 103 | 7 |
| `S3J_five_group_host_range_models.tsv` | 5 | 4,691 | 16, 19, 20, 43, 43 |
| `S3K_within_mode_lineage_sensitivity.tsv` | 5 | 875 | 6, 19, 43, 43, 12 |
| `S3L_country_breadth_and_evidence_models.tsv` | 6 | 351 | 19, 20, 11, 19, 20, 20 |

Section letters follow `Data_S1.xlsx` (Table S3): **K = within-mode lineage
leave-out sensitivity, L = country breadth and evidence models**. A previous
release had these two file names swapped relative to the Excel and manuscript
labelling; they were renamed on 2026-10-27 (see `UNRESOLVED_ITEMS.md` U25), the
builder `build_S3I_L.py` was updated with them, and the renamed files were
verified byte-identical (MD5) to the previous contents.

Because the widths differ, a plain read gives you a wrong table — silently:

```r
x <- read.delim("S3J_five_group_host_range_models.tsv")
# 4695 rows x 43 columns. Nothing errors and the row count is right —
# but the column names come from block 1 (16 columns) and blocks 2-5
# were padded with NA into that 43-column frame.
```

Read line by line and split on the block header instead:

```r
lines  <- readLines("S3J_five_group_host_range_models.tsv")
starts <- grep("^section_block\t", lines)          # one header per block
ends   <- c(starts[-1] - 1, length(lines))

blocks <- Map(function(a, b) {
  d <- read.delim(text = lines[a:b], check.names = FALSE, fill = TRUE)
  # fill = TRUE pads short blocks with unnamed columns; drop them.
  d[, !(names(d) == "" | is.na(names(d))), drop = FALSE]
}, starts, ends)

names(blocks) <- vapply(blocks, function(d) d$section_block[1], "")

# descriptive_summary            5 x 16
# host_range_model_terms       270 x 19
# host_range_primary_contrasts  36 x 20
# lineage_leave_one_out      4,140 x 43
# lineage_driver_ranking       240 x 43
```

`check.names = FALSE` matters: the source tables share column names such as
`term`, `estimate` and `p_value` across blocks, and default name repair would
rewrite them to `term.1`, `term.2`, … and break any code that joins on them.
`fill = TRUE` is what lets the 16-column first block sit in the same read loop
as the 43-column ones; the trailing unnamed columns it creates are the padding
and carry nothing.

If you would rather have the 17 source tables as separate files,
`build_S3I_L.py` in that directory regenerates them from the `04D` output, and
`S3I_L_README.md` maps each block back to its source file and row count.


> **`Data_S1.xlsx` is authoritative for the numbers.** It is the workbook the
> journal reviewed. The TSVs beside it are a convenience re-export: Table S1,
> S4 and S5 are direct copies, Table S2 and sections S3A-S3H are the frozen
> submission outputs, and sections S3I-S3L are regenerated from the released
> `04D` output and verified byte-identical. Where a TSV and the workbook
> disagree, take the workbook, and please tell us - a disagreement would be a
> defect rather than a rounding choice.

---

## Land Flora evidence tiers (fields `landflora_*`) — read A/B correctly

`data/curated/virus_host_association_dataset.tsv.gz` carries the Land Flora
records (Liu et al. 2026, *Molecular Plant* 19(10), doi:10.1016/j.molp.2026.07.023;
public dataset doi:10.6084/m9.figshare.30688322) with four annotation fields:
`landflora_record_id`, `landflora_evidence_tier`, `landflora_strict_plant_infection`
and `landflora_host_resolution_rank`.

**How the two groups were formed.** Land Flora records were recoded into two
groups using the original study's final lineage classification in **Data S4**.
Group A comprises 841 records whose classification labels contain
**"Plant-infecting"**, including labels qualified by **"for MP+ lineage"**;
group B comprises the remaining 1,675 included plant-associated records. These
groups represent a recoding of the source classification and do **not** establish
experimentally confirmed infection for every virus-plant association. The
sensitivity analysis retains only group A from Land Flora; records from the
other two sources remain unchanged.

**Sensitivity-analysis naming — one equivalence table:**

| Analysis | Name in Table S3H (`S3H_RdRP65_sensitivity.tsv`) and `Data_S1.xlsx` | Name in the frozen 04B output (`16_RdRP65_sensitivity_summary.tsv`) |
|---|---|---|
| Main analysis (groups A+B) | Land Flora A+B evidence set | `main_LF_AplusB` |
| Strict sensitivity (group A only) | Land Flora group-A-only sensitivity set | `strict_LF_A_only` |

`strict_LF_A_only` is a **historical internal name** retained in the frozen
output for provenance; it does not imply that every group-A record was
individually confirmed as a direct infection. The reader-facing name is
"Land Flora group-A-only sensitivity set" (renamed 2026-10-27; formerly
"Land Flora A-only direct-infection evidence set", which overstated what the
recoding supports). Group membership itself is unchanged — only the label was
corrected. Field values in the curated tables (`A_plant_infection_supported`,
`B_plant_associated_not_direct_infection`, `landflora_strict_plant_infection`)
are likewise kept verbatim as historical field names.
