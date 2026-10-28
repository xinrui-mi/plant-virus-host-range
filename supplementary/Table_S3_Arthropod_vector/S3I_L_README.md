# Table S3 — sections I to L

Sections A–H of Table S3 are shipped as one TSV per section (see the table in
`supplementary/README.md`). Sections I–L are the part of Table S3 that comes from
`04D_transmission_biology_extended_v3.R`, which writes 19 result tables
(`20_04D_*` … `38_04D_*`) into `analysis/03_arthropod_transmission_ecology/tables/`.

Those 19 tables have different column sets (9 to 41 columns), so each section is shipped as a
**stacked table with one self-contained block per source table**. Every block repeats two leading
columns and its own header row:

| Column | Meaning |
|---|---|
| `section_block` | Short label of the block, so rows can be read without consulting this file |
| `source_table` | Exact file name of the `04D` output the row came from |

All remaining columns are the original columns of that source table, copied verbatim. **No value
is recomputed, rounded, reformatted or padded.** Empty trailing cells in a block are the normal
consequence of `04D` writing blocks of unequal width; they are not missing data.

Rebuild with:

```bash
python3 supplementary/Table_S3_Arthropod_vector/build_S3I_L.py
```

The script only re-stacks; it fails loudly if a source table is absent.

---

## S3I — vector family × virus genus

| Block | Source | Rows | Figure |
|---|---|---|---|
| `vector_family_x_virus_genus` | `24_04D_vector_family_by_virus_genus.tsv` | 103 | context for Fig. 3D (not plotted directly) |

Cross-tabulation of vector family against virus genus with the canonical four-mode label and the
number of virus TaxIDs behind each combination. This is the table the manuscript's §3.4 statement
on lineage-specific transmission aggregation rests on.

## S3J — five-group host-range models

| Block | Source | Rows | Figure |
|---|---|---|---|
| `descriptive_summary` | `20_04D_five_group_descriptive_summary.tsv` | 5 | Fig. 3D |
| `host_range_model_terms` | `25_04D_five_group_host_range_model_terms.tsv` | 270 | Fig. 3D |
| `host_range_primary_contrasts` | `26_04D_five_group_host_range_primary_contrasts.tsv` | 36 | Fig. 3D |
| `lineage_leave_one_out` | `27_04D_lineage_leave_one_out_primary_contrasts.tsv` | 4,140 | — |
| `lineage_driver_ranking` | `34_04D_lineage_driver_ranked_summary.tsv` | 240 | — |

The five transmission groups are *no documented vector relationship*, *non-persistent*,
*semi-persistent*, *persistent non-propagative* and *persistent propagative*. The descriptive block
gives per-group host breadth; the model blocks give the adjusted host-range terms and the primary
contrasts between groups; the leave-one-lineage-out blocks test whether any single lineage drives
the contrast.

## S3K — within-mode lineage sensitivity

| Block | Source | Rows | Figure |
|---|---|---|---|
| `within_family_mode_overlap_audit` | `23_04D_within_family_mode_overlap_audit.tsv` | 74 | Fig. S3B |
| `leave_one_lineage_out_eligibility` | `35_04D_leave_one_lineage_out_eligibility_audit.tsv` | 345 | — |
| `within_mode_leaveout_contrasts` | `36_04D_within_mode_lineage_leaveout_primary_contrasts.tsv` | 207 | Fig. S3B |
| `within_mode_influence_ranking` | `37_04D_within_mode_lineage_influence_ranked_summary.tsv` | 180 | — |
| `within_mode_eligibility` | `38_04D_within_mode_lineage_eligibility_audit.tsv` | 69 | — |

This is the answer to the obvious objection to S3J: the five transmission groups differ in which
virus families and genera sit inside them, so a group contrast could be a lineage contrast. The
blocks audit how much mode overlap exists within a family, then remove one lineage at a time and
re-fit, and record for each removal whether the five groups were still all populated. Rows where
`all_five_groups_retained` is false are reported, not dropped.

## S3L — country breadth and evidence-support models

| Block | Source | Rows | Figure |
|---|---|---|---|
| `country_breadth_model_terms` | `28_04D_country_breadth_model_terms.tsv` | 249 | Fig. 3D |
| `country_breadth_primary_contrasts` | `29_04D_country_breadth_primary_contrasts.tsv` | 20 | Fig. 3D |
| `S0S2_vs_VIPDB_raw_fisher` | `30_04D_S0S2_vs_VIPDB_raw_fisher.tsv` | 1 | Fig. 3B |
| `S0S2_vs_VIPDB_adjusted_terms` | `31_04D_S0S2_vs_VIPDB_adjusted_model_terms.tsv` | 76 | Fig. 3B |
| `S0S2_vs_VIPDB_adjusted_primary` | `32_04D_S0S2_vs_VIPDB_adjusted_primary.tsv` | 1 | Fig. 3B |
| `binary_evidence_country_breadth` | `33_04D_binary_evidence_country_breadth_primary.tsv` | 4 | Fig. 3D |

These blocks quantify the association between arthropod-associated sequence evidence and
documented vector relationships, before and after adjustment for VC80 support and virus family
(raw Fisher exact, then the adjusted model terms and primary contrast). The country-breadth
blocks model geographic range against transmission group, with the M0–M4 adjustment ladder in the
term column.

---

## Verification performed

All 17 source-table blocks were re-read and compared cell by cell against the `04D` output:
**every data cell matches, all 17 blocks, 0 differences**. Row counts:

| Section | Blocks | Data rows |
|---|---|---|
| S3I | 1 | 103 |
| S3J | 5 | 4,691 |
| S3K | 5 | 875 |
| S3L | 6 | 351 |

## Caveat that applies to the whole of Table S3

The header rows are per block, so a naive `read.table(header=TRUE)` on any of these four files will
produce a single confusing frame. Read with `skip`/`nrows` per block, or filter on `section_block`
after reading with `header=FALSE` and assigning the correct header per block. This is a deliberate
trade-off: one wide rectangle would have required padding 41-column blocks and 9-column blocks into
a shared schema, and a row of empty cells is a worse lie than a repeated header.
