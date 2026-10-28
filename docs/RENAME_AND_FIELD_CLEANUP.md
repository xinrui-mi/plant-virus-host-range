# Rename and public-field cleanup record

A short log of every file renamed and every column dropped or reworded for the
public release, with the check that made each one safe. Nothing here changes an
analysis result.

## Renamed

| Before | After | Why safe |
|---|---|---|
| `data/external_evidence/prioritized_genera_external_evidence_FINAL.tsv` | `data/external_evidence/prioritized_genera_evidence_summary.tsv` | Reader-facing summary. `06A_candidate_external_context.R` does not read it (it writes `tables/02_candidate_external_evidence_matrix.tsv` itself), so no script breaks. References updated in `docs/MANUSCRIPT_OUTPUT_MAP.tsv` and `data/README.md`. |
| `data/external_evidence/priority_genus_external_evidence_provenance_FINAL.tsv` | `data/external_evidence/prioritized_genera_evidence_sources.tsv` | Same reasoning; also the old name was internally inconsistent (`priority_genus` vs `prioritized_genera`). Content changed too — see below. |

Both files lost their `FINAL` suffix because a public reader has no way to know
what it distinguishes.

## Not renamed, deliberately

| File | Reason |
|---|---|
| `data/analysis_data/20_three_source_sequence_master_FINAL.tsv.gz` | The analysis scripts reference this exact path. Renaming it would break `00_config.R` and every downstream read. |
| `data/analysis_data/20_three_source_association_master_FINAL.tsv.gz` | Same. |
| `data/analysis_data/20_arthropod_sequence_evidence_FINAL.tsv.gz` | Same. |
| `analysis/06_TSWV_evolution/plot/fig5/07Q_Fig7_…R` | Frozen artwork script; the name records which cluster run produced the submitted figure. |
| All `*_FINAL_v*/run_*.sh` wrappers | Version suffixes identify the frozen runs. Renaming to match figure numbers would detach the code from its provenance. |
| `analysis/*/pipeline_results/03L_FINAL_12family_host_spectrum/<Family>/…` | Read by `03L_FINAL_aggregate.R` and Fig. 2 by literal file name. |

A blanket `20_` → `21_` substitution was explicitly **not** performed: the `20_`
tables are the analysis masters the frozen code reads, and their provenance
columns were verified against the released `data/curated/` tables rather than
inferred from the file name.

## Columns dropped

`prioritized_genera_evidence_sources.tsv` (was the provenance file), 28 rows,
11 columns → 9:

| Column | Action | Justification |
|---|---|---|
| `Archived evidence` | **dropped** | Held only an internal Word document name (`eppo补充.docx`) — not a public source, not redistributable. |
| `Archive page(s)` | **dropped** | Page ranges into that internal document; meaningless to a reader without it. |

Every one of the 28 rows was checked **before** dropping: all retain a non-empty
`URL`, so no row loses its point of origin. Columns kept:
`Genus`, `Family`, `Named member / scope`, `Claim category`, `Claim supported`,
`Source type`, `Source title / identifier`, `URL`, `Verification status`.

The verification date survives inside `Verification status`
(`Verified 2026-09-10` on all 28 rows), so no separate date column was needed.

## Columns deliberately kept

| File / column | Kept because |
|---|---|
| `candidate_external_evidence.tsv.gz` → `external_annotation_source`, `external_notes` | These name the *type* of external annotation (EPPO / IPPC / literature / VIP-DB), not a private file. No internal path or personal document is referenced. |
| `…_note` columns (11 of them) | Per-claim free text that carries the actual reading of each source; removing it would strip the evidence. |

## Wording fixed, not renamed

| Location | Before | After |
|---|---|---|
| `analysis/05_vector_contextualization/tables/02_candidate_external_evidence_matrix.tsv` | "not exact Fig.6 accession-level (3)" | "not exact Fig. 4 accession-level (3)" |
| `analysis/05_vector_contextualization/analysis/06B_…R` (2 comments) | `Fig.6B`, `Fig.6` | `Fig. 4B`, `Fig. 4` |
| `README.md` | "TSWV sequences / **species** / countries  1,260 / 29 / 30" | "TSWV sequences / **VC80 clusters** / countries", with a note on which unit S5 and the gene datasets use |

## Verification after all of the above

- Every script that reads a renamed file was checked to not exist (renames were
  only applied where nothing reads the name).
- `grep -r "eppo补充\|Archived evidence\|Archive page"` over the tree → 0 hits.
- All eight plot wrappers and `03L_FINAL_aggregate.R` executed from a fresh
  clone: Fig. 1–5 pass; Fig. S1 fails on two never-shipped inputs, recorded as
  U22 in `UNRESOLVED_ITEMS.md`.
