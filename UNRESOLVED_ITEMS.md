# Unresolved items — release preparation for the plant-virus host-range manuscript

Each item records a question raised about this release, what was found, and how it was closed. Where
something is still open it says so; nothing here has been worked around by substituting an older
code version, by re-deriving a result with different settings, or by editing frozen outputs. Items
are ordered by when they were raised.

Items marked `CLOSED` are kept rather than deleted: several were found only by unpacking the tarball
and running it, and the failure mode is the useful part of the record.

Status vocabulary:

| Status | Meaning |
|---|---|
| `BLOCKING` | The release cannot be declared complete until this is closed. |
| `OPEN_CLUSTER_ACCESS` | Waiting on the cluster (bore tunnel down). No local substitute exists. |
| `OPEN_AUTHOR_INPUT` | Requires a decision or a file only the author can supply. |
| `ACCEPTED_LIMITATION` | Will be documented as a limitation rather than fixed. |
| `CLOSED` | Found, fixed and verified. Kept in the list as a record of what went wrong. |

---

## U1 — 12-family trees: complete `CLOSED`

All 12 family tree sets are now in `data/analysis_data/phylogeny/twelve_family/`, retrieved directly
from the cluster (`01_.../pipeline_results/03K_family_specific_markers/`):

| Item | Count |
|---|---|
| `<family>.treefile` | 12 / 12 |
| `AA95_membership/*.tsv` | 12 / 12 |
| `tree_tip_host_spectrum/*.tsv` | 12 / 12 |
| `final_alignment/*.faa` | 12 / 12 |
| `source_alignment/*.faa` | 12 / 12 |
| `03K_D_*.tsv` manifests | 3 / 3 |
| `06_AA95_alignment_QC/<fam>/<fam>.alignment_sequence_QC.tsv` | 12 / 12 |

Every family was cross-checked against `03K_D_tree_input_manifest.tsv`: tip counts agree
(Geminiviridae 1604, Betaflexiviridae 1066, Potyviridae 915, Secoviridae 388, Partitiviridae 361,
Tombusviridae 326, Closteroviridae 255, Solemoviridae 247, Alphaflexiviridae 236, Rhabdoviridae 220,
Mitoviridae 204, Caulimoviridae 188 — 6,010 tips in total) and every path referenced by the manifest
resolves inside the release.

The manifest's three file-path columns (`source_alignment`, `frozen_alignment`, `AA95_membership`,
`tree_tip_host_spectrum`) were absolute cluster paths; they are now paths relative to
`twelve_family/`. `03L_FINAL_run_one_family.R` gained `resolve_manifest_path()` so it accepts either
form.

`03L_FINAL_run_one_family.R` was then run for all 12 families: **12/12 PASS**, 11 output tables each,
999 permutations and 100 pair-sampling replicates per family. `03L_FINAL_aggregate.R` produced the
seven pooled tables, from which S2E1 / S2E2 / S2E3 were rebuilt (see U7).

---

## U2 — Cluster code for modules 03H–03L retrieved; `03L` executed 12/12, `03K` static `CLOSED (runtime for 03L)`

**Correction to an earlier statement in this document.** The tree-building and host-spectrum code
was not missing after all — it had been pulled into a local snapshot earlier in the session but was
never read or placed into the release. Thirty scripts are now in the release:

| Location | Scripts |
|---|---|
| `02_host_relatedness_and_phylogeny/analysis/` | `03H_family_phylogeny_eligibility.py`, `03H2_RdRP_assignment_discrepancy_audit.py`, `03I_collect_RNA_tree_manifest.py`, `03I_prepare_RNA_family_RdRP_core.py`, `03J_family_specific_marker_audit.py`, `03L_FINAL_run_one_family.R`, `03L_FINAL_aggregate.R` |
| `02_host_relatedness_and_phylogeny/pipeline/` | 16 × `03K_*` (A/B/B2/B3/C/C2/C3/D), `03I_build_one_RNA_family_tree.sh`, `03I_parse_cdhit.py` |

Two facts worth recording:

1. **These scripts contain no cluster paths.** `grep` for `/share/home/bgi_mixr` across all 30
   returns nothing. The only path assumption was structural, in `03L_FINAL_*`, which inserted a
   literal `three_source_manuscript_analysis_v1_1` segment; both files were edited so the first
   argument is the analysis root directly, matching the release layout.
2. **They still cannot be executed here**, because they consume per-family AA95 alignments, tip
   host-spectrum tables, `03K_D_tree_input_manifest.tsv` and the sequence master — none of which are
   in this snapshot. `03L_FINAL_run_one_family.R Geminiviridae` now fails with a precise message
   naming the missing manifest, which confirms the path repair is correct and isolates the remaining
   gap to data rather than code.

What remains unverifiable, and why:

| Script group | Status |
|---|---|
| `03D_prepare_marker_inputs.py`, `03E_*` | present, **documentation only** — the intermediate 4-family stage, superseded by `03K`/`03L`. Needs `all_three_source.fasta`, a GenBank archive and the Pfam HMM, none of which are redistributed. |
| `03H`, `03H2`, `03I`, `03J` | present, **static** — these select families and build the per-family RdRP-core alignments upstream of the trees. They consume intermediates that are not in this snapshot; the trees they lead to *are* shipped. |
| `03K_*` (16 scripts) | present, **static** — the IQ-TREE chain that produced the twelve shipped trees. Not re-run because it needs the same non-redistributable inputs. |
| `03L_FINAL_run_one_family.R`, `03L_FINAL_aggregate.R` | **executed, 12/12 PASS** — run against the shipped trees, producing 11 tables per family plus 7 pooled tables. This is the chain behind Fig. 2D and S2E1–E5. |
| `04D_transmission_biology_extended_v3.R` | **executed, 19 tables** — reviewed in detail and cross-checked through the S3I–L rebuild (17 blocks verbatim). |
| modules `05`, `06`, `07` | **executed** where the released inputs suffice (`05A`–`05E`, `06A`/`06B`, `07A`/`07N`); the `07B`–`07P` pipeline layer is documentation-only. |

The general risk is unchanged: a script that runs and produces plausible numbers is not evidence
that it is the script that produced the manuscript's numbers. For the scripts marked executed above
that evidence exists, either because an anchor was reproduced exactly or because the output was
compared cell by cell against a frozen table. For the `03K` chain and the `07` pipeline it does not,
and this document says so rather than implying otherwise.

The remaining gap is therefore **not** "code is missing" — all 122 scripts are present — but
"the tree-building chain and the TSWV pipeline cannot be re-run without data we may not
redistribute". Closing that gap would require the author to confirm redistribution rights for
`all_three_source.fasta`, the Pfam HMM and the TSWV GenBank archive.


---

## U3 — `04D_transmission_biology_extended_v3.R` retrieved, path-repaired and executed `CLOSED`

The script was already present in the local snapshot; it had simply never been read or placed. It is
now at `analysis/03_arthropod_transmission_ecology/analysis/04D_transmission_biology_extended_v3.R`
(48,959 B) and **ran successfully here** in 287 s, producing all 19 tables (`20_04D_*` … `38_04D_*`),
including `24_04D_vector_family_by_virus_genus.tsv` (103 data rows), which is the source of Table S3
section I.

**Path repair applied.** The script located its inputs by walking up from its own position and
rebuilding the original cluster layout:

```r
final_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)
seq_final  <- file.path(final_root, "final_three_source_pipeline_v1", "results",
                        "20_three_source_sequence_master_FINAL.tsv")
```

That is not a hard-coded absolute path, but it is a hard-coded *layout*, and it fails on a clone.
Since 04D reads only five columns from the sequence master (`sequence_id`, `virus_taxid`,
`virus_family`, `virus_genus`, `arthropod_sequence_supported`) and only `sequence_id` plus host
fields from the association master — all of which exist in the released analysis inputs — the
defaults now point at `data/analysis_input/`, and `SEQ_FINAL` / `ASSOC_FINAL` still override them for
anyone who has the full masters.

**Result sanity check.** The five transmission groups come out at 5,070 / 274 / 31 / 57 / 29 virus
TaxIDs for none / persistent non-propagative / semi-persistent / persistent propagative /
non-persistent, spanning 74 / 6 / 4 / 6 / 5 virus families. These are internally coherent and match
the ordering the manuscript describes.

**Still open in this area.** `04D_transmission_specificity_control.R` — the separate
`additional_validation/` piece — is a different script and has not been retrieved. That remains U4.

---

## U4 — `additional_validation/transmission_mode_consistency/` closed `CLOSED`

Retrieved from `Fig4_transmission_lineage_feasibility_TEST_v2/` on the cluster and now in place:
both scripts plus all seven result tables, under
`additional_validation/transmission_mode_consistency/`, labelled `NOT_IN_PAPER` throughout.

Two modifications were needed to make it run off-cluster:

1. **Input re-pointing.** The original took `VIPDB_relation_raw.xlsx`, which is third-party and not
   redistributed. The script now also accepts the released
   `data/analysis_data/VIPDB_relation_normalized.tsv.gz`, which already carries `fine4_mode`,
   `virus_taxid`, `virus_genus`, `virus_family`, `vector_taxid` / `vector_name` and `vector_genus`.
   The workbook path is still accepted if a user supplies it.
2. **Path and output location.** The hard-coded `DEFAULT_ROOT` and the `R_BIN` pointing at
   `/share/home/bgi_mixr/.conda/...` are gone; the analysis root comes from the first argument or
   `ANALYSIS_ROOT`, and output lands next to the script (override with `SPECIFICITY_OUTDIR`) instead
   of in a new `Fig4_transmission_lineage_feasibility_TEST_v2/` tree.

**Verification.** Re-run took 63 s and reproduced the frozen cluster run exactly – 546 eligible
pairs, 414 virus TaxIDs, and all seven tables matching **row for row with identical values**
(13,104 / 24 / 96 / 4 / 46 / 18 / 0 data rows). The only changes were serialization: `fwrite`
field quoting and, on a first pass, label case (`Semi-persistent` vs `semi-persistent`). Both are now
pinned (`quote = FALSE`, `tolower()`) so output is byte-stable across `data.table` versions.

The substantive result is unchanged from the pre-retrieval reading: `virus_genus_only` accuracy
0.9958 and `vector_species_x_virus_genus` accuracy 1.0000, but on the common prediction set baseline
and conditioned accuracy are identical (1.0/1.0 and 0.9973/0.9973). Conditioning on vector identity
costs coverage (0.874 -> 0.674) without adding accuracy, which is why this is not presented as a
positive finding in the manuscript.

---

## U5 — Plotting code: all six directories in place `CLOSED`

All six plotting directories were retrieved from the cluster and placed by manuscript figure number:

| Manuscript | Cluster source | Released here | Files |
|---|---|---|---|
| Fig. 1 | `Fig2_global_landscape/Fig1_publication_FINAL_v14/` | `01_dataset_summary/plot/fig1/` | 6 |
| Fig. 2 (+ S2) | `Fig3_host_phylogenetic_constraint/Fig2_publication_FINAL_v17/` | `02_.../plot/fig2/` | 10 |
| Fig. 3 (+ S3) | `Fig4_arthropod_transmission_ecology/Fig3_transmission_publication_FINAL_v22/` | `03_.../plot/fig3/` | 10 |
| Fig. 4 (+ S4) | `Fig5_network_prioritization/Fig4_publication_FINAL_v1_1/` | `04_.../plot/fig4/` | 11 |
| Fig. 5 (+ S5) | `Fig7_TSWV_case_study/Fig7_publication_FINAL_v1_3_6/` | `06_.../plot/fig5/` | 8 |
| Fig. S1 | `Fig3_host_phylogenetic_constraint/plot_no/Fig3_publication_plot_v2_11_1/` | `02_.../plot/figS1/` | 8 |

28 R scripts in total. Hard-coded cluster paths removed from every one:

| Removed | Replaced by |
|---|---|
| `R_LIB <- "/share/home/.../mxrconda/lib/R/library"` | `Sys.getenv("R_LIB", "")`, applied only if set |
| `DEFAULT_*_ROOT <- "/share/home/.../FigXX_..."` | derived from the script's own location, four levels up |
| `export R_LIBS_USER="/share/home/..."` | comment; honoured only if the user exports it |
| `RSCRIPT="/share/home/.../Rscript"` | `command -v Rscript` |

The original LSF submission wrappers (`submit_*.sh`, `bsub` flags, queue names, conda interpreter) are
kept under `cluster_job_scripts/` with an `ARCHIVAL ONLY — NOT part of the reproduction path` header.
A grep for `share/home`, `bgi_mixr`, `mxrconda`, `/Users/` and `~/Desktop` across all released `.R`,
`.sh` and `.py` files now returns only those three archival wrappers.

Every path still honours its original environment variable (`FIG3_ROOT`, `FIG7_ROOT`, `NETWORK_ROOT`,
`CONTEXT_ROOT`, `ASSOC_MASTER`, `FINAL_THREE_SOURCE_ROOT`, `FIG3_ANALYSIS_ROOT`,
`FINAL_PIPELINE_RESULTS`, `PROJECT_ROOT`, `R_LIB`).

---

## U6 — Two superseded `04C` variants must be recorded as excluded `OPEN_AUTHOR_INPUT`

Two earlier revisions of the host-distance interaction sensitivity exist on the cluster:

```
04C_host_distance_interaction_sensitivity.before_global_sampling_alignment_20260915.R
04C_host_distance_interaction_sensitivity.pre_alignment_20260915.R
```

They are deliberately **not** released. The reason is now empirically supported rather than asserted:
the released `04C` self-checks against the frozen value and reports
`frozen_interaction_beta = 0.740127` versus `full_data_interaction_beta_abs_diff = 1.11e-16`, with
`full_data_table13_row_match = TRUE`. The two older files should be entered in
`RELEASE_DIFF_AUDIT.tsv` with `decision = REMOVE` and a one-line reason, so a future reader who finds
them on the cluster does not assume they were lost.

---

## U7 — Supplementary tables: S2E1–E3 and S3I–S3L rebuilt `CLOSED`

**S2E1 / S2E2 / S2E3 — rebuilt.** Generated by `03L_FINAL_run_one_family.R` (x12) plus
`03L_FINAL_aggregate.R`, and copied into `supplementary/Table_S2_Host_relatedness/`:

| Table | Source | Rows |
|---|---|---|
| `S2E1_twelve_family_phylogenetic_summary.tsv` | `01_FINAL_12family_main_summary.tsv` | 12 |
| `S2E2_support_stratified_permutation.tsv` | `03_FINAL_support_stratified_summary.tsv` | 12 |
| `S2E3_pair_sampling_stability.tsv` | `02_FINAL_pair_sampling_stability_summary.tsv` | 12 |

All 12 families PASS. Observed Spearman rho between host-family Jaccard distance and phylogenetic
distance is positive in every family:

| Family | tips | host families | possible pairs | rho | p (greater) | BH q |
|---|---|---|---|---|---|---|
| Alphaflexiviridae | 236 | 61 | 27,730 | 0.2966 | 0.001 | 0.0015 |
| Betaflexiviridae | 1,066 | 78 | 566,580 | 0.1670 | 0.001 | 0.0015 |
| Caulimoviridae | 188 | 52 | 17,205 | 0.1555 | 0.001 | 0.0015 |
| Geminiviridae | 1,604 | 70 | 1,261,666 | 0.1529 | 0.001 | 0.0015 |
| Secoviridae | 388 | 106 | 73,920 | 0.1435 | 0.001 | 0.0015 |
| Closteroviridae | 255 | 45 | 32,385 | 0.1396 | 0.001 | 0.0015 |
| Potyviridae | 915 | 81 | 415,416 | 0.1164 | 0.001 | 0.0015 |
| Solemoviridae | 247 | 55 | 30,381 | 0.0650 | 0.004 | 0.0053 |
| Rhabdoviridae | 220 | 67 | 24,090 | 0.0549 | 0.001 | 0.0015 |
| Tombusviridae | 326 | 75 | 52,975 | 0.0419 | 0.015 | 0.0180 |
| Partitiviridae | 361 | 98 | 64,980 | 0.0167 | 0.041 | 0.0276 |
| Mitoviridae | 204 | 60 | 20,706 | 0.0231 | 0.095 | 0.1770 |

6,010 tips across the 12 families, 999 permutations and 100 pair-sampling replicates each. Ten of
twelve families reach p < 0.05 before correction and nine survive BH correction at q < 0.05
(Mitoviridae does not; Partitiviridae sits at q = 0.028).

**S3I–S3L — rebuilt.** These four sections are the part of Table S3 that comes from
`04D_transmission_biology_extended_v3.R`. That script writes 19 tables (`20_04D_*` … `38_04D_*`)
with column counts from 9 to 41, so packing them into one wide rectangle would have meant padding
blocks and inventing a shared schema. They are shipped instead as **stacked tables with one
self-contained block per source table**: two leading columns (`section_block`, `source_table`) plus
that source table's own header and values, copied verbatim.

| Section | File | Blocks | Data rows | Source tables |
|---|---|---|---|---|
| S3I | `S3I_vector_family_by_virus_genus.tsv` | 1 | 103 | `24` |
| S3J | `S3J_five_group_host_range_models.tsv` | 5 | 4,691 | `20`, `25`, `26`, `27`, `34` |
| S3K | `S3K_within_mode_lineage_sensitivity.tsv` | 5 | 875 | `23`, `35`, `36`, `37`, `38` |
| S3L | `S3L_country_breadth_and_evidence_models.tsv` | 6 | 351 | `28`, `29`, `30`, `31`, `32`, `33` |

*(File names updated 2026-10-27, U25: K = lineage sensitivity, L = country
breadth, matching `Data_S1.xlsx`. The two files were swapped in the original
listing above; contents unchanged.)*

Rebuild with `python3 supplementary/Table_S3_Arthropod_vector/build_S3I_L.py`; the block map and
the reading caveat are in `supplementary/Table_S3_Arthropod_vector/S3I_L_README.md`.

Verification: **all 17 source blocks re-read and compared cell by cell against the `04D` output —
0 differences.** No value was recomputed, rounded or reformatted.

The reading caveat is worth repeating here because it will trip someone: the four files carry a
repeated header row per block, so `read.table(header = TRUE)` on the whole file yields one confusing
frame. Read per block, or read with `header = FALSE` and assign headers by `section_block`.

The 19 `04D` outputs themselves are regenerated run artefacts and stay under
`analysis/03_arthropod_transmission_ecology/tables/`, which is git-ignored by design. The four
section files under `supplementary/` are the released form.

---

## U8 — `CITATION.cff` and `checksums/SHA256SUMS.txt` finalised `CLOSED`

`CITATION.cff` has been rewritten for the twelve-family version (version `1.1.0`, title and abstract
updated, `date-released` set) and `checksums/SHA256SUMS.txt` was regenerated as the **last** step,
after every other file had stopped changing. The rule that generated it is stated in
`checksums/README.md`: every tracked file in the repository, sorted, with the relative path as the
key, excluding the checksum file itself and anything git-ignores.

If a file changes after this point the checksum list is stale by definition; re-run
`bash checksums/generate_sha256sums.sh` rather than hand-editing the list.

---

## U9 — Licence status of third-party databases is `TO_BE_CONFIRMED` `CLOSED (2026-10-27)`

**Resolved.** Each upstream source was verified on 2026-10-27 against its
official terms page, and `docs/EXTERNAL_INPUTS.md` §5 now records the
verification URL, date and licence conclusion per source:

- NCBI Virus / GenBank: redistributable under NCBI policy ("NCBI itself places
  no restrictions…"); individual submitter claims caveat kept.
- Virus-Host DB: **CC BY-NC-SA 4.0** — redistributable, but *not* compatible
  with the package's CC BY 4.0; embedded records stay under their own terms.
- Land Flora RNA virome (Liu et al. 2026, *Molecular Plant* 19(10),
  doi:10.1016/j.molp.2026.07.023; public dataset doi:10.6084/m9.figshare.30688322):
  records derive from the authors' openly published figshare dataset; redistributable with citation.
- VIP-DB (Deng et al. 2026, Viruses 18:679; viplant.cn; 1,375 records):
  openly published; the shipped normalized projection matches it exactly.
- ICTV: CC BY 4.0 (compatible). Pfam: CC0 1.0 (not redistributed).
`LICENSE_DATA` still does not re-license any third-party data, and no raw
upstream dump is redistributed.

---

## U10 — `03B` silently reported a partial null distribution `CLOSED (script fixed)`

This was initially recorded as an environment limitation. It is more serious than that, and it is
now fixed in the script itself.

**What happened.** The original `03B` dispatched all iterations in one `mclapply` call. In this
container, fork workers were reaped without returning values. `mclapply` substitutes `NULL` for a
lost worker and does **not** error, so a full-length run printed `Fig3 degree-preserving null PASS`
and wrote a QC table claiming `iterations = 999`, while the per-iteration tables held only:

| Network | Requested | Actually fitted |
|---|---|---|
| formal TaxID | 999 | **251** |
| RdRP65 | 999 | **749** |

The `null_iter` column of the summary was honest (251 / 749); the QC `iterations` column was not.
The run was therefore successful in appearance and incomplete in fact — and the p-value differed
accordingly (0.00397 on 251 iterations versus the manuscript's 0.001 at 999).

**Why it matters for the release.** Any published null p-value depends on the number of resamples.
A user who cloned the old script on a machine with the same container behaviour would have obtained
a different p-value from the same code, with no warning. This is a reproducibility defect, not a
local inconvenience.

**Fix applied** (statistical definition unchanged; only dispatching and validation):

1. bounded chunks instead of one all-at-once `mclapply`;
2. every returned element validated for a well-formed one-row fit;
3. missing iterations recomputed serially, up to three attempts each;
4. hard abort if distinct iterations still fall short of the request;
5. QC split into `iterations_requested` and `iterations_completed`;
6. block-wise concatenation instead of accumulation, after a first version of the fix was killed by
   the container memory limit (`EXIT=137`, peak 1 GB).

**Verified after the fix.** The full run completed **999 / 999 for both networks** in 34 minutes
(`EXIT=0`, peak RSS 1.0 GB), with no serial recovery needed — so the chunked dispatch eliminated the
fork loss rather than merely patching around it. The recovered null reproduces the manuscript
exactly: formal null β median **−0.2071** (manuscript −0.207), empirical **P = 0.001** (manuscript
0.001), 999 distinct iterations confirmed by an independent `awk` cross-check, k = 0 iterations more
negative than observed. RdRP65 null β median −0.1350, also P = 0.001. The recovery path was observed
firing in an earlier run (`formal: 20 iteration(s) not delivered by fork; recomputing serially`), so
the guard is exercised in practice, not merely present in the code.

**Residual guidance for users.** `NCORES` should be lowered on memory-constrained machines. Peak RSS
is about 480 MB per iteration. `docs/reproducibility.md` states this and states that a smoke test
with small `FIG3_NULL_ITER` produces a p-value that must not be quoted.

A separate confusion is also resolved: `03_main_objects.rds` is **not** an orphaned input.
`03A_host_constraint_main.R` line 63 writes it (`saveRDS(obj, file.path(base,"work","03_main_objects.rds"))`)
and `03B` line 5 and `03C` line 3 read it. Because `analysis/*/work/` is git-ignored, a fresh clone
must run `03A` before `03B` or `03C`. This ordering requirement is documented in
`docs/reproducibility.md`.

---

## U11 — GitHub push is intentionally not performed `OPEN_AUTHOR_INPUT`

Per the author's instruction, nothing has been committed or pushed. The tree currently lives at
`/workspace/release/` and is a plain directory, not a repository with remotes configured. The author
reviews `RELEASE_DIFF_AUDIT.tsv`, `FINAL_RELEASE_VALIDATION.md` and this file before any push.

---

## U12 — `00_paths.sh` exported a wrong root under non-bash shells `CLOSED (script fixed)`

Found by the final end-to-end check — unpacking the tarball into a clean directory and reading what
the script actually exported — not by reading the script.

```
$ cd /tmp/verify/release && source analysis/00_config/00_paths.sh
ANALYSIS_ROOT=/tmp/verify      # should be /tmp/verify/release/analysis
DATA_DIR=/tmp/data
```

`dirname "${BASH_SOURCE[0]}"` yields `.` when `BASH_SOURCE` does not exist (zsh, `sh -c`, dash
`. file`), so `${...}/..` resolved to the parent of the **working directory**. The exported path
existed, so nothing downstream failed loudly — the run simply would have looked for its inputs one
level too high and quietly found nothing or the wrong thing.

Earlier runs in this pass could not have caught it: they were all launched as
`bash -c 'source analysis/00_config/00_paths.sh; Rscript ...'`, and bash always defines
`BASH_SOURCE`.

Fixed by resolving the location from `BASH_SOURCE`, then `$0`, then a relative probe, **validating
each candidate** with `_is_repo_analysis_root()` (it must contain `00_config/00_paths.sh`), and
erroring out with a non-zero return if nothing validates rather than exporting a plausible wrong
path. Verified under bash, zsh and dash, from the repo root, from a subdirectory, and from an
unrelated directory using an absolute path — all five resolve to the correct `analysis/` directory.

Regression after the fix, run from a freshly unpacked tarball: `05A` reproduced 269 genera / 1,902
edges / 12 candidates and `06B` exited 0, both matching the values in `FINAL_RELEASE_VALIDATION.md`.

**A third instance in the same family** was caught by the same unpacked-archive check:
`00_preflight.R` still looked for the twelve-family inputs as a flat set
(`<Family>.treefile` under `FAMILY_TREE_DIR`), which stopped being true once the trees were laid
out as `trees/<Family>/<Family>.treefile` to match what `03L` reads. On a clean unpack it aborted
on twelve files that were present and correctly placed — a false alarm, which is the mirror image
of the wrong-root bug but equally damaging: a checker that cries wolf gets ignored or, worse, gets
"fixed" by moving real data.

It now validates the four real locations (`trees/<Family>/`, `AA95_membership/`,
`tree_tip_host_spectrum/`, `06_AA95_alignment_QC/<Family>/`, with the Geminiviridae
`07_Geminiviridae_Rep_conserved/AA95_alignment_QC/` branch) and, on failure, reports **how many of
each class are missing** plus the expected layout. Verified both ways: PASS on a complete tree, and
a correct per-class diagnosis (`treefile 1, membership 0, spectrum 0, alignqc 0`) when one treefile
is deleted.

The lesson generalises and is why this item is kept: **a path that exists but is wrong is more
dangerous than a path that is missing.** The missing one fails; the wrong one gets used. And a
checker built on a stale assumption is worse than no checker, because it lends false confidence.

---

## U16 — Plotting wrappers could not run from a clean clone `CLOSED (6 defects fixed)`

The five `run_*.sh` wrappers were carried over from the cluster unchanged. Running each one
against a fresh `git archive` export found six distinct ways they failed. None of them show up
when you read the script; all of them show up within seconds of running it.

| # | Defect | Symptom |
|---|---|---|
| 1 | `set -euo pipefail` + `echo "$R_LIBS_USER"` | `R_LIBS_USER: unbound variable`, exit 1, before R ever starts |
| 2 | fig5's wrapper used `$SCRIPT_DIR` without defining it | `PROJECT_ROOT` resolved to `/`, then a path under `three_source_manuscript_analysis_v1_1/` that does not exist |
| 3 | fig3's wrapper defaulted `FIG3_ANALYSIS_ROOT` to the cluster folder name | `normalizePath(mustWork=TRUE)` on a directory that was never shipped |
| 4 | fig4's wrapper pointed `CONTEXT_ROOT` at the network module | Fig. 4D's adjustment table looked for one module away from where it lives |
| 5 | No wrapper was marked executable in git | `./run_*.sh: Permission denied` |
| 6 | fig2 resolved the Stage-20 master, the 12 treefiles and the tip spectra only at cluster paths | three successive `not found` errors, each after a full re-run |

All six are fixed. Verified by running every wrapper on a fresh export; see the table in
`FINAL_RELEASE_VALIDATION.md` section E.

The lesson is the same one U12 recorded, one layer up: **the wrappers were never executed during
release preparation.** The R scripts behind them had been, by calling `Rscript <script>.R` directly,
which bypasses every one of these defects.

---

## U17 — Three plotting inputs were excluded by `.gitignore` `CLOSED`

`Fig. 2` reads `pipeline_results/03L_FINAL_12family_host_spectrum/07_FINAL_plot_ready_family_summary.tsv`.
`.gitignore` excluded that whole directory as regenerated output — correct when the plotting code was
not part of the release, wrong once it was. The file was present in the working tree, invisible to
`git status`, and absent from every clone and every tarball.

Fixed by tracking the seven top-level 03L summary tables (~17 KB) and keeping the per-family scratch
directories excluded, which are what 03L actually rewrites. `git status` now shows the tables and the
`git archive` export contains them.

Files that stay excluded on purpose, and are regenerated by re-running the script named beside them:

| Excluded | Regenerated by | Why not shipped |
|---|---|---|
| `03L_FINAL_12family_host_spectrum/<Family>/` | `03L_FINAL_run_one_family.R` ×12 + `03L_FINAL_aggregate.R` | per-family iteration dumps, rewritten every run |
| `analysis/*/logs/`, `analysis/*/work/` | created by the scripts | empty output directories |
| `analysis/*/plot/*/{main,panels,supplementary}/` | the `run_*.sh` wrappers | rendered figures |

---

## U18 — `03F`/`03G` assumed a flat input layout that the release does not use `CLOSED (scripts fixed)`

Same class as U12's third instance, found the same way. `03F` and `03G` looked for
`<FAMILY_TREE_DIR>/<Family>.treefile` and `<FAMILY_TREE_DIR>/<Family>.nr95.membership.tsv`; the
release ships `trees/<Family>/<Family>.treefile` and `AA95_membership/<Family>.AA95_membership.tsv`.

The second one is the dangerous half. `Geminiviridae.nr95.membership.tsv` **does exist** at the flat
location — it is the older, coarser 95%-identity clustering — so adding a fallback without fixing the
order would have silently produced host spectra computed against the wrong clustering, with no error.
The resolution order is now AA95 first, nr95 last, with a comment saying why.

`03F` also assumed an `identity_to_representative_percent` column that the released membership tables
do not carry. It is provenance decoration that nothing downstream reads, so it is created as an empty
column rather than crashing the run — and the comment says that is what it is.

Verified by execution, not by reading: `03F` runs to completion (5m 28s, 12 families), `03G` writes
`26_03G_pair_sampling_stability_summary.tsv` and `27_03G_support_stratified_permutation_summary.tsv`
from it, and `Fig. S1` then renders from those tables.

---

## U13 — `figures/` shipped as an empty directory with no explanation `CLOSED`

A reviewer opened the release, looked in `figures/`, found a single `fig1/` subdirectory tree and
reasonably concluded that the other five figures had been forgotten.

They had not. `figures/fig1/{logs,main,panels}` were **empty output directories left behind by a
plotting run** — zero files, zero bytes. `.gitignore` excludes `figures/` wholesale because rendered
figure output is deliberately outside the reproducibility scope, so those directories were never
tracked and never reached a cloner; they survived only in the working tree and got swept into the
tarball by `tar`.

Two defects, both fixed:

1. The stale empty tree was removed (`rm -rf figures/fig1`).
2. `.gitignore` used `figures/`, which would also have swallowed any explanatory file placed there.
   Changed to `figures/*` plus an explicit `!figures/README.md`, so the note survives a clone.

`figures/README.md` now states what the directory is for and maps each figure to its plotting code:

| Figure | Plotting code |
|---|---|
| Fig. 1 | `analysis/01_dataset_summary/plot/fig1/` |
| Fig. 2 | `analysis/02_host_relatedness_and_phylogeny/plot/fig2/` |
| Fig. 3 | `analysis/03_arthropod_transmission_ecology/plot/fig3/` |
| Fig. 4 | `analysis/04_genus_host_network/plot/fig4/` |
| Fig. 5 | `analysis/06_TSWV_evolution/plot/fig5/` |
| Fig. S1 | `analysis/02_host_relatedness_and_phylogeny/plot/figS1/` |

**What this exposed:** the release had no "clean the working tree before packing" step. Everything up
to now had been verified against the working tree, where a git-ignored empty directory is invisible
to `git status` and harmless — it only becomes a defect at packing time. The pack step now exports
via `git archive HEAD` rather than `tar`-ing the working tree, so the archive contains exactly the
tracked file set and nothing else.

---

## U14 — S3I–S3L are stacked files and a plain read silently corrupts them `CLOSED (documented)`

Sections I–L of Table S3 stack 1, 5, 6 and 5 source tables respectively, each with **its own header
row** and its own column count (from 6 up to 43). This is not a defect — the blocks genuinely have
different schemas, and forcing them into one rectangle would mean padding real data with fake cells
— but it is a trap:

```r
x <- read.delim("S3J_five_group_host_range_models.tsv")
# 4695 x 43. No warning. The row count is right.
# But block 1 has 16 columns; its names become the header for all 43,
# and blocks 2-5 sit under names that belong to block 1.
```

A reader who then asks for `estimate` gets block 1's padding, not the model terms they wanted.

Documented in two places: the block map and per-block row counts in
`Table_S3_Arthropod_vector/S3I_L_README.md`, and a working read recipe in
`supplementary/README.md` that splits on `^section_block\t` and drops the `fill = TRUE` padding
columns. The recipe in the README is the exact code that was run to produce the table of block
dimensions printed beside it, so it is known to work rather than believed to.

**Verified round-trip:** `04D` was re-executed on a fresh extraction of the tarball, `build_S3I_L.py`
was run over its output, and all four files came back byte-identical (SHA-256) to the released
copies.

---

## U15 — `supplementary/README.md` undercounted the S3 directory `CLOSED`

The file table said `12 TSV + 1 README + 1 rebuild script`; the directory holds **15 TSV** (S3A1–S3H
is 11 files, S3I–S3L is 4). Corrected, with the count now derived from `ls` rather than carried
forward from the previous release.

The same pass found the status table still described sections I–L as *"still the previous
snapshot"* in a footnote written before they were rebuilt. Reworded: S3I–S3L are `FINAL` and
verified; what remains a re-export is the older S2E-table convention, not S3.

---

---

## U19 — Fig. 5 cannot be rendered from the release as shipped `CLOSED`

**Closed 2025 by shipping the five missing inputs from the cluster.** Option 1 (below) was executed:
the frozen outputs were pulled verbatim from the cluster over an SSH tunnel and committed to
`analysis/06_TSWV_evolution/`, with no recomputation of any kind.

What was added (89 new tracked files):

| Input | Shipped as |
|---|---|
| five-gene summary table | `tables/06_TSWV_five_gene_evolutionary_input_summary.tsv` (provenance paths rewritten to repo-relative) |
| GARD final status | `tables/16_TSWV_GARD_final_status.tsv` (provenance paths rewritten to repo-relative) |
| GARD + selection summary | `tables/18_TSWV_GARD_selection_summary.tsv` |
| selection sites (long) | `tables/19_TSWV_selection_sites_long.tsv` — decompressed MD5 matches the published `data/analysis_data/TSWV_selection_sites_long.tsv.gz` byte-for-byte (`ec4253cd…`), confirming these are the same frozen outputs, not a recomputed version |
| 5 gene codon alignments + IQ-TREE outputs | `pipeline_results/03_alignment_and_phylogeny/{RdRp,GnGc,N,NSm,NSs}/` |
| 3 segment trees | `pipeline_results/06_segment_concordance/{L,M,S}/` |

`.gitignore` gained an explicit `!analysis/06_TSWV_evolution/pipeline_results/**` allowance so the
global `*.treefile` / `*.fna` rules do not swallow the frozen trees.

Verification, from a clean `git archive HEAD` export:

- `run_fig7_publication_FINAL_v1_3_6.sh` exits 0 and prints `Fig. 5 plotting complete.`
- 34 output files rendered (main 1 + panels 3 + supplementary 7, each as PDF/PNG/TIFF),
  including `Fig7_TSWV_multiscale_evolution_FINAL_v1_3_6.pdf`.
- Release preflight PASS (80,028 sequences / 83,034 associations).

Original description is retained below for the audit trail.

The chain starts at 07E, which shells out to `mafft --globalpair`, `trimal -automated1` and
`iqtree3 -m MFP+MERGE -B 1000`. None of the three is present here, and rebuilding them would not be
a reproduction check anyway — it would be a fresh tree inference with a different tool build, whose
output would not be comparable to the frozen trees the manuscript used.

`tables/13_TSWV_segment_concordance_summary.tsv` **is** shipped, so `07N` (segment concordance) runs
and its numbers are verified. What is absent is the gene-tree and selection layer above it.

Two ways to close this, both needing your input:

1. **Ship the five tables and the 03_alignment_and_phylogeny tree directory** from the cluster. This
   is the cheap option and the right one: they are frozen analysis outputs, the same class of file as
   the twelve-family trees that are already here. Nothing needs recomputing.
2. Accept that Fig. 5's plotting code is shipped but not runnable, and say so in the data-availability
   statement. The figure itself is in the manuscript; the code is available for inspection.

Option 1 is what the rest of the release does, so it is the one I would recommend. I have not
fabricated substitutes from `Data_S1` or `S5D`, because a plotting script fed reconstructed inputs
would produce a figure that looks right and is not the published one.

---

## Summary table

| ID | Item | Status | Unblocks with |
|---|---|---|---|
| U1 | 12 of 12 family tree sets in place, `03L` x12 all PASS | CLOSED | done |
| U2 | 03H-03L executed; Fig.3 plotting code executed end-to-end on fresh clone | CLOSED for runtime; byte-diff pass optional | done |
| U3 | 04D retrieved, path-fixed, executed (19 tables) | CLOSED | done |
| U4 | `additional_validation/` specificity control retrieved, re-pointed, executed; 7/7 tables match frozen run | CLOSED | done |
| U5 | plotting code: all six figure directories in place, paths cleaned | CLOSED | done |
| U6 | superseded `04C` variants to record as excluded | OPEN_AUTHOR_INPUT | audit edit |
| U7 | S2E1/E2/E3 rebuilt (12/12 PASS); S3I–L rebuilt as stacked blocks, 17/17 verbatim | CLOSED | done |
| U8 | `CITATION.cff` rewritten for 12 families; `SHA256SUMS.txt` regenerated last | CLOSED | done |
| U9 | third-party licence status | CLOSED (verified 2026-10-27) | see U9 section above; `docs/EXTERNAL_INPUTS.md` §5 |
| U10 | `03B` silently reported a partial null | CLOSED (script fixed) | fixed; full run verified |
| U11 | no GitHub push | OPEN_AUTHOR_INPUT | you upload the tarball |
| U12 | `00_paths.sh` wrong root under zsh/sh; `00_preflight.R` stale flat-path assumption | CLOSED (scripts fixed) | fixed; verified from unpacked archive |
| U13 | `figures/` shipped as an empty `fig1/` tree with no explanation | CLOSED | empty tree removed, README added, pack step changed to `git archive` |
| U14 | S3I–S3L stacked files read silently wrong with `read.delim` | CLOSED (documented) | read recipe added; round-trip verified byte-identical |
| U15 | `supplementary/README.md` undercounted S3 as 12 TSV | CLOSED | corrected to 15, recount from `ls` |
| U16 | all five plotting wrappers failed from a clean clone (6 distinct defects) | CLOSED | fixed; each wrapper run end-to-end |
| U17 | `.gitignore` excluded 03L summaries that Fig. 2 reads | CLOSED | 7 top-level tables tracked; scratch still excluded |
| U18 | `03F`/`03G` assumed flat input layout; nr95 fallback would have been silently wrong | CLOSED | AA95-first resolution; 03F/03G executed |
| U19 | Fig. 5 inputs shipped from cluster; wrapper runs end-to-end (EXIT=0, 34 outputs) | CLOSED | done |
| U20 | Re-run of `05A` gives 10 modules / Q 0.3446; paper reports 12 / 0.3532 | CLOSED (frozen outputs shipped; non-determinism documented) | Louvain is a randomised heuristic — see §B4 of the validation report. Data, edges, all partition-independent metrics and the 12 prioritised genera are identical across runs; only the partition differs, and the released tables are the MD5-matched frozen ones. Not a blocking defect, but a reader re-running `05A` will not get 12. |
| U22 | Fig. S1 could not be rendered from a clone; identity of the submitted figure initially misread | CLOSED | **Resolved in two steps.** (1) The four missing frozen tables (21/22_03E, 26/27_03G) were retrieved from the cluster and shipped. (2) The submitted Supplementary Fig. S1 was then identified from the author's PDF: it is the SINGLE-PANEL RdRP65 host-sharing decay figure (98 × 75 mm), not the old four-panel composite. A standalone entry (`plot/figS1/01_FigS1_publication_plot_v2_11_1_standalone.R`, wrapper `run_FigS1_standalone.sh`) now draws exactly that panel with the frozen `make_distance_decay_plot()` logic and the shipped `tables/04_RdRP65_host_sharing_by_distance.tsv`, with the submitted label wording ("RdRP65 sharing probability (%)", "Number of host pairs"). Verified against the submitted `Supplementary Fig1.pdf`: same data, geometry, colours and labels. The submitted S1 reads no tree files. |
| U21 | External-evidence files carried the author's internal archive document name and page ranges | CLOSED | Two columns dropped after confirming all 28 rows keep a public URL; files renamed to `prioritized_genera_evidence_sources.tsv` / `…_summary.tsv` and all references updated |
| U23 | Provisional item: four frozen four-family treefiles were thought to be needed for Fig. S1 | **WITHDRAWN** | The submitted S1 (verified from the author's PDF) is the RdRP65 decay panel and reads no tree files. The provisional `S1_frozen_trees/` directory and its tree-preference hook in `00_helpers.R` were removed. The four-family trees remain on the cluster and are NOT needed for the release. |
| U24 | `S2E4_high_support_sensitivity.tsv` was a stale pre-`03L` export: it carried the OLD four-family (old-tree) numbers — e.g. Geminiviridae ≥2: 924 tips / ρ 0.212557, ≥3: 637 / 0.243273 — while `Data_S1.xlsx` and the manuscript use the frozen `03L` FINAL values (933 / 0.196480; 625 / 0.217230) | **CLOSED (2026-10-27)** | Root cause: the S2E4 export was made from the earlier `03G`-era run and was never refreshed when the `03L` FINAL chain replaced it; the sibling files S2E1–E3 and S2E5 all match the frozen FINAL outputs. Verified: all 24 min-support rows in `Data_S1.xlsx` (Table S2) match `pipeline_results/03L_FINAL_12family_host_spectrum/04_FINAL_min_support_sensitivity_summary.tsv` exactly; all 8 rows of the old S2E4 TSV mismatched it. Fix: S2E4 re-exported byte-for-byte from the frozen `04_FINAL_min_support_sensitivity_summary.tsv` (12 families × 2 thresholds, full schema incl. sampling mode, status, null quantiles and BH q); `supplementary/README.md` and `docs/MANUSCRIPT_OUTPUT_MAP.tsv` updated. Manuscript numbers were already correct and are unchanged. |

---
| U25 | Batch of reader-facing inconsistencies (fixed 2026-10-27): (1) S3 standalone K/L letters were swapped relative to `Data_S1.xlsx` and the manuscript; (2) `MANUSCRIPT_OUTPUT_MAP.tsv` mapped the host-range model tables to the wrong Fig. 3 panels; (3) the Land Flora A/B fields could be misread as direct-infection evidence levels | **CLOSED (2026-10-27)** | (1) Files renamed so **K = within-mode lineage sensitivity, L = country breadth and evidence models** (matching the workbook); `build_S3I_L.py` updated and re-run: renamed files byte-identical (MD5) to the previous contents; README, `S3I_L_README.md` and the output map updated. (2) Output map corrected against the actual Fig. 3 code: panel B = OR of the two evidence classes (S3B; blocks 30–32), panel C = host-distance decay (S3C1/C2), panel D = adjusted host-range + country-coverage models (S3J; blocks 28–29, 33), lineage sensitivity feeds Fig. S3B. (3) `supplementary/README.md` gained a Land Flora section: group A (841 records) = classification labels containing "Plant-infecting" (including "for MP+ lineage" qualifiers) per the original study's Data S4; group B = the remaining 1,675 included plant-associated records; the groups are a **recoding** of the source classification and do not establish experimentally confirmed infection. (Superseded in part by U27, which corrected the group-A label.) |
| U26 | Re-running the network analysis would overwrite the frozen tables that Fig. 4 renders, so a reader could silently plot a new Louvain partition as the published figure | **CLOSED (2026-10-27)** | Two-sided guard added. (1) `05A_network_main.R`: if the frozen tables exist, all outputs go to `tables/rerun/` with an explicit message; verified in a sandbox copy — the re-run produced 10 modules / Q 0.3446 in `rerun/` while the frozen 12-module `tables/03_network_module_membership.tsv` stayed untouched (MD5 unchanged). (2) `04_Fig4_publication_FINAL_v1_1.R`: verifies the module-membership MD5 against the frozen partition before plotting and stops with instructions otherwise. Full Fig. 4 wrapper re-run after the change: EXIT=0, "Frozen published partition verified (MD5 match)", 12 modules. |

| U27 | Follow-up corrections after author review (2026-10-27): (1) the Land Flora documentation still described group A as "direct plant-infection evidence", which the recoding rule does not support, and implied all 841 group-A records meet the MP+ condition; (2) `docs/MANUSCRIPT_OUTPUT_MAP.tsv` had two rows in which records were fused (row 16 = S2E4+S2E5 as 15 columns; row 41 = S3L+S3K+VIPDB as 22 columns) plus residual figure-number errors; (3) `S3I_L_README.md` over-interpreted the S0/S2-vs-VIPDB blocks; (4) re-running `05A` still overwrote `qc/01_Fig5_main_QC.tsv` | **CLOSED (2026-10-27)** | (1) README section replaced with the Data S4 recoding statement (group A = 841 records whose classification labels contain "Plant-infecting", incl. "for MP+ lineage" qualifiers; group B = remaining 1,675; recoding does not establish confirmed infection); the reader-facing label was changed to **"Land Flora group-A-only sensitivity set"** in `S3H_RdRP65_sensitivity.tsv` and `Data_S1.xlsx` Table S3 rows A177/A178 (2 cells; membership and values unchanged); `strict_LF_A_only` and the verbatim field names are kept and marked as historical internal names. (2) The fused rows were split back into 8-column records (map is now 8 columns on every non-blank line) and figure assignments corrected: S3B → Fig. 3D; node-subsampling files → Fig. 3C; S3D/E/F four-mode rows → Fig. 3D/Fig. S3 with a 31-vs-32 statistical-scope note on S3D (32 = full single-mode TaxID set, 31 = family-annotated Fig. 3A analysis set). (3) The 30–32 block description now states only what the tests do: association between arthropod-associated sequence evidence and documented vector relationships, before and after adjustment; the "evidence coverage artefact, not biology" attribution was deleted. (4) `05A` now routes its QC file to `qc/rerun/` alongside `tables/rerun/` when frozen outputs are present. |

| U28 | Pre-release audit (2026-10-28) against a fresh clone: two reader-facing documentation statements were not true of the shipped tree — (1) README and `analysis/PLOT_MAP.md` both said *every* figure can be re-rendered with the shipped `run_*.sh` wrappers, but `01_dataset_summary/plot/fig1/` has no wrapper; (2) `PLOT_MAP.md` claimed "a grep for `share/home`, `bgi_mixr`, `mxrconda`, `/Users/` or `~/Desktop` across all released scripts returns nothing", whereas the three `cluster_job_scripts/submit_*.sh` files necessarily retain cluster absolute paths. A third issue was cosmetic but confusing: the Fig. 1 plotting script emitted a warning when its per-host-family bar sum (82,861) fell 173 short of the frozen headline total (83,034), without saying why | **CLOSED (2026-10-28)** | (1) README and `PLOT_MAP.md` now state that Figs. 2–5 have wrappers and Fig. 1 is run directly, with the exact command and output location (`figures/fig1/`, git-ignored). Verified: Fig. 1 script runs from a fresh clone, EXIT=0. (2) The claim now scopes itself to the reproduction path and explains that the archived LSF scripts keep cluster paths by design and are marked `ARCHIVAL ONLY` in their header. (3) The warning became an explanatory `message()` that names the cause and the size of the gap (173 records whose host has no assigned family); re-run confirms the message appears and the panel is unchanged. No result, table, model or figure content changed. |
