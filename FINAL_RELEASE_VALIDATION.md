# Final release validation

Scope: the release assembled at `/workspace/release/` for the manuscript
*Phylogenetic structure, lineage differences and transmission ecology of global plant virus host
associations*. This report records what was **actually executed and observed** in this container,
what was verified only by static inspection, and what could not be verified at all. Nothing in this
document is inferred from filenames, timestamps or checksums alone; the method for choosing the
final code version is described in `RELEASE_DIFF_AUDIT.tsv` and rests on script input/output paths,
execution logs, frozen result files and the final Data S1.

Three verdicts are used and they are not interchangeable:

| Verdict | Meaning |
|---|---|
| **PASS (runtime)** | The script was executed here, exited 0, and its outputs were compared against the frozen manuscript anchors. |
| **PASS (static)** | The script was read and its logic traced against frozen outputs and the manuscript text, but it was not executed (missing input or external dependency). |
| **NOT VERIFIED** | Neither executed nor traceable to a frozen output in this snapshot. |

The governing constraint throughout: an old public code version must never be allowed to overwrite a
final result, and a missing file is recorded as missing rather than replaced by something plausible.

---

## A. Environment and method

| Item | Value |
|---|---|
| OS | Ubuntu 22.04, linux/amd64 |
| R | 4.3.3 |
| R packages | data.table 1.14.10, Matrix 1.6.5, MASS 7.3.60.0.1, ape 5.7.1, igraph 1.6.0 |
| Python | 3.11.0rc1, pandas / numpy / openpyxl available |
| Cores | 32 physical, 123 GB RAM |
| Path resolution | `source analysis/00_config/00_paths.sh`, then `Rscript ...`; no path is hard-coded |
| Analysis root | `/workspace/release/analysis` |
| Data root | `/workspace/release/data` |

Every command in section B has the form:

```bash
cd /workspace/release
bash -c 'source analysis/00_config/00_paths.sh; Rscript <script>'
```

`00_config.R` refuses to start if `ANALYSIS_ROOT` is unset, so a run that silently used the wrong
root would fail loudly rather than produce numbers against the wrong inputs.

---

## B. Module-by-module results

### B1. Module 01 — dataset summary (`Results 3.1`, Fig. 1)

| Script | Verdict | Observed |
|---|---|---|
| `01_Fig1_data_framework_QC.R` | **PASS (runtime)**, 27 s | 80,028 sequence records; 83,034 association records; 11,744 unique `vc80_final_id`; 11,708 VC80 analysis units; 27,549 Stage-8 RdRP core sequences; 5 ineligible; 0 conflicting; 33 / 38 / 50,807 / 72,494 for the evidence-tier splits |
| `02_Fig2_global_landscape.R` | **PASS (runtime)** | virus families 75; host families 296; countries 184; continents 7; RNA-classified 50,807; rarefaction 100 permutations |

Every one of these values matches the frozen anchor and the manuscript's §3.1 and §3.2 text. The
rarefaction loop is seeded (`SEED = 20260813`) and produced a stable summary across the 100
permutations.

> Deliberate note on 75 vs 76: the count returned here is 75 for virus families carrying an
> *analysis-eligible* label, while the raw label vocabulary contains 76 distinct strings. This
> distinction is already handled in the manuscript and is not a discrepancy.

### B2. Module 02 — host relatedness and phylogeny (`Results 3.3`, Fig. 3, Table S2)

| Script | Verdict | Observed |
|---|---|---|
| `03A_host_constraint_main.R` | **PASS (runtime)**, 8.8 s | assoc rows 83,034; rows with complete host + formal TaxID 71,976; host TaxIDs 2,282; formal virus TaxIDs 5,547; host pairs 2,602,621; RdRP65 lineages 1,353; RdRP65 host pairs 1,418,770; per-family slopes 24 |
| `03B_degree_preserving_null.R` | **PASS (runtime)**, 34 min | formal null **999/999** iterations; RdRP65 null **999/999**; `median_swap_completion = 1` both; formal null β median **−0.2071** [−0.2596, −0.1459], empirical p = **0.001**; RdRP65 null β median −0.1350, p = 0.001 |
| `03C_sensitivity.R` | **PASS (runtime)**, 88 s | host70 median β −1.1962 [−1.2479, −1.1451]; virus70 −1.1822 [−1.3699, −1.1236]; exclude top-1% host degree −1.2703; exclude top-1% virus degree −1.5054; 50/50 valid iterations per subsample scheme |
| `03L_FINAL_run_one_family.R` (x12) | **PASS (runtime)**, 10 min total, 12/12 families | 11 tables per family, 999 permutations and 100 pair-sampling replicates each. Tip counts match `03K_D_tree_input_manifest.tsv` exactly; 6,010 tips across 12 families. Observed rho positive in all 12 families (0.0167–0.2966). **Three distinct significance counts, not interchangeable:** main test 10/12 at p<0.05 and **11/12 at BH q<0.05** (`S2E1`); support-stratified sensitivity **9/12 at q<0.05** (`S2E2`); minimum-support sensitivity **12/12 at q<0.05 when each tip needs ≥2 supporting sequences** (10/12 at ≥3; `04_FINAL_min_support_sensitivity_summary.tsv`). |
| `03L_FINAL_aggregate.R` | **PASS (runtime)**, < 5 s | 7 pooled tables; S2E1/S2E2/S2E3 rebuilt from them. Note `S2E2`/`S2E3` come from this 03L chain (it reads the per-family `04_pair_sampling_stability_summary` / `05_support_stratified_summary`), **not** from `03G`, which is an independent earlier implementation. |
| `03F_phylogeny_host_spectrum_structure.R` | **NOT VERIFIED** | Superseded for the published numbers. The released Fig. 2D / Table S2E1–E5 results come from the `03L` chain above, which was executed here (12/12 PASS). `03F`/`03G` remain in the release as the earlier implementation of the same idea and are documented as such. |
| `03G_phylogeny_host_spectrum_sensitivity.R` | **NOT VERIFIED** | Same as `03F`. |

The two numbers that carry the manuscript's main quantitative claim reproduce **exactly**:

| Quantity | Manuscript | Recomputed here | Agreement |
|---|---|---|---|
| Global decay slope β | −1.2014780 | −1.20147801 | exact to 8 d.p. |
| Odds ratio per distance step | 0.3007494 | 0.30074937 | exact to 8 d.p. |
| RdRP65-only β | −0.5981323 | −0.59813227 | exact to 7 d.p. |
| RdRP65-only OR | 0.5498376 | 0.54983763 | exact to 7 d.p. |

The categorical-distance sensitivity (no linearity assumption, distance 0 as reference) is monotone
across all six contrasts: β = −1.118, −2.626, −4.720, −4.902, −9.028, −11.009 for distances 1–6, all
with p < 1e-14. This matters because it shows the decay is not an artefact of forcing a single
linear slope onto a non-linear relationship.

The degree-preserving null reproduces the manuscript's reported figure exactly:

| Quantity | Manuscript (`MANUSCRIPT_OUTPUT_MAP.tsv`) | Recomputed here |
|---|---|---|
| null β median | −0.207 | **−0.2071** |
| empirical P | 0.001 | **0.001** |
| iterations | 999 | **999 / 999 completed** |

Getting that agreement required first repairing the script, because an earlier run of the
unmodified version reported success while silently using a fraction of the requested resamples.
The full account is below; it is the most important methodological finding of this pass.

A first attempt at `03B` did fail in this container, and the reason matters beyond this machine.
The original script used a single `mclapply` call over all iterations. When fork workers are reaped
without returning a value, `mclapply` substitutes `NULL` and — crucially — does **not** raise an
error. One full-length run therefore finished with `Fig3 degree-preserving null PASS` and a QC table
claiming 999 iterations, while the actual per-iteration tables contained only **251** fitted
iterations for the formal network and **749** for RdRP65:

```
Warning message:
In mclapply(seq_len(iter), worker, mc.cores = min(NCORES, iter),  :
  scheduled cores 3, 4, 5, 7, 8, 9, 10, 11, 12 did not deliver results,
  all values of the jobs will be affected
```

This is exactly the failure mode the audit is meant to catch: a run that looks successful, writes a
summary with a plausible p-value, and is silently based on a quarter of the requested resamples. The
`null_iter` column recorded the truth (251) while the `iterations` column recorded the request (999),
and nothing reconciled them.

`03B` was therefore changed — the statistical definition of the null is untouched, only the
dispatching and validation were fixed:

1. iterations are dispatched in bounded chunks instead of one all-at-once `mclapply`;
2. every returned element is checked for a well-formed one-row fit, and any missing iteration is
   recomputed serially with up to three attempts;
3. the run **aborts** if the number of distinct iterations still falls short of the request;
4. the QC table now reports `iterations_requested` and `iterations_completed` as separate columns,
   so a short run can never again be mistaken for a complete one;
5. per-iteration structures are concatenated block by block rather than accumulated, because a
   first version of the fix was killed by the container memory limit (`EXIT=137`, peak 1 GB).

The reproduction evidence in section B for `03B` refers to the run of this corrected script, not to
the earlier short run. The distinction between a passing test and a complete test is the single most
important lesson from this pass, and it is why the QC table carries both counts.

`03_main_objects.rds` was initially suspected of being an orphaned input. It is not:
`03A` line 63 writes it and `03B` line 5 / `03C` line 3 read it. Because `analysis/*/work/` is
git-ignored, a fresh clone must run `03A` first. That ordering is now documented in
`docs/reproducibility.md`.

### B3. Module 03 — arthropod transmission ecology (`Results 3.4`, Fig. 4, Table S3)

| Script | Verdict | Observed |
|---|---|---|
| `04A_VIPDB_normalize.R` | **PASS (runtime)**, fallback exercised | raw VIPDB workbook absent by design; released normalized outputs restored, 1,375 relations |
| `04B_arthropod_transmission_analysis.R` | **PASS (runtime)**, 45 s | formal virus TaxIDs 5,979; sequence-supported 340; VIP-DB documented vector **410**; four-mode single TaxID 392; single group 401; four-mode multi 0; RdRP65 rows main+strict 2,923 |
| `04C_host_distance_interaction_sensitivity.R` | **PASS (runtime)**, 54 s | observed Δβ 0.740127; 100/100 iterations for both subsampling schemes; 0 failed |
| `04D_transmission_biology_extended_v3.R` | **PASS (runtime)**, 287 s | 19 tables produced (`20_04D_*` … `38_04D_*`), including `24_04D_vector_family_by_virus_genus.tsv` (103 rows, the source of Table S3 section I). Five transmission groups: 5,070 / 274 / 31 / 57 / 29 virus TaxIDs across 74 / 6 / 4 / 6 / 5 families. |
| `additional_validation/.../04D_transmission_specificity_control.R` | **PASS (runtime)**, 63 s, `NOT_IN_PAPER` | 546 eligible pairs / 414 virus TaxIDs. All 7 output tables match the frozen cluster run row for row (13,104 / 24 / 96 / 4 / 46 / 18 / 0 rows, identical values). Re-pointed from the non-redistributed `VIPDB_relation_raw.xlsx` to the released `data/analysis_data/VIPDB_relation_normalized.tsv.gz`. Strict-rule accuracy with lineage alone: 1.000 / 0.997; adding vector identity does not improve it. |

`04A` deserves an explicit note because it does *not* re-run the normalisation when the upstream
workbook is missing: it restores byte-identical released outputs and says so in the log. This is a
deliberate design decision so that downstream modules remain runnable on a clone, and it does not
alter any normalisation logic. A user who wants to rebuild from the raw workbook sets
`EXTERNAL_DATA_DIR`.

`04C` carries its own defence against being the wrong version, and it passed:

| Self-check | Value |
|---|---|
| `full_data_table13_row_match` | TRUE |
| `full_data_probability_max_abs_diff` | 4.44e-16 |
| `frozen_interaction_beta` | 0.740127 |
| `full_data_interaction_beta_abs_diff` | 1.11e-16 |
| `requested_iterations` / `successful_iterations` | 100 / 100 for both schemes |

The interaction Δβ is positive in 100% of iterations for both host-node and virus-TaxID
subsampling, which is what the manuscript's §3.4 statement requires. The 1.1e-16 agreement with the
frozen value is at floating-point precision, so this is the script that produced the published
number.

The 410-vs-340 distinction was checked directly rather than inferred: 410 is the count of virus
TaxIDs with a VIP-DB documented vector, 340 the count with sequence support, and 4,401 the count of
plant-virus **sequences** with any S0/S1/S2 evidence. The Data S1 `Contents & Notes!A18` correction
was made on exactly this basis and verified cell by cell: **1 changed cell out of the whole
workbook**, with zero changes in Table S1–S5.

### B4. Module 04 — genus–host network (`Results 3.5`, Fig. 5)

| Script | Verdict | Observed |
|---|---|---|
| `05A_network_main.R` | **PASS (runtime)** — but see the module-count block below | 269 virus genera; 197 host families; 1,902 weighted edges; **12 candidates** (identical set to the frozen run); modules / Q depend on the Louvain draw — 10 / 0.3446 on a re-run here, 12 / 0.3532 in the frozen run. The released `tables/03`–`06` are the **frozen** outputs, MD5-matched to the cluster, so Fig. 4 draws 12 modules and Q = 0.3532. |
| `05B_network_robustness.R` | **PASS (runtime)**, 11.2 s | 200 iterations at 80% resampling; all 12 candidates recur; median recurrence 1.00, min 0.43 |
| `05C_RdRP65_network_sensitivity.R` | **PASS (runtime)**, 2.8 s | RdRP65 network rebuilt and compared |
| `05D_projection_randomization.R` | **PASS (runtime)**, 19.4 s | transitivity 0.7274 vs random 0.6426; path length 1.7514 vs 1.5797; small-world σ 1.021; 999 randomisations |
| `05E_cross_genus_VC80_sensitivity.R` | **PASS (runtime)**, 0.9 s | candidate union 12; Jaccard 1.0; retention 1.0; node classification agreement 1.0 |

The thresholds were read from `00_config.R` rather than assumed: `FIG5_HOST_Q = 0.95`,
`FIG5_BRIDGE_Q = 0.75`, `FIG5_SUPPORT_MIN = 10L`, which in this network resolve to
`host_breadth_cutoff = 32`, `bridge_cutoff = 0.3456523`, `support_cutoff = 10`.

`05B`'s minimum recurrence of 0.43 for one candidate is the honest number and is reported as such;
it is not rounded up to 1.0 in the README. `05E`'s Jaccard of 1.0 means the candidate set is
insensitive to the genus-purity rule, which is the point of that sensitivity.

#### Network reproduction boundary — what holds and what does not

The module assignment comes from `cluster_louvain()` in igraph, called with the master seed. Louvain
is a **randomised, multi-start** heuristic: it returns the best modularity it happens to find, not a
unique partition. Everything the partition does not touch is bit-stable; the partition itself is not.

Measured on the shipped edge list (1,902 edges, identical MD5 on cluster and in the release):

| Observation | Result |
|---|---|
| 500 successive draws, seed varied | **500 distinct partitions**, module counts 9 / 10 / 11 / 12 / 13 |
| Module count distribution over those draws | 9: 4, **10: 90**, 11: 252, **12: 128**, 13: 26 |
| Re-run with the shipped seed here | 10 modules, Q = 0.344585 |
| Frozen cluster run | **12 modules, Q = 0.353235901330736** |
| Is the frozen Q reachable in 500 local draws? | **No** — absence over 500 draws, not a near-miss |
| Partition-independent quantities: betweenness / strength / degree, 269 genera | **max abs diff 3.7e-16 / 0 / 0** |
| Prioritised candidate set | **identical, 12 of 12** |

The release makes four specific statements and no more:

1. **Fig. 4 has been re-rendered from the frozen result tables and passes** (exit 0; the figure
   carries 12 modules, Q = 0.3532).
2. **Node-level metrics and the candidate set are consistent**: betweenness, strength and degree are
   partition-independent and match the frozen run to floating-point noise; the 12 prioritised genera
   are identical.
3. **The Louvain community partition has not been reproduced value-for-value**: a re-run of `05A` on
   the same edges and seed lands on a different (equally valid) optimum. The exact mechanism for why
   the frozen Q does not reappear — igraph build, BLAS, iteration order — has **not** been
   identified, and this release does not attribute it to any one cause.
4. **What is published is the frozen partition the paper actually uses**: the four module tables are
   the cluster's frozen outputs, MD5-matched (below). A reader re-running `05A` should expect a
   different module count; that is a property of the algorithm, not a defect in the data or code.

Consequence for the release: because a re-run cannot land on the published partition, the four tables
the figure reads are shipped as the **frozen cluster outputs**, MD5-matched:

| Released file | MD5 | Cluster MD5 |
|---|---|---|
| `analysis/04_genus_host_network/tables/03_network_module_membership.tsv` | `171e8d46ce9604ea7726958c78185dc6` | identical |
| `…/04_current_ICTV_genus_network_metrics.tsv` | `ea5fe9762bd0e6e5038e4af7bb21ca04` | identical |
| `…/05_candidate_thresholds.tsv` | `3592eccb354d1ac7f1946c7741c94542` | identical |
| `…/06_network_candidates.tsv` | `4b8b6f05d53df81fee4710bdd150c064` | identical |

Fig. 4 was re-rendered against them and exits 0, so the released figure carries 12 modules and
Q = 0.3532. `05A` remains runnable and remains the code that produced the frozen tables; running it
here simply lands on a different, equally valid Louvain optimum. The README states this boundary
rather than presenting 10 or 12 as a fixed deterministic output.

### B5. Module 05 — vector contextualisation (`Results 3.5`)

| Script | Verdict | Observed |
|---|---|---|
| `06A_candidate_external_context.R` | **PASS (runtime)**, 1.3 s | no GATE triggered; all 12 candidate genera had curated external evidence available |
| `06B_candidate_vector_support_adjusted_sensitivity.R` | **PASS (runtime)**, 0.4 s | M0 OR 53.25 (p 1.70e-4); M1 OR 2.94 (p 0.342); M2 OR 3.23; M3 OR 4.62; M4 OR 5.20; all converged |

The attenuation from OR 53.25 unadjusted to OR 2.94–5.20 after adjusting for evidence support is the
result the manuscript reports, and the script prints its own guardrail:

> These are sampling/support sensitivities, not causal estimates. Attenuation after support
> adjustment means the raw candidate-vector enrichment is partly associated with public evidence
> intensity and should not be presented as independent biological validation of candidate status.

That guardrail is preserved verbatim in the release.

### B6. Module 06 — TSWV case study (`Results 3.6`, Fig. 5 case study)

| Script | Verdict | Observed |
|---|---|---|
| `07A_candidate_species_decomposition.R` | **PASS (runtime)**, 4.1 s | 12 candidate genera; 1,081 current species; 285 species with a VIPDB vector record |
| `07N_segment_concordance_reassortment.R` | **PASS (runtime)**, 9.0 s | 120 matched isolates; 999 permutations; RF 172 / 176 / 184 (L-M / L-S / M-S); normalised RF 0.735 / 0.752 / 0.786; patristic p = 0.001 for all three |
| `07B`, `07F`, `07J`, `07K`, `07L`, `07O`, `07P` | **NOT VERIFIED** | depend on `pipeline_results/` intermediates (haplotype maps, GARD output, alignments) that are not part of the released data tree |

`07N` reports 202 / 154 / 75 incompatible high-support split pairs and a split Jaccard of 0.18–0.20
between segment pairs, together with its own guardrail that discordance is consistent with
reassortment but is not proof of a specific reassortment event. Both the statistic and the caveat
are kept.

The `07B/07F/07J/07K/07L/07O/07P` scripts belong to the *pipeline* layer (they consume
`pipeline_results/`, i.e. the rebuilt haplotype and GARD outputs), not to the released re-run layer.
They are published for transparency but are not expected to run on a clone. That is stated in
`docs/reproducibility.md` under level L1 rather than left implicit.

---

## C. Cross-cutting checks

### C1. Path hygiene

`grep` for hard-coded locations across the released tree returns nothing outside documentation that
deliberately quotes a cluster path when describing provenance. `00_paths.sh` resolves every root
relative to its own location, so the repository runs from any directory.

One real defect was found and fixed during this pass: inside `00_paths.sh`, the dependent external
paths were built as `${EXTERNAL_DATA_DIR}/file` while `EXTERNAL_DATA_DIR` defaults to empty, so
`VIPDB_RAW` became `/VIPDB_relation_raw.xlsx` — a path that "exists as a string" and therefore
defeated the `-z` test in `04A`. The variables are now only expanded when `EXTERNAL_DATA_DIR` is
non-empty, and otherwise left empty so callers detect "not supplied" cleanly. Verified after the fix:

```
EXTERNAL=[]   VIPDB_RAW=[]   PFAM=[]
```

This is exactly the class of bug the author asked to be eliminated: code that silently depends on
the author's machine layout.

A **second, worse instance of the same class** was found by the final end-to-end check — that is,
by unpacking the tarball into a clean directory and looking at what `00_paths.sh` actually exported,
rather than by reading it:

```
$ cd /tmp/verify/release && source analysis/00_config/00_paths.sh
ANALYSIS_ROOT=/tmp/verify          # wrong: one level above the repo
DATA_DIR=/tmp/data                 # therefore also wrong
```

The cause is `dirname "${BASH_SOURCE[0]}"`. Under bash it works. Under zsh, or under `sh -c`, or
under `. file` in dash, `BASH_SOURCE` does not exist, the expansion is empty, `dirname ""` is `.`,
and `${...}/..` becomes the **parent of the current working directory**. So the export was not
missing or obviously broken — it was a plausible path that happens to exist, which is precisely the
failure mode that survives review.

It is also worth being precise about how this was missed earlier: every previous run in this report
was launched as `bash -c 'source analysis/00_config/00_paths.sh; Rscript ...'`, and bash always had
`BASH_SOURCE`, so the defect could not appear. The wrong-root case only shows up in a shell the
reproduction instructions did not happen to use.

`00_paths.sh` now resolves its own location by trying `BASH_SOURCE`, then `$0`, then a relative
probe, and **validates each candidate** with `_is_repo_analysis_root()` — it must contain
`00_config/00_paths.sh`. If no candidate validates, the script prints an actionable error and
returns non-zero instead of exporting a wrong root. Verified across shells and working directories:

| Shell | Invocation | `ANALYSIS_ROOT` |
|---|---|---|
| bash | `source analysis/00_config/00_paths.sh` | `/workspace/release/analysis` |
| zsh | `source analysis/00_config/00_paths.sh` | `/workspace/release/analysis` |
| dash | `. analysis/00_config/00_paths.sh` | `/workspace/release/analysis` |
| bash | from `data/`, `source ../analysis/...` | `/workspace/release/analysis` |
| bash | from `/tmp`, absolute path | `/workspace/release/analysis` |

Regression after the fix, run from a **freshly unpacked tarball in `/tmp`**, not from the working
copy: `05A_network_main.R` printed `Fig5 main network PASS` with 269 genera / 1,902 edges / 12
candidates, and `06B_candidate_vector_support_adjusted_sensitivity.R` exited 0. The genera, edge and
candidate counts match §B4; the module count from that re-run does **not** match the frozen 12
(Louvain draw — see the module-count block in §B4). The released tables are the frozen ones, so the
figure is unaffected.

`00_preflight.R` had the **same category of staleness** and was caught by the same check. It still
looked for the twelve-family inputs as a flat set
(`FAMILY_TREE_DIR/<Family>.treefile`, `<Family>.nr95.membership.tsv`,
`<Family>.alignment_sequence_QC.tsv`), which stopped being true once the trees were laid out the way
`03L` actually reads them. On a clean unpack it therefore aborted on twelve files that were present
and correctly placed. It now checks the four real locations — `trees/<Family>/`,
`AA95_membership/`, `tree_tip_host_spectrum/`, and the per-family `06_AA95_alignment_QC/<Family>/`
(→ `07_Geminiviridae_Rep_conserved/AA95_alignment_QC/` for Geminiviridae) — and its failure message
reports **how many of each class are missing** rather than just listing paths. Verified both ways:
`Release preflight PASS` on a complete tree, and a per-class count plus the expected layout when one
treefile is deleted.

This is the argument for verifying from an unpacked archive instead of the working directory. Three
defects in this pass — the hard-coded `analysis_root`, the shell-dependent root resolution, and now
the stale preflight — were all invisible from the directory the code was written in.

### C2. Released inputs are internally consistent

Independently re-read with Python, not through the R code path:

| Check | Result |
|---|---|
| `sequence_analysis_input.tsv.gz` | 80,028 rows × 18 columns |
| `association_analysis_input.tsv.gz` | 83,034 rows × 30 columns |
| `vc80_final_id` unique count | 11,744 |
| distinct virus family labels | 76 |
| distinct virus TaxIDs | 5,979 |

These agree with the R-side counts, so the two independent reads corroborate each other.

### C3. Data S1

The corrected workbook (1,380,000 B,
SHA256 `eb4b5ef67fdcf4e0c806a2ac644248537ff0d4fff3196ab955ac6ba673c840ce`) was compared against the
original **cell by cell** across all sheets: **1 differing cell total**, `Contents & Notes!A18`,
with the unit correction from "4,401 TaxIDs" to "4,401 plant-virus sequences". Table S1–S5 show zero
changes. `supplementary/README.md` declares the `.xlsx` to be the numeric authority and the TSVs to
be convenience re-exports.

### C4. Manifest discipline

`docs/FROZEN_DEPENDENCY_CLOSURE.tsv` (48 rows) records, per script, its role, frozen inputs, frozen
outputs, manuscript anchor and a trace status. Two statuses exist:
`TRACED_FROM_CODE` (confirmed against code present in this snapshot) and `PENDING_CLUSTER_ACCESS`
(not retrievable). No row is marked traced on the strength of a filename, a timestamp or a hash
alone. `RELEASE_DIFF_AUDIT.tsv` (52 rows, 10 columns, including `change_type`) records the same
decision from the opposite direction — old path versus final path. Every row now carries a
`DONE`-prefixed status note or a stated reason why it cannot, so no row is left as a silent
`PENDING`.

### C5. Supplementary tables rebuilt from the release's own outputs

Two groups of supplementary tables did not exist in a final form and were rebuilt here. Unlike
Table S3A–H, which are the frozen submission outputs, these are **generated**: the check is that
the generator output is traceable to a script in this repository and that every value matches that
script's output.

| Table | Source script | Verification |
|---|---|---|
| S2E1 / S2E2 / S2E3 | `03L_FINAL_aggregate.R` over 12 × `03L_FINAL_run_one_family.R` | 12/12 families PASS; tip counts match `03K_D_tree_input_manifest.tsv` |
| S3I / S3J / S3K / S3L | `build_S3I_L.py` over 17 `04D_*` output tables | 17/17 blocks compared cell by cell, 0 differences |

The S3 sections are stacked rather than joined because their source blocks range from 9 to 41
columns. Each block carries `section_block` + `source_table` + its own header, and **no value is
recomputed, rounded or padded into a foreign schema**. The cost is that `read.table(header=TRUE)`
on the whole file is meaningless; the block map is in
`supplementary/Table_S3_Arthropod_vector/S3I_L_README.md` and the reading recipe is in
`supplementary/README.md`.

**Round-trip check on the packaged tarball.** The `04D` script was executed on a fresh extraction
of the release, then `build_S3I_L.py` was run over its output. All four S3 files came back
**byte-identical** to the released copies (SHA-256 equal for S3I, S3J, S3K and S3L), i.e. the
tables in the release are exactly what the released code produces from the released inputs:

```
04D_transmission_biology_extended_v3.R   EXIT=0   (4m 45s)
build_S3I_L.py  S3I: 1 block(s),  103 rows,  7 cols
                S3J: 5 block(s), 4691 rows, 43 cols
                S3K: 6 block(s),  351 rows, 20 cols
                S3L: 5 block(s),  875 rows, 43 cols
IDENTICAL  S3I_vector_family_by_virus_genus.tsv
IDENTICAL  S3J_five_group_host_range_models.tsv
IDENTICAL  S3K_country_breadth_and_evidence_models.tsv
IDENTICAL  S3L_within_mode_lineage_sensitivity.tsv
```

**Update 2026-10-27 (U25).** The K/L section letters above followed the *old*
file names; `Data_S1.xlsx` labels **K = within-mode lineage sensitivity,
L = country breadth**. The two files were renamed accordingly
(`S3K_within_mode_lineage_sensitivity.tsv`, `S3L_country_breadth_and_evidence_models.tsv`),
the builder was updated, and the renamed files were re-verified byte-identical
to the round-trip output above (same MD5s, names swapped). Contents unchanged.


### C6. Every figure rendered from a clean export

Each plotting wrapper was executed against a fresh `git archive` of the release, from the wrapper
itself rather than from the R script — that distinction is the point, because calling `Rscript`
directly had hidden six wrapper defects (U16).

| Figure | Wrapper | Result |
|---|---|---|
| Fig. 1 | `01_Fig1_data_framework_QC.R` (analysis script, no wrapper) | PASS — tables and QC |
| Fig. 2 | `run_fig2_publication_FINAL_v12.sh` | PASS — main + Fig. S2 at 3x4 and 2x6 |
| Fig. 3 | `run_fig3_publication_FINAL_v22.sh` | PASS — 34 files (main, 4 panels, supplement; PDF/PNG/TIFF) |
| Fig. 4 | `run_Fig4_publication_FINAL_v1_1.sh` | PASS — 34 files (main, panels A-D, Fig. S4) |
| Fig. S1 | `run_fig3_publication_plot_v2_11_1.sh` | PASS — main + supplementary atlases |
| Fig. 5 | `run_fig7_publication_FINAL_v1_3_6.sh` | PASS — 34 files (main, panels, 5-gene cladograms + segment trees; PDF/PNG/TIFF). Closed U19 by shipping the five frozen inputs from the cluster; `19_TSWV_selection_sites_long.tsv` matches the published gz byte-for-byte (MD5 `ec4253cd…`) |

Fig. 1 and Fig. 3 are worth stating explicitly: their rendered output was checked against the
manuscript's numbers, not just against exit code 0. Fig. 3 reports S0 = 5 / S1 = 60 / S2 = 338 /
union = 340 / Overall = 410, matching Data S1.

The two figures that needed generated inputs before they could run (Fig. 2 reads the `03L` summary
tables, Fig. S1 reads the `03G` sensitivity tables) were verified by running the generating scripts
first, in order, on the same export — `03F` (5m 28s, 12 families) then `03G`, and separately `03L`.
That is the full path a reader would follow.

### C7. Fig. S1 verified against the submitted PDF

The author supplied the submitted `Supplementary Fig1.pdf`. It identifies Fig. S1
as the **single-panel RdRP65 host-sharing decay figure** (98 × 75 mm) — not the
old four-panel composite this report had provisionally assumed.

| Check | Result |
|---|---|
| Page size, submitted PDF | 277 × 212 pt = 98 × 75 mm — exactly the frozen code's `save_triplet(width_mm = 98, height_mm = 75)` |
| Data, decay line, bubble scaling, axis categories | identical to the frozen code's internal copy of the same panel (`FigS3A_RdRP65_host_taxonomic_decay`) |
| y-axis title | submitted: "RdRP65 sharing probability (%)" (frozen copy: "Probability of sharing RdRP65 units (%)") |
| size-legend title | submitted: "Number of host pairs" (frozen copy: "Host-pair count") |
| Standalone entry | `plot/figS1/01_FigS1_publication_plot_v2_11_1_standalone.R` + `run_FigS1_standalone.sh`: frozen helper, frozen input `tables/04_RdRP65_host_sharing_by_distance.tsv`, submitted label wording; output `main/FigS1_RdRP65_host_sharing_decay.{pdf,png,tiff}`; EXIT=0 from a clone |
| Tree files | **none** — the submitted S1 reads no trees; the provisional four-tree requirement (former U23) was withdrawn and `S1_frozen_trees/` removed |

**This is the only figure verified directly against submitted artwork.**
Fig. 1–5 and Figs. S2–S5 were verified at runtime (exit codes, logged values,
panel counts) but have **not** been compared pixel-by-pixel against submitted
artwork; visual parity for those figures rests on identical frozen inputs and
unchanged plotting code.

---

## D. What this report does **not** claim

1. **It does not claim the release is complete.** The 12-family phylogenies, the `03K_D_*`
   manifests, all six figure directories, S2E1–E3 and S3I–L are now in place. Still open: a
   script-level diff pass over the 05–07 pipeline scripts against the cluster frozen copies **(now
   done for the TSWV layer: all 19 pipeline/analysis scripts + 5 plotting files diffed against the
   frozen cluster copies over an SSH tunnel; 19 of 23 byte-identical after normalising the
   `Fig7_TSWV_case_study` → `06_TSWV_evolution` rename, the rest differ only in path resolution and
   missing-input handling — see `RELEASE_DIFF_AUDIT.tsv`)**, and
   `20_three_source_sequence_master_FINAL.tsv`, which Fig. 1 map labelling and
   `03K_D_finalize_inputs.py` would prefer but which is not redistributed. See
   `UNRESOLVED_ITEMS.md`.
2. **It does not claim every released script is the cluster-frozen version.** Those re-verified here
   are marked PASS (runtime). `03F`/`03G`, the `03E` intermediate stage and the `07B`–`07P`
   pipeline layer are `PASS (static)` or `NOT VERIFIED` for the reasons given per row; where a
   script is superseded (e.g. `03F`/`03G` by the `03L` chain) the report says so instead of
   quietly dropping the row.
3. **It does not claim one-command reproduction from raw archives.** Level L1 in
   `docs/reproducibility.md` is explicitly a manual, documented procedure.
4. **It does not claim any third-party data is re-licensable.** Upstream licences are verified
   per source in `docs/EXTERNAL_INPUTS.md` §5 (checked 2026-10-27, with URLs); embedded
   Virus-Host DB records remain CC BY-NC-SA 4.0, outside the package's CC BY 4.0 grant.
5. **It does not claim scientific correctness beyond reproduction.** Reproducing a number is not
   evidence that the model is appropriate; that judgement belongs to the manuscript and its
   reviewers.

---

## E. Bottom line

| Layer | Result |
|---|---|
| Modules executed end to end and matched to frozen anchors | 01, 02 (`03A`/`03B`/`03C`/`03L`×12/`03L_aggregate`), 03 (`04A`/`04B`/`04C`/`04D` + specificity control), 04 (`05A`–`05E`), 05 (`06A`/`06B`), 06 (`07A`/`07N`) |
| Scripts total executed successfully | 40 (incl. 12 per-family phylogeny runs and the aggregate step) |
| Hard numeric anchors reproduced exactly | 4 of 4 (β, OR, RdRP65 β, RdRP65 OR) |
| `NOT_IN_PAPER` tables reproduced exactly | 7 of 7 (specificity control) |
| 12-family phylogenies rebuilt | 12 of 12, 6,010 tips, 11 tables each |
| Supplementary tables rebuilt and verified | S2E1–E3 from `03L`; S3I–L restacked from `04D`, 17 of 17 blocks verbatim, and byte-identical on re-run |
| Integrity | 435 SHA-256 entries, 435 OK / 0 FAILED on a fresh extraction of the tarball |
| Syntax | 124 scripts (62 R, 38 py, 24 sh), 0 parse / compile / `bash -n` errors |
| Path hygiene | `setwd(` 0, `/Users/` 0, `mxrconda` 0 in executable code; the only cluster paths left are in 3 plotting scripts marked `ARCHIVAL ONLY` |
| Figures rendered from a clean clone | **6 of 6 pass.** Fig. 1, 2, 3, 4, 5 and the standalone Fig. S1 all exit 0. Fig. S1 is verified against the submitted PDF (C7) |
| Plotting wrappers executed as shipped | **6 of 6 pass** (Fig. 2, 3, 4, 5 wrappers + Fig. 1 direct + Fig. S1 standalone wrapper) |
| Twelve-family block rebuildable from a clone | **yes** — `03L_FINAL_aggregate.R` runs exit 0 with zero diff, after the per-family tables were un-ignored |
| Modules fully verifiable | 01, 02 (incl. `03L` ×12), 03, 04, 05, 06 |
| Release status | **READY FOR PUBLIC RELEASE.** Fig. 1–5 and Fig. S1 all render from a clone; Fig. S1 is verified against the submitted `Supplementary Fig1.pdf`; the twelve-family block rebuilds with zero diff; network data is frozen and MD5-matched; the 05–07 script diff pass is done; no open items remain (U23 withdrawn — the submitted S1 reads no tree files). See `PUBLIC_RELEASE_SIGNOFF.md`. |

