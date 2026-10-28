# `pipeline_results/` — how the released inputs are wired up

`03L_FINAL_run_one_family.R` expects its frozen inputs in a specific three-level layout that mirrors
the cluster. Rather than duplicate several megabytes of trees and alignments, the release provides
that layout with directories and one relative symlink.

```
analysis/02_host_relatedness_and_phylogeny/pipeline_results/03K_family_specific_markers/
├── 06_AA95_alignment_QC/
│   └── <Family>/<Family>.alignment_sequence_QC.tsv          # 12 families, real files
├── 07_Geminiviridae_Rep_conserved/
│   └── AA95_alignment_QC/
│       └── Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv
└── 08_final_12family_trees -> ../../../../data/analysis_data/phylogeny/twelve_family
```

Everything the trees side needs lives once, in `data/analysis_data/phylogeny/twelve_family/`
(see `data/README.md` for the file-by-file description). The symlink is **relative**, so the
repository works from any clone location and after any directory move; it is stored in the tarball as
a relative link, not as an absolute path.

## Why not just point the script at `data/`?

Because keeping the cluster layout means the published script is the same file that ran on the
cluster, with no release-only branch in the code. The only change the script needed was
`resolve_manifest_path()`, which lets manifest path columns be relative to the tree root instead of
absolute. That helper is documented at its definition in `03L_FINAL_run_one_family.R`.

## Rebuilding from scratch

If you would rather not rely on the symlink, set the environment variable and run directly:

```bash
export FIG3_ROOT=<repo>/analysis/02_host_relatedness_and_phylogeny
Rscript analysis/02_host_relatedness_and_phylogeny/analysis/03L_FINAL_run_one_family.R \
  <repo>/analysis Tombusviridae
```

The script resolves `03K_family_specific_markers/08_final_12family_trees` from `tree_root`, so as long
as that path leads to the tree directory the run proceeds.

## What gets written here

`03L_FINAL_run_one_family.R` writes to `03L_FINAL_12family_host_spectrum/<Family>/` — eleven tables
per family (main summary, permutation iterations, pair-sampling stability long and summary,
support-stratified summary and iterations, min-support sensitivity summary and iterations, alignment
filtering bias detail and summary, family QC). `03L_FINAL_aggregate.R` then writes seven pooled tables
one level up. All of it is regenerable and is excluded from the download archive by `.gitignore`
rules; the released copies of the three tables that appear in the paper live in
`supplementary/Table_S2_Host_relatedness/`.
