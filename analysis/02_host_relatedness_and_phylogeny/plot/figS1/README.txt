Fig3_publication_plot_v2_11
===========================

Status
------
Final publication-layout candidate for Fig.3.

The scientific structure is frozen:
A | B
C full width
D full-width summary strip

No formal analysis is rerun.
No IQ-TREE reconstruction is rerun.
No long branch is automatically removed.

FONT
----
v2.11 switches the manuscript figure font to:

Open Sans

The cluster has already confirmed:
fc-match "Open Sans"
  -> OpenSans-Regular.ttf: "Open Sans" "Regular"

R systemfonts has also resolved the actual environment font file.

The run script therefore sets:

FIG_FONT="Open Sans"

00_plot_theme.R stops with an explicit error if Open Sans is requested
but fontconfig resolves to a different family.

Typography remains unified:
- panel labels: bold
- axis titles: bold
- axis text: regular
- legend titles: bold
- legend text: regular

Panel C
-------
- topology-only display compression unchanged
- Geminiviridae 1463 tips -> 310 display units
- real branch lengths retained
- midpoint rooting is for visualization only
- dominant host order = terminal branch colour
- Host-family breadth = endpoint point size
- actual Geminiviridae breadth range = 1–17
- main legend landmarks = 1, 5, 10, 15
- underlying mapping remains continuous through 17
- C title is shifted right to reserve a clean gutter for the C panel label
- folded-clade n labels remain removed
- publication-friendly branch-length scale bar retained

Panel D
-------
- full-width short summary strip
- panel D tag receives a larger inset box so it is not clipped
- point = observed rho
- line = pair-sampling 95% interval
- no rho text
- no P text

Supplementary trees
-------------------
All four formal trees remain COMPLETE and UNCOMPRESSED.

Each family:
- has its own Top 8 dominant host orders + Other
- retains a family-specific legend
- uses the same global colour for the same host order across families
- uses the same point-size mapping up to the global maximum breadth = 17
- uses regular publication-facing breadth legend landmarks rather than
  irregular 1/2/4/8/max sequences

With the current data this means approximately:
- Geminiviridae: 1, 5, 10, 15
- Caulimoviridae: 1, 2, 3, 4, 5
- Potyviridae: 1, 5, 10
- Tombusviridae: 1, 5, 10

Scale bars:
- use the nice 1–2–5 rule
- are moved farther below dense trees using a tip-count-dependent y offset
- Supplementary panels receive additional lower margin to avoid tree overlap

Long-branch QC
--------------
Generated:
logs/FigS3_terminal_branch_top1pct_QC.tsv

For the longest 1% of terminal branches in each family, the table reports:
- branch length
- median terminal branch length
- 95th percentile terminal branch length
- top-1% cutoff
- branch/median ratio
- available genus / host-order / breadth annotation

This is diagnostic only.

Run
---
conda activate mxrconda
bash Fig3_publication_plot_v2_11/run_fig3_publication_plot_v2_11.sh

Expected log
------------
Font requested: Open Sans
Resolved figure font: Open Sans

Outputs
-------
<Fig3_host_phylogenetic_constraint>/publication_figures_v2_11/

main/
panels/
supplementary/
logs/


v2.11.1 HOTFIX
---------------

Only two issues are changed relative to v2.11.

1. Long-branch QC bug fixed

The previous code used data.table::fifelse() with:
- a scalar TRUE/FALSE test
- a vector-valued yes branch

This caused:
"Length of 'yes' is ... but must be 1 or length of 'test' (1)."

It is now implemented with a scalar R if/else block:

if (is.finite(med_len) && med_len > 0) {
  out[, branch_to_median_ratio := edge_length / med_len]
} else {
  out[, branch_to_median_ratio := NA_real_]
}

No tree or statistic is changed.

2. Panel C label spacing

- the C tag inset box is widened so the Open Sans "C" glyph is not clipped
- the Geminiviridae panel title is shifted slightly farther right
- no phylogeny geometry, legend mapping, folding, scale bar, or panel layout is changed

Expected font:
Font requested: Open Sans
Resolved figure font: Open Sans

Run:
bash Fig3_publication_plot_v2_11_1/run_fig3_publication_plot_v2_11_1.sh
