Fig2_publication_FINAL_v10
=========================

Purpose
-------
Publication-visualization update only. This package does NOT rerun any model,
permutation, marker extraction, alignment, or phylogenetic inference.

Main Fig.2 V3
-------------
A. Global host-taxonomic-distance decay: retained from the accepted V2.1 design.
B. Degree-preserving null distribution: retained from the accepted V2.1 design.
C. Integrated family-level summary:
   - left: family-specific logistic beta + frozen 95% CI for all 24 families;
   - blue points/intervals: the final 12-family phylogeny set;
   - orange dashed line: global beta;
   - right: observed Spearman rho + pair-sampling interval for the 12 eligible families;
   - filled/open rho points: support-stratified FDR significance.
   Both halves share the same 24-family row order, so the old C/D become one panel.

Supplementary Fig.2 V3
----------------------
The final 12 family trees are plotted as complete rectangular phylograms.
Endpoint/terminal-branch colour = host order annotation.
Endpoint size = host family breadth.
Point-size range is deliberately smaller than the old single-tree main-panel design.

Two atlas layouts are exported for visual comparison:
1. FigS2_12family_phylogenetic_atlas_2x6_v3.pdf
2. FigS2_12family_phylogenetic_atlas_3x4_v3.pdf

Trait mapping is strict
-----------------------
The script scans final 03K/03L TSV outputs for each family and selects per-tip
trait sources by actual overlap with the FINAL tree tip labels.
It requires:
- a scalar host order field (dominant_host_order / host_order / top_host_order), and
- host family breadth (host_family_n / host_family_breadth / equivalent).

It will NOT silently substitute a plural host-order set or draw an all-white atlas.
The exact source file, ID column, trait column, and mapped coverage are written to:
logs/FigS2_FINAL_v3_trait_source_audit.tsv

Run
---
conda activate mxrconda
bash Fig2_publication_FINAL_v10/run_fig2_publication_FINAL_v3_3.sh

Output root
-----------
<FIG3_ROOT>/publication_figures_FINAL_Fig2_v10/

Primary files
-------------
main/Fig2_main_publication_FINAL_v3.pdf
supplementary/FigS2_12family_phylogenetic_atlas_2x6_v3.pdf
supplementary/FigS2_12family_phylogenetic_atlas_3x4_v3.pdf
logs/FigS2_FINAL_v3_trait_source_audit.tsv


V3.1 corrections
----------------
1. The frozen 24-family decay-model set and final 12-family phylogeny set overlap in 11 families, not 12.
   Mitoviridae is present in the final phylogeny set but absent from the frozen 24-family beta table.
   The integrated panel therefore uses the 25-family union. Mitoviridae has a blank beta row and a rho estimate.
2. Tree host-order colour no longer requires a scalar dominant_host_order column.
   Priority is: existing scalar host-order annotation; otherwise an existing host_orders/order-set result is used.
   For a set containing exactly one host order, that order is coloured. Multi-order sets are mapped to Other.
   No dominant order is invented from an unordered set.
3. Endpoint-circle size range is reduced for the atlas.
4. Trait source and mapping method are written to FigS2_FINAL_v3_3_trait_source_audit.tsv.


V3.2 exact final-tree annotation fix
------------------------------------
The final 03K tree_tip_host_spectrum tables store representative_sequence_id,
host_taxids, host_family_n and host_families, but not host_order. V3.2 therefore
uses host_family_n directly for endpoint size and maps the already-frozen host_taxids
to host_order using only the frozen Stage-20 association master taxonomy fields.
This is a plotting annotation lookup only; no statistical analysis, permutation,
clustering, marker extraction or tree inference is rerun. Tips spanning >1 host order
are displayed as Other (grey), rather than assigning an arbitrary dominant order.

V3.3 resolver fix
- Geminiviridae host-spectrum is explicitly read from 07_Geminiviridae_Rep_conserved, matching the frozen 03K-D manifest.
- The other 11 families are explicitly read from their final 06_AA95_alignment_QC directories.
- Superseded duplicate Geminiviridae host-spectrum files are ignored.

V4 revisions
- Panel C and Panel D are vertically stacked: 24-family beta on top; 12-family rho below in the inherited family order.
- Long panel-C caption is removed.
- Supplementary Fig.2 now prioritizes the 3x4 atlas layout, removes repeated rho values from titles, uses a shared Top-8 host-order legend, and reduces endpoint point sizes.

V5 revisions
- Panel C and Panel D return to a side-by-side structure. Panel C spans the full 25-row union (13 beta-only families on top; 12 phylogeny-panel families on bottom), whereas Panel D occupies only the lower shared-height block and aligns row-by-row with those 12 families.
- The in-figure D legend retains the title "Support-stratified FDR" and sits in the upper-right blank block.
- Supplementary Fig.2 keeps the 3x4 atlas as the primary layout and tightens the host-order/host-family-breadth legend with a higher-contrast shared palette.

V6: C/D redesigned as independent aligned forest plots sharing a 25-family order.

V7 revisions
- Main Fig.2 family-level section is redesigned as one integrated three-column plot: beta forest (left), family labels (center), rho forest (right).
- The family axis is shared visually through the central label column, avoiding the separated D-panel layout from v6.
- Supplementary Fig.S2 is carried forward unchanged from v6 except for output filenames.

V8 revisions
- Main Fig.2 family-level section uses a beta forest plot with a compact rho bar-annotation strip aligned by family rows.
- The rho strip encodes observed rho as horizontal bars and support-stratified significance as filled/open endpoint dots.
- Supplementary Fig.S2 is carried forward from v7 except for output filenames.

V10 revisions
- Main Fig.2 family-level section is redesigned as a side-by-side paired display: left beta forest for 24 family-specific decay coefficients; right rho dot plot aligned on the same 25-family order.
- In panel C, families included in the phylogeny analysis are highlighted in blue; other beta-only families are shown as open circles with grey intervals.
- In panel D, only the 12 families with phylogenetic analysis are plotted, using filled/open dots to denote support-stratified significance, with empty rows preserved for the remaining families to maintain alignment.
- Supplementary Fig.S2 is carried forward from v8 except for output filenames.
