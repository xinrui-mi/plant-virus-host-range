Fig3_transmission_publication_FINAL_v22
======================================

Purpose
-------
Publication plotting package for revised Fig.3 / Fig.S3 using frozen
Fig4_arthropod_transmission_ecology outputs and the final
04D_transmission_biology_extended_v3 definitions/results.

Main updates in v17
-------------------
1. Panel A keeps the no-box three-column structure, but the evidence-group
   titles are centred and moved closer to the row-label column.
2. Panel A displayed counts continue to be read from frozen files at runtime.
3. The S0/S1/S2 rows now use the same lighter blue as the old Fig.4 package,
   while S0-S2 union remains the darker blue.
4. The B/C row is rebalanced to a slightly narrower C panel while preserving a
   one-line Fig.3C legend.
5. Fig.S3A x-axis labels are restored to the frozen small-text scale.

Typography
----------
Exactly inherited from Fig4_publication_FINAL_v25/v26.

Run
---
conda activate mxrconda
bash Fig3_transmission_publication_FINAL_v22/run_fig3_publication_FINAL_v22.sh

Output
------
<Fig4_arthropod_transmission_ecology>/publication_figures_NEW_Fig3_FINAL_v22/
  main/
  panels/
  supplementary/
  logs/

V20 display fix
---------------
- Panel A information block widened and the right-aligned group-title anchor shifted right so both group titles are fully visible inside the shared shaded boxes.
- Typography, analytical content, colours, B/C/D panels, and supplementary panels are otherwise unchanged.


V21 note: Panel A switches from shaded grouping boxes to bracket-based group annotations. Brackets link only S0/S1/S2 and the four documented transmission modes; S0-S2 union and Overall remain outside the brackets.


V22 note: Panel A bracket annotations are shifted left, detailed row labels are shifted slightly right, and short leader lines connect group titles to bracket centers.
