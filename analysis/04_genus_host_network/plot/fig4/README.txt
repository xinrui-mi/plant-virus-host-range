Fig4_publication_FINAL_v1_1
========================

Purpose
-------
Final publication-only plotting package for manuscript Fig.4 and Fig.S4.
It merges the frozen current network figure package with the already-fitted
VIP-DB documented-vector adjustment panel from the former contextualization figure.

NO ANALYSIS IS RE-RUN
---------------------
- No network reconstruction.
- No Louvain rerun.
- No bridge-score recomputation.
- No candidate redefinition.
- No logistic regression refitting.
- All candidate cutoffs are read from frozen Fig5_network_prioritization tables.
- Fig.4D reads M0/M1/M3 ORs and 95% CIs directly from:
  Fig6_candidate_contextualization/tables/12_candidate_VIPDB_support_adjusted_models.tsv

Final main figure
-----------------
A: virus-family x host-order aggregated network overview (unchanged).
B: network module composition (unchanged).
C: host-family breadth x network bridge score candidate screen (unchanged).
D: documented vector relationship association for prioritized vs other genera:
   unadjusted versus adjustment for VC80 support or sequence count.

Fig.4D
------
- Facets: VC80 support; Sequence count.
- Rows: Unadjusted; Adjusted.
- Points/lines: frozen odds ratio and 95% CI.
- Dashed vertical line: OR = 1.
- X axis is log10.
- No values are copied manually from a previous PDF.

Supplementary Fig.S4
--------------------
A: 80% VC80 resampling recovery of the 12 prioritized genera.
B: threshold-perturbation robustness.
C: genus-pure VC80 sensitivity.
Formal virus genus names are italicized; generic group labels such as
"Prioritized genera" and "Other genera" remain roman.

Typography
----------
Open Sans; frozen manuscript scale:
TEXT_TINY   5.7 pt
TEXT_SMALL  6.8 pt
TEXT_NORMAL 8.2 pt

Important Panel C disclosure
----------------------------
Panel C retains the frozen monotonic DISPLAY-ONLY transform that compresses
host-family breadth values below the empirical 95th-percentile candidate cutoff.
Values at/above the cutoff are uncompressed. Candidate classification, thresholds,
and all numerical analyses use untransformed values. State this in the figure caption.

Run
---
bash run_Fig4_publication_FINAL_v1_1.sh

Cluster submit
--------------
bash submit_Fig4_publication_FINAL_v1_1.sh

Outputs
-------
<NETWORK_ROOT>/publication_figures_Fig4_FINAL_v1_1/
  main/
    Fig4_main_publication_FINAL_v1_1.pdf/.png/.tiff
  FINAL_individual_panels/
    Fig4A_network_overview.*
    Fig4B_network_module_composition.*
    Fig4C_host_breadth_bridge_candidates.*
    Fig4D_documented_vector_relationship_adjustment.*
  supplementary/
    FigS4_candidate_recovery.*
    FigS4_threshold_robustness.*
    FigS4_single_genus_VC80_sensitivity.*
    FigS4_network_robustness_FINAL_v1_1.*
  logs/
