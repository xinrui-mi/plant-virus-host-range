Fig7 publication FINAL v1.3.6
===================================

STATUS
------
Publication plotting package only.
All Fig.7 analytical outputs are read-only frozen inputs. This version does
not rerun or alter phylogenetic reconstruction, segment concordance, GARD,
BUSTED, MEME or FUBAR analyses.

FROZEN ANALYTICAL INPUTS
------------------------
Tables:
  06_TSWV_five_gene_evolutionary_input_summary.tsv
  13_TSWV_segment_concordance_summary.tsv
  16_TSWV_GARD_final_status.tsv
  18_TSWV_GARD_selection_summary.tsv
  19_TSWV_selection_sites_long.tsv

Formal coding-region trees:
  pipeline_results/03_alignment_and_phylogeny/
    RdRp/RdRp.codon_partitioned.treefile
    GnGc/GnGc.codon_partitioned.treefile
    NSm/NSm.codon_partitioned.treefile
    N/N.codon_partitioned.treefile
    NSs/NSs.codon_partitioned.treefile

Formal L/M/S segment trees:
  pipeline_results/06_segment_concordance/
    L/L.treefile
    M/M.treefile
    S/S.treefile

MAIN FIG.7
----------
A. TSWV genome organization and persistent propagative transmission.
   Genome organization is schematic and is not drawn to scale.
   The five coding regions are RdRp, glycoprotein precursor (Gn/Gc), NSm,
   N and NSs. NSm is labelled as Movement protein.

B. L/M/S segment-level phylogenetic discordance.
   Quantitative summaries:
     Patristic rho
     Normalized Robinson-Foulds (RF) distance
     High-support split Jaccard
   High-support split Jaccard is the Jaccard similarity between the sets of
   high-support bipartitions.
   Segment-tree discordance is interpreted as consistent with reassortment;
   it is not presented as direct proof of a reassortment event.

C. Five-coding-region recombination and selection landscape.
   Left:
     relative coding-region position (0-100%)
     MEME episodic sites
     FUBAR diversifying sites
     GARD breakpoints
   Right:
     MEME episodic selected-site fraction (%)
     FUBAR diversifying selected-site fraction (%)
     FUBAR purifying selected-site fraction (%)

   Selection thresholds used by the frozen analysis:
     MEME:  P <= 0.10
     FUBAR: posterior probability >= 0.90

   GARD evaluates within-coding-region recombination. It must not be described
   as confirming segment reassortment.

BUSTED
------
BUSTED is a coding-region-wide episodic diversifying-selection test. Its
formal results are retained in the analysis tables and should be reported in
Results / supplementary results rather than added as a redundant visual column
in Panel C.

SUPPLEMENTARY FIG.S7
--------------------
Five full topology-oriented coding-region cladograms:
  RdRp, Gn/Gc, NSm, N, NSs

No midpoint rooting is performed. Branch lengths are suppressed for the
cladogram presentation; displayed root position is not biologically
interpreted. Subtle UFBoot >=95 support nodes are retained.

SUPPLEMENTARY FIG.S8
--------------------
Full branch-length L/M/S phylogenies supporting main Panel B. Fig.S8 retains
branch lengths, substitutions/site scale bars, UFBoot >=95 support and detailed
topology-comparison metrics.

TYPOGRAPHY
----------
Frozen manuscript typography, identical to finalized Fig.2-Fig.6:
  Open Sans
  TEXT_TINY   = 5.7 pt
  TEXT_SMALL  = 6.8 pt
  TEXT_NORMAL = 8.2 pt
  TEXT_LARGE  = 9.2 pt

In v1.3.6, segment-metric row labels and the compact Panel-C key use 6.8 pt,
matching manuscript category / legend text.

PUBLICATION PLOTTING INTEGRITY CHECKS
-------------------------------------
These are plotting-input checks only; they do not modify the analysis.
The script requires:
  1. L_vs_M, L_vs_S and M_vs_S each occur exactly once;
  2. the three pairwise segment comparisons share one matched-isolate count;
  3. MEME episodic fraction remains within 0-10%;
  4. FUBAR diversifying fraction remains within 0-2%;
  5. FUBAR purifying fraction remains within 0-80%.

The percentage checks prevent silent clipping if frozen inputs are changed in
the future.

FINAL v1.3.6 POLISH
-------------------
1. Segment-metric row labels: 5.7 -> 6.8 pt.
2. Panel-C compact key text: 5.7 -> 6.8 pt.
3. NSm functional label: Movement -> Movement protein.
4. Removed the duplicated site_sel / bp_site data-preparation block.
5. Added publication plotting integrity checks for segment pairs, matched
   isolates and frozen percentage display ranges.
6. No analytical value, tree topology, GARD result or selection result changed.

RUN
---
cd /share/home/bgi_mixr/project/biye/0_newstart/0_rawseq/final_three_source/three_source_manuscript_analysis_v1_1/Fig7_TSWV_case_study

bash Fig7_publication_FINAL_v1_3_6/run_fig7_publication_FINAL_v1_3_6.sh

LSF
---
bash Fig7_publication_FINAL_v1_3_6/submit_fig7_publication_FINAL_v1_3_6.sh

OUTPUTS
-------
publication_figures_FINAL_v1_3_6/
  main/
  panels/
  supplementary/
  logs/

PDF / PNG 600 dpi / TIFF 600 dpi LZW.
