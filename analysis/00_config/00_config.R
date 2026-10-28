# ---------------------------------------------------------------------------
# Shared configuration for the released manuscript analysis code.
#
# The repository is self-contained: every dataset required by the six analysis
# modules is shipped under <repo>/data. Upstream raw databases that cannot be
# redistributed (NCBI Virus dumps, Virus-Host DB dumps, the raw VIPDB workbook,
# Pfam HMM, source GenBank records) are NOT part of the release and are
# referenced only through EXTERNAL_DATA_DIR; they are required solely by the
# phylogeny-rebuilding pipelines.
#
# From the repository root, run:
#
#     source analysis/00_config/00_paths.sh     # bash: exports the roots
#     Rscript analysis/<module>/analysis/<script>.R
#
# or in an interactive R session:
#
#     Sys.setenv(ANALYSIS_ROOT = normalizePath("analysis"))
#     source("analysis/<module>/analysis/<script>.R")
# ---------------------------------------------------------------------------

ANALYSIS_ROOT <- Sys.getenv("ANALYSIS_ROOT", unset = "")
if (!nzchar(ANALYSIS_ROOT)) {
  stop("ANALYSIS_ROOT is not set. From the repository root run:\n  source analysis/00_config/00_paths.sh")
}
ANALYSIS_ROOT <- normalizePath(ANALYSIS_ROOT, winslash = "/", mustWork = TRUE)
REPO_ROOT <- dirname(ANALYSIS_ROOT)
PROJECT_ROOT <- REPO_ROOT          # alias kept for compatibility with older module scripts
PIPE_ROOT <- file.path(REPO_ROOT, "data")

DATA_DIR <- file.path(REPO_ROOT, "data")
ANALYSIS_DATA <- file.path(DATA_DIR, "analysis_data")

# --- Released analysis-ready inputs (deterministic projections of the final master tables)
SEQ_FINAL <- file.path(DATA_DIR, "analysis_input", "sequence_analysis_input.tsv.gz")
ASSOC_FINAL <- file.path(DATA_DIR, "analysis_input", "association_analysis_input.tsv.gz")

# --- Released derived analysis inputs
ARTH_SEQ_FINAL <- file.path(ANALYSIS_DATA, "arthropod_sequence_evidence_FINAL.tsv.gz")
VIPDB_NORMALIZED <- file.path(ANALYSIS_DATA, "VIPDB_relation_normalized.tsv.gz")
VIPDB_VIRUS_LEVEL <- file.path(ANALYSIS_DATA, "VIPDB_virus_level_summary.tsv.gz")
CANDIDATE_EVIDENCE <- file.path(DATA_DIR, "external_evidence", "candidate_external_evidence.tsv.gz")

# --- Family-specific marker phylogenies (manuscript Fig. 2 / Table S2 section E)
# twelve_family/ holds the 12 final .treefile, .contree, .iqtree reports,
# the AA95 filtering alignments, the tree input manifest and the tip-to-host
# tables. See data/analysis_data/phylogeny/twelve_family/README.md.
PHYLO_DIR <- file.path(ANALYSIS_DATA, "phylogeny")
FAMILY_TREE_DIR <- file.path(PHYLO_DIR, "twelve_family")
TSWV_TREE_DIR <- file.path(PHYLO_DIR, "TSWV")
TSWV_LMS_META <- file.path(TSWV_TREE_DIR, "LMS_matched_isolate_metadata.tsv")

# --- Non-redistributable upstream inputs (optional; only needed to rebuild trees)
EXTERNAL_DATA_DIR <- Sys.getenv("EXTERNAL_DATA_DIR", unset = "")
ALL_THREE_FASTA <- file.path(EXTERNAL_DATA_DIR, "all_three_source.fasta")
RDRP_CORE_FASTA <- file.path(EXTERNAL_DATA_DIR, "08_RdRP_core.faa")
VIPDB_RAW <- file.path(EXTERNAL_DATA_DIR, "VIPDB_relation_raw.xlsx")
DNA_GENBANK <- file.path(EXTERNAL_DATA_DIR, "dna_allaccession.gb")
PFAM_HMM <- file.path(EXTERNAL_DATA_DIR, "Pfam-A.hmm")
TSWV_OLD_ROOT <- file.path(EXTERNAL_DATA_DIR, "TSWV_final_pipeline")
TSWV_SOURCE_GB <- file.path(TSWV_OLD_ROOT, "input", "sequence.gb")

# ---------------------------------------------------------------------------
# Analysis parameters (frozen; identical to the submitted analysis).
# Do not change these values if you intend to reproduce the published numbers.
# ---------------------------------------------------------------------------
SEED <- 20260813L                      # master seed for every stochastic step
FIG2_RAREFACTION_ITER <- 100L          # fixed-permutation evidence accumulation
FIG3_NULL_ITER <- as.integer(Sys.getenv("FIG3_NULL_ITER", "999"))            # degree-preserving edge-swap nulls
FIG3_RDRP_NULL_ITER <- as.integer(Sys.getenv("FIG3_RDRP_NULL_ITER", "999"))  # RdRP65 nulls
FIG3_SWAP_MULTIPLIER <- 5L
FIG3_MAX_SWAPS <- 300000L
FIG3_SENSITIVITY_ITER <- 50L           # 70% host/virus subsampling replicates
FIG3_SUBSAMPLE_FRACTION <- 0.70
FIG3_TREE_PERM_ITER <- 999L            # family-tree permutation nulls
FIG5_ROBUSTNESS_ITER <- as.integer(Sys.getenv("FIG5_ROBUSTNESS_ITER", "200"))
FIG5_ROBUSTNESS_FRACTION <- 0.80
FIG5_HOST_Q <- 0.95                    # host-breadth threshold quantile
FIG5_BRIDGE_Q <- 0.75                  # bridge-score threshold quantile
FIG5_SUPPORT_MIN <- 10L                # minimum QC-pass VC80 support
FIG5_PROJECTION_NULL_ITER <- as.integer(Sys.getenv("FIG5_PROJECTION_NULL_ITER", "999"))
TSWV_GARD_MAX_N <- as.integer(Sys.getenv("TSWV_GARD_MAX_N", "180"))
TSWV_GARD_MAX_BREAKPOINTS <- as.integer(Sys.getenv("TSWV_GARD_MAX_BREAKPOINTS", "20"))
TSWV_CDS_MIN_RELATIVE_LENGTH <- as.numeric(Sys.getenv("TSWV_CDS_MIN_RELATIVE_LENGTH", "0.90"))
TSWV_CDS_MAX_RELATIVE_LENGTH <- as.numeric(Sys.getenv("TSWV_CDS_MAX_RELATIVE_LENGTH", "1.10"))
NCORES <- as.integer(Sys.getenv("NCORES", "12"))
