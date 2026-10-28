# Resolve the repository layout from this script's own location.
# Released path: <repo>/analysis/<module>/plot/<fig>/00_paths.R
local({
  this <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE)), error = function(e) NA_character_)
  if (is.na(this)) {
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) this <- dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = TRUE))
  }
  if (!is.na(this)) {
    assign("RELEASE_FIG_DIR", this, envir = .GlobalEnv)
    assign("RELEASE_REPO_ROOT", normalizePath(file.path(this, "..", "..", "..", ".."), winslash = "/", mustWork = FALSE), envir = .GlobalEnv)
    assign("RELEASE_ANALYSIS_ROOT", file.path(get("RELEASE_REPO_ROOT", envir = .GlobalEnv), "analysis"), envir = .GlobalEnv)
  }
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = if (exists("RELEASE_REPO_ROOT")) RELEASE_REPO_ROOT else "")
ANALYSIS_ROOT <- file.path(PROJECT_ROOT, "analysis")
FIG7_ROOT <- Sys.getenv("FIG7_ROOT", unset = file.path(ANALYSIS_ROOT, "06_TSWV_evolution"))

# Frozen Fig.7 analytical inputs.
GENE_SUMMARY <- file.path(FIG7_ROOT, "tables", "06_TSWV_five_gene_evolutionary_input_summary.tsv")
SEGMENT_SUMMARY <- file.path(FIG7_ROOT, "tables", "13_TSWV_segment_concordance_summary.tsv")
GARD_STATUS <- file.path(FIG7_ROOT, "tables", "16_TSWV_GARD_final_status.tsv")
SELECTION_SUMMARY <- file.path(FIG7_ROOT, "tables", "18_TSWV_GARD_selection_summary.tsv")
SELECTION_SITES <- file.path(FIG7_ROOT, "tables", "19_TSWV_selection_sites_long.tsv")

GENE_TREE_ROOT <- file.path(
  FIG7_ROOT, "pipeline_results", "03_alignment_and_phylogeny"
)

SEGMENT_TREE_ROOT <- file.path(
  FIG7_ROOT, "pipeline_results", "06_segment_concordance"
)

PUB_ROOT <- file.path(FIG7_ROOT, "publication_figures_FINAL_v1_3_6")
OUT_MAIN <- file.path(PUB_ROOT, "main")
OUT_SUPP <- file.path(PUB_ROOT, "supplementary")
OUT_PANEL <- file.path(PUB_ROOT, "panels")
OUT_LOG <- file.path(PUB_ROOT, "logs")

for (d in c(PUB_ROOT, OUT_MAIN, OUT_SUPP, OUT_PANEL, OUT_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Deterministic root guard.
project_norm <- normalizePath(PROJECT_ROOT, mustWork = FALSE)
for (p in c(
  FIG7_ROOT, GENE_SUMMARY, SEGMENT_SUMMARY,
  GARD_STATUS, SELECTION_SUMMARY, SELECTION_SITES,
  GENE_TREE_ROOT, SEGMENT_TREE_ROOT, PUB_ROOT
)) {
  pp <- normalizePath(
    if (file.exists(p) || dir.exists(p)) p else dirname(p),
    mustWork = FALSE
  )
  if (!startsWith(pp, project_norm)) {
    stop("Path guard failed; path outside final_three_source: ", p)
  }
}
