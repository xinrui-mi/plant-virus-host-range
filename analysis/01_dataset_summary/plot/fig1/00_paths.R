# ---------------------------------------------------------------------------
# Path configuration for the Fig. 1 plotting step.
#
# The plotting code reads the frozen analysis tables produced by
# `analysis/01_dataset_summary/analysis/01_Fig1_data_framework_QC.R` and
# `02_Fig2_global_landscape.R`; it does not recompute anything.
#
# Nothing here is machine specific. Paths resolve from ANALYSIS_ROOT, which is
# exported by `analysis/00_config/00_paths.sh`, and every value can be
# overridden with an environment variable if the tables live elsewhere.
# ---------------------------------------------------------------------------

# Resolve the repository layout from this script's own location, so the plotting
# step runs from a fresh clone without exporting anything first.
# Released path: <repo>/analysis/01_dataset_summary/plot/fig1/00_paths.R
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

# ANALYSIS_ROOT is still honoured if the user exports it (that is what
# analysis/00_config/00_paths.sh does); otherwise take it from above.
if (!nzchar(Sys.getenv("ANALYSIS_ROOT"))) {
  if (exists("RELEASE_ANALYSIS_ROOT")) {
    Sys.setenv(ANALYSIS_ROOT = RELEASE_ANALYSIS_ROOT)
  } else {
    stop("ANALYSIS_ROOT is not set and the repository layout could not be inferred.\n",
         "Export ANALYSIS_ROOT or source analysis/00_config/00_paths.sh first.")
  }
}
if (!nzchar(Sys.getenv("REPO_ROOT")) && exists("RELEASE_REPO_ROOT")) {
  Sys.setenv(REPO_ROOT = RELEASE_REPO_ROOT)
}

# Directory holding the frozen Fig. 2 input tables (01_virus_family_evidence.tsv, ...).
FIG2_ROOT <- Sys.getenv(
  "FIG2_ROOT",
  unset = file.path(Sys.getenv("ANALYSIS_ROOT"), "01_dataset_summary")
)
FIG2_ROOT <- normalizePath(FIG2_ROOT, winslash = "/", mustWork = FALSE)
TABLE_DIR <- file.path(FIG2_ROOT, "tables")

# Master sequence table. On the cluster this was the full three-source sequence
# master; the plotting code only uses it to derive a virus_family -> DNA/RNA
# lookup (derive_family_type_lookup) and to report the sequence-record count, so
# the released analysis input — which carries both fields at the same row count —
# is used here. Set SEQ_MASTER to override.
SEQ_MASTER <- Sys.getenv(
  "SEQ_MASTER",
  unset = file.path(Sys.getenv("REPO_ROOT"), "data", "analysis_input",
                    "sequence_analysis_input.tsv.gz")
)
SEQ_MASTER <- normalizePath(SEQ_MASTER, winslash = "/", mustWork = FALSE)

# Rendered output. Kept under a git-ignored directory at the repository root so a
# run never writes into the analysis tree.
OUT_ROOT  <- file.path(Sys.getenv("REPO_ROOT"), "figures", "fig1")
OUT_MAIN  <- file.path(OUT_ROOT, "main")
OUT_PANEL <- file.path(OUT_ROOT, "panels")
OUT_LOG   <- file.path(OUT_ROOT, "logs")

for (d in c(OUT_ROOT, OUT_MAIN, OUT_PANEL, OUT_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
