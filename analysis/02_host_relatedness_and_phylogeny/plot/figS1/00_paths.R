# Optional: an extra R library to prepend. Unset on a normal install.
R_LIB <- Sys.getenv("R_LIB", unset = "")
if (nzchar(R_LIB) && dir.exists(R_LIB)) .libPaths(unique(c(R_LIB, .libPaths())))

# Released layout: this script lives in <repo>/analysis/02_host_relatedness_and_phylogeny/plot/<fig>/,
# so the analysis module root is three levels up. On the original cluster this was
# <project>/three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint.
local({
  this <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE)), error = function(e) NA_character_)
  if (is.na(this)) {
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) this <- dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = TRUE))
  }
  if (!is.na(this)) assign("RELEASE_FIG_DIR", this, envir = .GlobalEnv)
})
DEFAULT_FIG3_ROOT <- if (exists("RELEASE_FIG_DIR")) normalizePath(file.path(RELEASE_FIG_DIR, "..", ".."), winslash = "/", mustWork = FALSE) else ""
FIG3_ROOT <- Sys.getenv("FIG3_ROOT", unset = DEFAULT_FIG3_ROOT)

if (!dir.exists(FIG3_ROOT)) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(wd, "tables")) && dir.exists(file.path(wd, "analysis"))) {
    FIG3_ROOT <- wd
  } else {
    stop(
      "Fig.3 root not found.\n",
      "Set FIG3_ROOT to the analysis module root (the directory holding tables/ and pipeline_results/).\n",
      "Expected default: ", DEFAULT_FIG3_ROOT
    )
  }
}

FIG3_ROOT <- normalizePath(FIG3_ROOT, winslash = "/", mustWork = TRUE)
TABLE_DIR    <- file.path(FIG3_ROOT, "tables")
WORK_DIR     <- file.path(FIG3_ROOT, "work")
PIPE_DIR     <- file.path(FIG3_ROOT, "pipeline")
ANALYSIS_DIR <- file.path(FIG3_ROOT, "analysis")

OUT_ROOT  <- file.path(FIG3_ROOT, "publication_figures_v2_11_1")
OUT_PANEL <- file.path(OUT_ROOT, "panels")
OUT_MAIN  <- file.path(OUT_ROOT, "main")
OUT_SUPP  <- file.path(OUT_ROOT, "supplementary")
OUT_LOG   <- file.path(OUT_ROOT, "logs")

for (d in c(OUT_ROOT, OUT_PANEL, OUT_MAIN, OUT_SUPP, OUT_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("FIG3_ROOT: ", FIG3_ROOT)
message("R libraries: ", paste(.libPaths(), collapse = " | "))
