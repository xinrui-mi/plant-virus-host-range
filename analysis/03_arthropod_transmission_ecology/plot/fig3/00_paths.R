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
# On the cluster this was <project>/final_three_source; here the equivalent is
# the repository root, whose `analysis/` holds every module.
FINAL_THREE_SOURCE_ROOT <- Sys.getenv(
  "FINAL_THREE_SOURCE_ROOT",
  unset = if (exists("RELEASE_REPO_ROOT")) RELEASE_REPO_ROOT else ""
)
FINAL_THREE_SOURCE_ROOT <- normalizePath(FINAL_THREE_SOURCE_ROOT, mustWork = TRUE)

# Historical folder name retained. This package renders the NEW manuscript
# Fig.3 from the frozen arthropod / transmission-ecology analyses.
FIG3_ANALYSIS_ROOT <- Sys.getenv(
  "FIG3_ANALYSIS_ROOT",
  unset = if (exists("RELEASE_ANALYSIS_ROOT")) file.path(RELEASE_ANALYSIS_ROOT, "03_arthropod_transmission_ecology")
          else file.path(FINAL_THREE_SOURCE_ROOT, "analysis", "03_arthropod_transmission_ecology")
)
FIG3_ANALYSIS_ROOT <- normalizePath(FIG3_ANALYSIS_ROOT, mustWork = TRUE)

PIPELINE_RESULTS <- Sys.getenv(
  "FINAL_PIPELINE_RESULTS",
  unset = file.path(FINAL_THREE_SOURCE_ROOT, "data", "analysis_data")
)
PIPELINE_RESULTS <- normalizePath(PIPELINE_RESULTS, mustWork = TRUE)

TABLE_DIR <- normalizePath(file.path(FIG3_ANALYSIS_ROOT, "tables"), mustWork = TRUE)

assert_under_final_root <- function(path, mustWork = TRUE) {
  p <- normalizePath(path, mustWork = mustWork)
  root0 <- normalizePath(FINAL_THREE_SOURCE_ROOT, mustWork = TRUE)
  root <- paste0(root0, .Platform$file.sep)
  if (!(p == root0 || startsWith(p, root))) {
    stop(
      "REFUSED input outside current final_three_source root:\n",
      p, "\nAllowed root:\n", root0
    )
  }
  p
}

FIG3_ANALYSIS_ROOT <- assert_under_final_root(FIG3_ANALYSIS_ROOT)
PIPELINE_RESULTS <- assert_under_final_root(PIPELINE_RESULTS)
TABLE_DIR <- assert_under_final_root(TABLE_DIR)

OUT_DIR <- file.path(FIG3_ANALYSIS_ROOT, "publication_figures_NEW_Fig3_FINAL_v22")
MAIN_DIR <- file.path(OUT_DIR, "main")
PANEL_DIR <- file.path(OUT_DIR, "panels")
SUPP_DIR <- file.path(OUT_DIR, "supplementary")
LOG_DIR <- file.path(OUT_DIR, "logs")

for (d in c(OUT_DIR, MAIN_DIR, PANEL_DIR, SUPP_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("FINAL_THREE_SOURCE_ROOT: ", FINAL_THREE_SOURCE_ROOT)
message("PIPELINE_RESULTS: ", PIPELINE_RESULTS)
message("FIG3_ANALYSIS_ROOT: ", FIG3_ANALYSIS_ROOT)
message("TABLE_DIR: ", TABLE_DIR)
