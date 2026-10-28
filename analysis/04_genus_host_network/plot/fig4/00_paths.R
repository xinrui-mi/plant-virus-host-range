# Optional: an extra R library to prepend. Unset on a normal install.
R_LIB <- Sys.getenv("R_LIB", unset = "")
if (nzchar(R_LIB) && dir.exists(R_LIB)) .libPaths(unique(c(R_LIB, .libPaths())))

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

DEFAULT_NETWORK_ROOT <- if (exists("RELEASE_ANALYSIS_ROOT")) file.path(RELEASE_ANALYSIS_ROOT, "04_genus_host_network") else ""
DEFAULT_CONTEXT_ROOT <- if (exists("RELEASE_ANALYSIS_ROOT")) file.path(RELEASE_ANALYSIS_ROOT, "05_vector_contextualization") else ""
# The cluster used the full three-source association master. The plotting code
# only consumes it through make_taxonomy_maps(), i.e. it needs virus_genus,
# virus_family, host_family and host_order — all present in the released
# analysis input, so the release points here instead.
DEFAULT_ASSOC_MASTER <- if (exists("RELEASE_REPO_ROOT")) file.path(RELEASE_REPO_ROOT, "data", "analysis_input", "association_analysis_input.tsv.gz") else ""

NETWORK_ROOT <- Sys.getenv("NETWORK_ROOT", unset = Sys.getenv("FIG5_ROOT", unset = DEFAULT_NETWORK_ROOT))
CONTEXT_ROOT <- Sys.getenv("CONTEXT_ROOT", unset = Sys.getenv("FIG6_ROOT", unset = DEFAULT_CONTEXT_ROOT))
ASSOC_MASTER <- Sys.getenv("ASSOC_MASTER", unset = DEFAULT_ASSOC_MASTER)

if (!dir.exists(NETWORK_ROOT)) stop("Frozen network-analysis root not found: ", NETWORK_ROOT)
if (!dir.exists(CONTEXT_ROOT)) stop("Frozen candidate-contextualization root not found: ", CONTEXT_ROOT)
if (!file.exists(ASSOC_MASTER)) stop("Frozen three-source association master not found: ", ASSOC_MASTER)

NETWORK_ROOT <- normalizePath(NETWORK_ROOT, winslash = "/", mustWork = TRUE)
CONTEXT_ROOT <- normalizePath(CONTEXT_ROOT, winslash = "/", mustWork = TRUE)
ASSOC_MASTER <- normalizePath(ASSOC_MASTER, winslash = "/", mustWork = TRUE)

TABLE_DIR <- file.path(NETWORK_ROOT, "tables")
CONTEXT_TABLE_DIR <- file.path(CONTEXT_ROOT, "tables")
VECTOR_MODEL_TABLE <- file.path(CONTEXT_TABLE_DIR, "12_candidate_VIPDB_support_adjusted_models.tsv")
if (!file.exists(VECTOR_MODEL_TABLE)) stop("Frozen vector-adjustment model table not found: ", VECTOR_MODEL_TABLE)

OUT_ROOT  <- file.path(NETWORK_ROOT, "publication_figures_Fig4_FINAL_v1_1")
OUT_MAIN  <- file.path(OUT_ROOT, "main")
OUT_PANEL <- file.path(OUT_ROOT, "FINAL_individual_panels")
OUT_SUPP  <- file.path(OUT_ROOT, "supplementary")
OUT_LOG   <- file.path(OUT_ROOT, "logs")
for (d in c(OUT_ROOT, OUT_MAIN, OUT_PANEL, OUT_SUPP, OUT_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("NETWORK_ROOT: ", NETWORK_ROOT)
message("CONTEXT_ROOT: ", CONTEXT_ROOT)
message("ASSOC_MASTER: ", ASSOC_MASTER)
message("VECTOR_MODEL_TABLE: ", VECTOR_MODEL_TABLE)
message("R libraries: ", paste(.libPaths(), collapse = " | "))
