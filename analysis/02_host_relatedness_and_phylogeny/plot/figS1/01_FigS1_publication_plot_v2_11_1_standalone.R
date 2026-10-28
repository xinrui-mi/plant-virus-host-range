#!/usr/bin/env Rscript
# ============================================================================
# STANDALONE Supplementary Fig. S1 plotting entry
#
# The submitted Supplementary Fig. S1 is the SINGLE-PANEL RdRP65 host-sharing
# decay plot (98 x 75 mm): orange decay line over host taxonomic distance
# (same species -> same kingdom), bubble size = number of host pairs.
#
# It is drawn with the frozen plotting logic of
# 03_Fig3_publication_plot_v2_11_1.R (cluster run Fig3_publication_plot_v2_11_1):
# same helper `make_distance_decay_plot()`, same frozen input table, same page
# size. The only differences from the frozen code's internal copy (which it
# saved as FigS3A_RdRP65_host_taxonomic_decay) are the two text labels that the
# submitted figure uses:
#   y-axis title : "RdRP65 sharing probability (%)"   (frozen copy: "Probability
#                  of sharing RdRP65 units (%)")
#   size legend  : "Number of host pairs"             (frozen copy: "Host-pair
#                  count")
# Data, geometry, colours and bubble scaling are identical.
#
# Required frozen input (tables/):
#   04_RdRP65_host_sharing_by_distance.tsv
#
# No phylogeny is used: the submitted S1 reads no tree files.
#
# Run:  Rscript 01_FigS1_publication_plot_v2_11_1_standalone.R
# or:   bash run_FigS1_standalone.sh
# Output: publication_figures_v2_11_1/main/FigS1_RdRP65_host_sharing_decay.{pdf,png,tiff}
# ============================================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

source(file.path(SCRIPT_DIR, "00_paths.R"))
source(file.path(SCRIPT_DIR, "00_check_packages.R"))
source(file.path(SCRIPT_DIR, "00_plot_theme.R"))
source(file.path(SCRIPT_DIR, "00_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

log_file <- file.path(OUT_LOG, "FigS1_publication_plot_v2_11_1.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Fig. S1 standalone plotting started: ", format(Sys.time()), "\n", sep = "")
cat("FIG3_ROOT: ", FIG3_ROOT, "\n", sep = "")
cat("Resolved figure font: ", FONT_FAMILY, "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")

# ============================================================
# Frozen input: RdRP65 host-sharing decay
# ============================================================
f04 <- file.path(TABLE_DIR, "04_RdRP65_host_sharing_by_distance.tsv")

if (!file.exists(f04)) {
  stop("Required frozen Fig. S1 input is missing:\n  ", f04)
}

input_manifest <- data.table(
  role = "Fig. S1 RdRP65 host-sharing decay",
  file = f04,
  exists = TRUE
)
fwrite(input_manifest, file.path(OUT_LOG, "FigS1_v2_11_1_input_manifest.tsv"), sep = "\t")

s <- read_tsv(f04)
cat("RdRP65 columns: ", paste(names(s), collapse = ", "), "\n", sep = "")

# ============================================================
# Fig. S1 — single panel, labels as in the submitted figure
# ============================================================
pS1 <- make_distance_decay_plot(
  s,
  line_colour = COL_ORANGE,
  y_title = "RdRP65 sharing probability (%)",
  show_pair_size = TRUE
)

# Match the submitted figure's legend title (frozen helper sets the size
# scale name to "Host-pair count"; the submitted S1 says "Number of host
# pairs"). The scale name must be overridden via guides(), not labels().
pS1 <- pS1 + guides(size = guide_legend(title = "Number of host pairs"))

save_triplet(
  pS1,
  "FigS1_RdRP65_host_sharing_decay",
  OUT_MAIN,
  width_mm = 98,
  height_mm = 75
)

cat("Fig. S1 standalone plotting completed: ", format(Sys.time()), "\n", sep = "")
cat("Main output: ", OUT_MAIN, "\n", sep = "")
