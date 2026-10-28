# Preflight check for the released repository: verifies that every released
# input required by the six analysis modules is present and exposes the columns
# the final code depends on.
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_config.R"))
source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_helpers.R"))

required_files <- c(
  SEQ_FINAL, ASSOC_FINAL, ARTH_SEQ_FINAL,
  VIPDB_NORMALIZED, VIPDB_VIRUS_LEVEL, CANDIDATE_EVIDENCE, TSWV_LMS_META
)
FAMILIES_12 <- c(
  "Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
  "Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
  "Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
)

# The twelve-family inputs are NOT a flat set of files. Each is laid out where
# the script that consumes it expects to find it, and preflight has to check the
# real locations or it reports a false alarm:
#
#   trees/<Family>/<Family>.treefile            FAMILY_TREE_DIR/trees/...
#   AA95_membership/*.tsv                       FAMILY_TREE_DIR/AA95_membership/
#   tree_tip_host_spectrum/*.tsv                FAMILY_TREE_DIR/tree_tip_host_spectrum/
#   06_AA95_alignment_QC/<Family>/*.tsv         under 03K_family_specific_markers/,
#     or 07_Geminiviridae_Rep_conserved/...       because Geminiviridae uses the
#                                               Rep-conserved branch (03L lines 552-564)
fam_tree_file <- function(f) {
  file.path(FAMILY_TREE_DIR, "trees", f, paste0(f, ".treefile"))
}
fam_membership_file <- function(f) {
  # Geminiviridae's membership file carries the .Rep_conserved infix.
  cand <- c(
    file.path(FAMILY_TREE_DIR, "AA95_membership", paste0(f, ".AA95_membership.tsv")),
    file.path(FAMILY_TREE_DIR, "AA95_membership", paste0(f, ".Rep_conserved.AA95_membership.tsv"))
  )
  cand[which(file.exists(cand))[1]]
}
fam_tip_spectrum_file <- function(f) {
  file.path(FAMILY_TREE_DIR, "tree_tip_host_spectrum",
            paste0(f, ".tree_tip_host_spectrum.tsv"))
}
fam_alignment_qc_file <- function(f) {
  pipeline_root <- file.path(ANALYSIS_ROOT, "02_host_relatedness_and_phylogeny",
                             "pipeline_results", "03K_family_specific_markers")
  if (identical(f, "Geminiviridae")) {
    file.path(pipeline_root, "07_Geminiviridae_Rep_conserved", "AA95_alignment_QC",
              "Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv")
  } else {
    file.path(pipeline_root, "06_AA95_alignment_QC", f,
              paste0(f, ".alignment_sequence_QC.tsv"))
  }
}

fam_paths <- list(
  treefile  = vapply(FAMILIES_12, fam_tree_file, character(1), USE.NAMES = FALSE),
  membership = vapply(FAMILIES_12, fam_membership_file, character(1), USE.NAMES = FALSE),
  spectrum  = vapply(FAMILIES_12, fam_tip_spectrum_file, character(1), USE.NAMES = FALSE),
  alignqc   = vapply(FAMILIES_12, fam_alignment_qc_file, character(1), USE.NAMES = FALSE)
)

required_files <- c(required_files,
  unlist(fam_paths, use.names = FALSE),
  file.path(FAMILY_TREE_DIR, c("03K_D_tree_input_manifest.tsv",
                               "03K_D_final_eligibility_FROZEN.tsv",
                               "03K_D_tree_QC.tsv")),
  file.path(TSWV_TREE_DIR, c("L.treefile","M.treefile","S.treefile")))

missing <- required_files[!file.exists(required_files)]
if (length(missing)) {
  # Report per input class so the message says which of the four is short.
  short <- vapply(names(fam_paths), function(k) {
    sum(!file.exists(fam_paths[[k]]))
  }, integer(1))
  fam_hint <- if (any(short > 0)) {
    paste0("\n\nTwelve-family inputs missing (of 12): ",
           paste(sprintf("%s %d", names(short), short), collapse = ", "),
           "\nExpected layout under ", FAMILY_TREE_DIR, ":\n",
           "  trees/<Family>/<Family>.treefile\n",
           "  AA95_membership/<Family>[.Rep_conserved].AA95_membership.tsv\n",
           "  tree_tip_host_spectrum/<Family>.tree_tip_host_spectrum.tsv\n",
           "and the per-family alignment QC under\n  ",
           file.path("analysis/02_host_relatedness_and_phylogeny/pipeline_results",
                     "03K_family_specific_markers/06_AA95_alignment_QC/<Family>/"),
           " (Geminiviridae: 07_Geminiviridae_Rep_conserved/AA95_alignment_QC/).",
           "\nModules 01 and 03-06 do not need the trees; only 03L does.")
  } else ""
  stop("Missing released inputs:\n  ", paste(missing, collapse = "\n  "), fam_hint)
}


SEQ_COLS <- c("sequence_id","virus_realm","virus_kingdom","virus_phylum","virus_class","virus_order",
              "virus_family","virus_genus","virus_species","virus_taxid","virus_ictv_match_rank",
              "virus_type","genome","vc80_final_id","rdrp_operational_cluster_id",
              "arthropod_sequence_supported","source_label_4way","source_LandFlora")
ASSOC_COLS <- c("sequence_id","virus_realm","virus_kingdom","virus_phylum","virus_class","virus_order",
                "virus_family","virus_genus","virus_species","virus_taxid","virus_ictv_match_rank",
                "virus_type","genome","host_kingdom","host_phylum","host_class","host_order",
                "host_family","host_genus","host_species","host_taxid","country","continent",
                "vc80_final_id","rdrp_operational_cluster_id","arthropod_sequence_supported",
                "landflora_record_id","landflora_strict_plant_infection",
                "source_label_4way","source_LandFlora")

chk <- function(path, cols, label) {
  hdr <- names(fread_any(path, sep = "\t", nrows = 0, showProgress = FALSE))
  miss <- setdiff(cols, hdr)
  if (length(miss)) stop(label, " is missing required columns: ", paste(miss, collapse = ", "))
  TRUE
}
invisible(chk(SEQ_FINAL, SEQ_COLS, "sequence_analysis_input"))
invisible(chk(ASSOC_FINAL, ASSOC_COLS, "association_analysis_input"))

seq_n <- nrow(fread_any(SEQ_FINAL, sep = "\t", select = "sequence_id", showProgress = FALSE))
assoc_n <- nrow(fread_any(ASSOC_FINAL, sep = "\t", select = "sequence_id", showProgress = FALSE))
if (seq_n != 80028L) stop("sequence_analysis_input has ", seq_n, " rows; expected 80028")
if (assoc_n != 83034L) stop("association_analysis_input has ", assoc_n, " rows; expected 83034")

cat("Release preflight PASS\n")
cat("sequence records:", seq_n, "\n")
cat("association records:", assoc_n, "\n")
