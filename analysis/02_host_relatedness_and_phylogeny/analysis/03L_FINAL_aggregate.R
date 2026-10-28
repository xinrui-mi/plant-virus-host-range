#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1L) stop("Usage: 03L_FINAL_aggregate.R <analysis_root>")

# The argument is the analysis root itself (`analysis/` inside the released
# repository). The release keeps the same internal layout as the original
# cluster directory, so no extra path segment is inserted.
base <- file.path(normalizePath(args[1], mustWork=TRUE),
                  "02_host_relatedness_and_phylogeny")
outroot <- file.path(
  base,
  "pipeline_results",
  "03L_FINAL_12family_host_spectrum"
)

FAMS <- c(
  "Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
  "Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
  "Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
)

read_all <- function(filename) {
  out <- lapply(FAMS, function(fam) {
    p <- file.path(outroot, fam, filename)
    if (!file.exists(p)) stop("Missing final 03L file: ", p)
    fread(p, sep="\t")
  })
  rbindlist(out, fill=TRUE)
}

main <- read_all("01_main_summary.tsv")
main[, BH_q_greater := p.adjust(empirical_p_greater, method="BH")]
main[, BH_q_two_sided_abs := p.adjust(empirical_p_two_sided_abs, method="BH")]
setorder(main, -observed_spearman_rho_host_jaccard)
fwrite(
  main,
  file.path(outroot, "01_FINAL_12family_main_summary.tsv"),
  sep="\t"
)

stab <- read_all("04_pair_sampling_stability_summary.tsv")
setorder(stab, -rho_median)
fwrite(
  stab,
  file.path(outroot, "02_FINAL_pair_sampling_stability_summary.tsv"),
  sep="\t"
)

support <- read_all("05_support_stratified_summary.tsv")
support[, BH_q_greater := p.adjust(empirical_p_greater, method="BH")]
support[, BH_q_two_sided_abs := p.adjust(empirical_p_two_sided_abs, method="BH")]
setorder(support, -observed_rho)
fwrite(
  support,
  file.path(outroot, "03_FINAL_support_stratified_summary.tsv"),
  sep="\t"
)

thr <- read_all("07_min_support_sensitivity_summary.tsv")
thr[, BH_q_greater := p.adjust(empirical_p_greater, method="BH"),
    by=min_cluster_member_n]
thr[, BH_q_two_sided_abs := p.adjust(empirical_p_two_sided_abs, method="BH"),
    by=min_cluster_member_n]
setorder(thr, min_cluster_member_n, -observed_rho)
fwrite(
  thr,
  file.path(outroot, "04_FINAL_min_support_sensitivity_summary.tsv"),
  sep="\t"
)

bias <- read_all("10_alignment_filtering_bias_summary.tsv")
setorder(bias, virus_family, action, drop_reason)
fwrite(
  bias,
  file.path(outroot, "05_FINAL_alignment_filtering_bias_summary.tsv"),
  sep="\t"
)

qc <- read_all("11_family_QC.tsv")
fwrite(
  qc,
  file.path(outroot, "06_FINAL_family_QC.tsv"),
  sep="\t"
)

plot_ready <- merge(
  main[, .(
    virus_family,
    final_tree_tip_n,
    host_annotated_tip_n,
    host_family_universe_n,
    distinct_host_spectrum_n,
    total_possible_pair_n,
    tested_pair_n,
    pair_sampling_mode,
    observed_rho=observed_spearman_rho_host_jaccard,
    empirical_p_greater,
    BH_q_greater
  )],
  stab[, .(
    virus_family,
    pair_sampling_replicate_n=replicate_n,
    pair_sampling_rho_median=rho_median,
    pair_sampling_q025=rho_q025,
    pair_sampling_q975=rho_q975
  )],
  by="virus_family",
  all.x=TRUE
)

plot_ready <- merge(
  plot_ready,
  support[, .(
    virus_family,
    support_stratified_p=empirical_p_greater,
    support_stratified_BH_q=BH_q_greater
  )],
  by="virus_family",
  all.x=TRUE
)

plot_ready[, direction := fifelse(
  observed_rho > 0, "positive",
  fifelse(observed_rho < 0, "negative", "zero")
)]
plot_ready[, main_FDR_significant := is.finite(BH_q_greater) & BH_q_greater < 0.05]
plot_ready[, support_stratified_FDR_significant :=
             is.finite(support_stratified_BH_q) &
             support_stratified_BH_q < 0.05]

setorder(plot_ready, -observed_rho)
fwrite(
  plot_ready,
  file.path(outroot, "07_FINAL_plot_ready_family_summary.tsv"),
  sep="\t"
)

cat("03L FINAL aggregate PASS\n")
cat("family_n=", nrow(main), "\n", sep="")
cat(
  "positive_rho_n=",
  sum(main$observed_spearman_rho_host_jaccard > 0, na.rm=TRUE),
  "\n",
  sep=""
)
cat(
  "main_BH_q_lt_0.05_n=",
  sum(main$BH_q_greater < 0.05, na.rm=TRUE),
  "\n",
  sep=""
)
cat(
  "support_stratified_BH_q_lt_0.05_n=",
  sum(support$BH_q_greater < 0.05, na.rm=TRUE),
  "\n",
  sep=""
)

print(main)
