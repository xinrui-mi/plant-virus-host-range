#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L) {
  stop("Usage: 03L_FINAL_run_one_family.R <analysis_root> <virus_family>")
}

# The first argument is the analysis root itself (the directory holding
# 00_config/, 01_dataset_summary/, 02_host_relatedness_and_phylogeny/, ...),
# i.e. `analysis/` inside the released repository. On the original cluster this
# was <project_root>/three_source_manuscript_analysis_v1_1; the release uses the
# same internal layout, so no extra path segment is inserted.
analysis_root <- normalizePath(args[1], mustWork=TRUE)
fam <- args[2]

Sys.setenv(ANALYSIS_ROOT=analysis_root)
source(file.path(analysis_root, "00_config", "00_config.R"))

base <- file.path(analysis_root, "02_host_relatedness_and_phylogeny")
kroot <- file.path(base, "pipeline_results", "03K_family_specific_markers")
tree_root <- file.path(kroot, "08_final_12family_trees")
outroot <- file.path(base, "pipeline_results", "03L_FINAL_12family_host_spectrum")
outdir <- file.path(outroot, fam)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

FAMS <- c(
  "Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
  "Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
  "Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
)
if (!fam %in% FAMS) stop("Unknown family: ", fam)
fam_index <- match(fam, FAMS)

# ---------------------------
# Frozen analysis parameters
# ---------------------------
MAX_PAIRS <- 50000L
PERM_ITER <- 999L
PAIR_STABILITY_REPS <- 100L
SUPPORT_THRESHOLDS <- c(2L, 3L)
SUPPORT_BREAKS <- c(-Inf, 1, 2, 4, 9, Inf)
SUPPORT_LABELS <- c("1", "2", "3-4", "5-9", "10+")
SEED0 <- if (exists("SEED")) as.integer(SEED) else 20260813L
FAM_SEED <- SEED0 + 900000L + fam_index * 10000L

safe_cor <- function(x, y) {
  z <- suppressWarnings(cor(x, y, method="spearman", use="complete.obs"))
  as.numeric(z)
}

safe_quant <- function(x, p) {
  as.numeric(quantile(x, p, na.rm=TRUE, names=FALSE, type=7))
}

read_manifest_row <- function() {
  manifest_path <- file.path(tree_root, "03K_D_tree_input_manifest.tsv")
  if (!file.exists(manifest_path)) stop("Missing tree manifest: ", manifest_path)
  man <- fread(manifest_path, sep="\t")
  z <- man[virus_family == fam]
  if (nrow(z) != 1L) stop("Tree manifest row not unique for ", fam)
  z
}

# The released manifest stores its three file columns as paths relative to
# `tree_root` (the cluster original used absolute paths). Resolve them here so
# the script works both ways: an absolute path is used as-is, a relative one is
# anchored to `tree_root`.
resolve_manifest_path <- function(p) {
  p <- as.character(p)[1]
  if (!nzchar(p) || is.na(p)) return(p)
  if (file.exists(p)) return(p)
  alt <- file.path(tree_root, p)
  if (file.exists(alt)) return(alt)
  p
}

parse_host_set <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  z <- unique(trimws(strsplit(x, ";", fixed=TRUE)[[1]]))
  z[nzchar(z)]
}

make_all_pairs <- function(n) {
  which(upper.tri(matrix(FALSE, nrow=n, ncol=n)), arr.ind=TRUE)
}

choose_pairs <- function(pairmat, seed) {
  if (nrow(pairmat) <= MAX_PAIRS) {
    list(
      pairs=pairmat,
      mode="all_pairs",
      tested_n=nrow(pairmat)
    )
  } else {
    set.seed(seed)
    idx <- sample.int(nrow(pairmat), MAX_PAIRS, replace=FALSE)
    list(
      pairs=pairmat[idx,,drop=FALSE],
      mode="random_pair_sample_50000",
      tested_n=MAX_PAIRS
    )
  }
}

rho_jaccard <- function(pairmat, dmat, jmat) {
  safe_cor(dmat[pairmat], jmat[pairmat])
}

rho_breadth <- function(pairmat, dmat, breadth) {
  bd <- abs(breadth[pairmat[,1]] - breadth[pairmat[,2]])
  safe_cor(dmat[pairmat], bd)
}

# ---------------------------
# Final frozen inputs
# ---------------------------
mr <- read_manifest_row()
treefile <- file.path(tree_root, "trees", fam, paste0(fam, ".treefile"))
tipfile <- resolve_manifest_path(mr$tree_tip_host_spectrum[1])
memfile <- resolve_manifest_path(mr$AA95_membership[1])

for (p in c(treefile, tipfile, memfile)) {
  if (!file.exists(p)) stop("Missing frozen input: ", p)
}

tr <- read.tree(treefile)
if (is.null(tr$tip.label) || length(tr$tip.label) < 10L) {
  stop("Final tree has too few tips for ", fam)
}

tip <- fread(tipfile, sep="\t", quote="", fill=TRUE)
need_tip <- c(
  "representative_sequence_id","cluster_member_n",
  "host_annotated","host_families"
)
miss_tip <- setdiff(need_tip, names(tip))
if (length(miss_tip)) {
  stop("Missing tip-table columns for ", fam, ": ", paste(miss_tip, collapse=","))
}

tip[, representative_sequence_id := as.character(representative_sequence_id)]
tip[, cluster_member_n := as.integer(cluster_member_n)]
tip[, host_annotated := as.character(host_annotated)]
tip[, host_families := fifelse(is.na(host_families), "", as.character(host_families))]

# Only tree tips with at least one host-family annotation enter host-spectrum analysis.
eligible_tip <- tip[
  host_annotated == "YES" &
  nzchar(host_families) &
  representative_sequence_id %in% tr$tip.label
]

ordered_tips <- tr$tip.label[tr$tip.label %in% eligible_tip$representative_sequence_id]
if (length(ordered_tips) < 10L) {
  stop("Fewer than 10 host-annotated final tree tips for ", fam)
}

tip2 <- eligible_tip[match(ordered_tips, representative_sequence_id)]
if (anyNA(tip2$representative_sequence_id)) {
  stop("Tree-tip/host-spectrum mapping failed for ", fam)
}
if (!all(tip2$representative_sequence_id == ordered_tips)) {
  stop("Tree-tip order mismatch for ", fam)
}

tr2 <- keep.tip(tr, ordered_tips)

# ---------------------------
# Host-family Jaccard matrix
# ---------------------------
host_sets <- lapply(tip2$host_families, parse_host_set)
host_universe <- sort(unique(unlist(host_sets, use.names=FALSE)))
if (length(host_universe) < 2L) {
  stop("Fewer than two host families in final host-spectrum set for ", fam)
}

H <- matrix(
  0L,
  nrow=length(host_sets),
  ncol=length(host_universe),
  dimnames=list(ordered_tips, host_universe)
)
hidx <- setNames(seq_along(host_universe), host_universe)
for (i in seq_along(host_sets)) {
  if (length(host_sets[[i]])) {
    H[i, unname(hidx[host_sets[[i]]])] <- 1L
  }
}

breadth <- rowSums(H)
intersect_mat <- tcrossprod(H)
union_mat <- outer(breadth, breadth, "+") - intersect_mat

J <- matrix(NA_real_, nrow=nrow(H), ncol=nrow(H))
ok <- union_mat > 0
J[ok] <- 1 - intersect_mat[ok] / union_mat[ok]
diag(J) <- 0

# ---------------------------
# Patristic distance matrix
# ---------------------------
D <- cophenetic.phylo(tr2)
D <- D[ordered_tips, ordered_tips, drop=FALSE]

all_pairs <- make_all_pairs(length(ordered_tips))
total_pairs <- nrow(all_pairs)
if (total_pairs < 1L) stop("No valid tip pairs for ", fam)

chosen <- choose_pairs(all_pairs, FAM_SEED + 1L)
pp <- chosen$pairs

obs_rho <- rho_jaccard(pp, D, J)
obs_breadth <- rho_breadth(pp, D, breadth)
if (!is.finite(obs_rho)) {
  stop("Primary rho is not finite for ", fam)
}

# ---------------------------
# Main 999-permutation null
# ---------------------------
perm_rho <- rep(NA_real_, PERM_ITER)
perm_breadth <- rep(NA_real_, PERM_ITER)

set.seed(FAM_SEED + 100L)
for (it in seq_len(PERM_ITER)) {
  perm <- sample.int(length(ordered_tips))
  jp <- J[cbind(perm[pp[,1]], perm[pp[,2]])]
  bp <- abs(breadth[perm[pp[,1]]] - breadth[perm[pp[,2]]])
  dv <- D[pp]

  perm_rho[it] <- safe_cor(dv, jp)
  perm_breadth[it] <- safe_cor(dv, bp)
}

finite_perm <- is.finite(perm_rho)
if (!any(finite_perm)) stop("No finite permutation rho values for ", fam)

p_greater <- (sum(perm_rho[finite_perm] >= obs_rho) + 1) /
             (sum(finite_perm) + 1)
p_two <- (sum(abs(perm_rho[finite_perm]) >= abs(obs_rho)) + 1) /
         (sum(finite_perm) + 1)

finite_b <- is.finite(perm_breadth)
p_breadth_greater <- if (is.finite(obs_breadth) && any(finite_b)) {
  (sum(perm_breadth[finite_b] >= obs_breadth) + 1) /
  (sum(finite_b) + 1)
} else {
  NA_real_
}

main_summary <- data.table(
  virus_family=fam,
  final_tree_tip_n=length(tr$tip.label),
  host_annotated_tip_n=length(ordered_tips),
  host_family_universe_n=length(host_universe),
  distinct_host_spectrum_n=uniqueN(tip2$host_families),
  total_possible_pair_n=total_pairs,
  tested_pair_n=chosen$tested_n,
  pair_sampling_mode=chosen$mode,
  observed_spearman_rho_host_jaccard=obs_rho,
  null_median_rho=as.numeric(median(perm_rho[finite_perm])),
  null_q025=as.numeric(safe_quant(perm_rho[finite_perm], .025)),
  null_q975=as.numeric(safe_quant(perm_rho[finite_perm], .975)),
  empirical_p_greater=p_greater,
  empirical_p_two_sided_abs=p_two,
  observed_spearman_rho_breadth_difference=obs_breadth,
  breadth_null_median_rho=if(any(finite_b)) as.numeric(median(perm_breadth[finite_b])) else NA_real_,
  breadth_empirical_p_greater=p_breadth_greater,
  permutations=PERM_ITER
)

fwrite(main_summary, file.path(outdir, "01_main_summary.tsv"), sep="\t")
fwrite(
  data.table(
    virus_family=fam,
    iteration=seq_len(PERM_ITER),
    perm_rho_host_jaccard=perm_rho,
    perm_rho_breadth_difference=perm_breadth
  ),
  file.path(outdir, "02_main_permutation_iterations.tsv"),
  sep="\t"
)

# ---------------------------
# Pair-sampling stability
# ---------------------------
stab_n <- if (total_pairs > MAX_PAIRS) PAIR_STABILITY_REPS else 1L
stab_list <- vector("list", stab_n)

if (total_pairs > MAX_PAIRS) {
  for (rr in seq_len(PAIR_STABILITY_REPS)) {
    set.seed(FAM_SEED + 1000L + rr)
    idx <- sample.int(total_pairs, MAX_PAIRS, replace=FALSE)
    ppr <- all_pairs[idx,,drop=FALSE]
    stab_list[[rr]] <- data.table(
      virus_family=fam,
      replicate=rr,
      tested_pair_n=MAX_PAIRS,
      rho=rho_jaccard(ppr, D, J)
    )
  }
  stab_mode <- "random_pair_sample_50000"
} else {
  stab_list[[1]] <- data.table(
    virus_family=fam,
    replicate=1L,
    tested_pair_n=total_pairs,
    rho=obs_rho
  )
  stab_mode <- "all_pairs"
}

stab_long <- rbindlist(stab_list)
finite_stab <- stab_long[is.finite(rho)]
if (!nrow(finite_stab)) stop("No finite pair-sampling rho values for ", fam)

stab_summary <- data.table(
  virus_family=fam,
  total_possible_pair_n=total_pairs,
  tested_pair_n=if(total_pairs > MAX_PAIRS) MAX_PAIRS else total_pairs,
  sampling_mode=stab_mode,
  replicate_n=nrow(stab_long),
  rho_median=as.numeric(median(finite_stab$rho)),
  rho_mean=as.numeric(mean(finite_stab$rho)),
  rho_sd=if(nrow(finite_stab)>1L) as.numeric(sd(finite_stab$rho)) else NA_real_,
  rho_q025=as.numeric(safe_quant(finite_stab$rho, .025)),
  rho_q25=as.numeric(safe_quant(finite_stab$rho, .25)),
  rho_q75=as.numeric(safe_quant(finite_stab$rho, .75)),
  rho_q975=as.numeric(safe_quant(finite_stab$rho, .975)),
  rho_min=as.numeric(min(finite_stab$rho)),
  rho_max=as.numeric(max(finite_stab$rho))
)

fwrite(stab_long, file.path(outdir, "03_pair_sampling_stability_long.tsv"), sep="\t")
fwrite(stab_summary, file.path(outdir, "04_pair_sampling_stability_summary.tsv"), sep="\t")

# ---------------------------
# Support-stratified permutation
# ---------------------------
strata <- cut(
  tip2$cluster_member_n,
  breaks=SUPPORT_BREAKS,
  labels=SUPPORT_LABELS,
  right=TRUE,
  include.lowest=TRUE
)
strata <- droplevels(strata)

sperm_rho <- rep(NA_real_, PERM_ITER)
set.seed(FAM_SEED + 2000L)

for (it in seq_len(PERM_ITER)) {
  perm <- seq_along(ordered_tips)
  for (lev in levels(strata)) {
    ii <- which(strata == lev)
    if (length(ii) > 1L) {
      perm[ii] <- sample(ii, length(ii), replace=FALSE)
    }
  }

  jp <- J[cbind(perm[pp[,1]], perm[pp[,2]])]
  sperm_rho[it] <- safe_cor(D[pp], jp)
}

finite_sperm <- is.finite(sperm_rho)
if (!any(finite_sperm)) stop("No finite support-stratified permutation values for ", fam)

sp_greater <- (sum(sperm_rho[finite_sperm] >= obs_rho) + 1) /
              (sum(finite_sperm) + 1)
sp_two <- (sum(abs(sperm_rho[finite_sperm]) >= abs(obs_rho)) + 1) /
          (sum(finite_sperm) + 1)

support_summary <- data.table(
  virus_family=fam,
  host_annotated_tip_n=length(ordered_tips),
  support_stratum_n=nlevels(strata),
  support_strata_present=paste(levels(strata), collapse=";"),
  tested_pair_n=chosen$tested_n,
  observed_rho=obs_rho,
  null_median_rho=as.numeric(median(sperm_rho[finite_sperm])),
  null_q025=as.numeric(safe_quant(sperm_rho[finite_sperm], .025)),
  null_q975=as.numeric(safe_quant(sperm_rho[finite_sperm], .975)),
  empirical_p_greater=sp_greater,
  empirical_p_two_sided_abs=sp_two,
  permutations=PERM_ITER
)

fwrite(support_summary, file.path(outdir, "05_support_stratified_summary.tsv"), sep="\t")
fwrite(
  data.table(
    virus_family=fam,
    iteration=seq_len(PERM_ITER),
    perm_rho=sperm_rho
  ),
  file.path(outdir, "06_support_stratified_iterations.tsv"),
  sep="\t"
)

# ---------------------------
# >=2 / >=3 AA95 member support sensitivity
# ---------------------------
threshold_summary <- list()
threshold_iterations <- list()

for (thr in SUPPORT_THRESHOLDS) {
  keep <- which(tip2$cluster_member_n >= thr)

  if (length(keep) < 10L) {
    threshold_summary[[as.character(thr)]] <- data.table(
      virus_family=fam,
      min_cluster_member_n=thr,
      tip_n=length(keep),
      total_possible_pair_n=NA_real_,
      tested_pair_n=NA_integer_,
      sampling_mode=NA_character_,
      status="FAIL_TOO_FEW_TIPS",
      observed_rho=NA_real_,
      null_median_rho=NA_real_,
      null_q025=NA_real_,
      null_q975=NA_real_,
      empirical_p_greater=NA_real_,
      empirical_p_two_sided_abs=NA_real_,
      permutations=PERM_ITER
    )
    next
  }

  Dt <- D[keep, keep, drop=FALSE]
  Jt <- J[keep, keep, drop=FALSE]
  apt <- make_all_pairs(length(keep))
  cpt <- choose_pairs(apt, FAM_SEED + 3000L + thr)
  ppt <- cpt$pairs
  ob <- rho_jaccard(ppt, Dt, Jt)

  pr <- rep(NA_real_, PERM_ITER)
  set.seed(FAM_SEED + 4000L + thr)

  for (it in seq_len(PERM_ITER)) {
    perm <- sample.int(length(keep))
    jp <- Jt[cbind(perm[ppt[,1]], perm[ppt[,2]])]
    pr[it] <- safe_cor(Dt[ppt], jp)
  }

  fin <- is.finite(pr)
  if (!is.finite(ob) || !any(fin)) {
    threshold_summary[[as.character(thr)]] <- data.table(
      virus_family=fam,
      min_cluster_member_n=thr,
      tip_n=length(keep),
      total_possible_pair_n=nrow(apt),
      tested_pair_n=cpt$tested_n,
      sampling_mode=cpt$mode,
      status="FAIL_NONFINITE_STATISTIC",
      observed_rho=ob,
      null_median_rho=NA_real_,
      null_q025=NA_real_,
      null_q975=NA_real_,
      empirical_p_greater=NA_real_,
      empirical_p_two_sided_abs=NA_real_,
      permutations=PERM_ITER
    )
    next
  }

  threshold_iterations[[as.character(thr)]] <- data.table(
    virus_family=fam,
    min_cluster_member_n=thr,
    iteration=seq_len(PERM_ITER),
    perm_rho=pr
  )

  threshold_summary[[as.character(thr)]] <- data.table(
    virus_family=fam,
    min_cluster_member_n=thr,
    tip_n=length(keep),
    total_possible_pair_n=nrow(apt),
    tested_pair_n=cpt$tested_n,
    sampling_mode=cpt$mode,
    status="PASS",
    observed_rho=ob,
    null_median_rho=as.numeric(median(pr[fin])),
    null_q025=as.numeric(safe_quant(pr[fin], .025)),
    null_q975=as.numeric(safe_quant(pr[fin], .975)),
    empirical_p_greater=(sum(pr[fin] >= ob) + 1)/(sum(fin) + 1),
    empirical_p_two_sided_abs=(sum(abs(pr[fin]) >= abs(ob)) + 1)/(sum(fin) + 1),
    permutations=PERM_ITER
  )
}

fwrite(
  rbindlist(threshold_summary, fill=TRUE),
  file.path(outdir, "07_min_support_sensitivity_summary.tsv"),
  sep="\t"
)

if (length(threshold_iterations)) {
  fwrite(
    rbindlist(threshold_iterations, fill=TRUE),
    file.path(outdir, "08_min_support_sensitivity_iterations.tsv"),
    sep="\t"
  )
}

# ---------------------------
# Alignment-filtering bias audit
# ---------------------------
mem <- fread(memfile, sep="\t")
need_mem <- c("representative_sequence_id","member_sequence_id","cluster_member_n")
miss_mem <- setdiff(need_mem, names(mem))
if (length(miss_mem)) {
  stop("Missing AA95 membership columns for ", fam, ": ", paste(miss_mem, collapse=","))
}

mem[, representative_sequence_id := as.character(representative_sequence_id)]
mem[, member_sequence_id := as.character(member_sequence_id)]
mem[, cluster_member_n := as.integer(cluster_member_n)]

assoc <- fread(
  ASSOC_FINAL,
  sep="\t",
  select=c("sequence_id","host_family"),
  showProgress=FALSE
)
assoc[, sequence_id := as.character(sequence_id)]
assoc[, host_family := as.character(host_family)]
assoc <- assoc[
  !is.na(sequence_id) & nzchar(sequence_id) &
  !is.na(host_family) & nzchar(host_family),
  .(sequence_id, host_family)
]
assoc <- unique(assoc)
setnames(assoc, "sequence_id", "member_sequence_id")

mh <- merge(
  mem[, .(representative_sequence_id, member_sequence_id, cluster_member_n)],
  assoc,
  by="member_sequence_id",
  all.x=TRUE,
  sort=FALSE
)

fullrep <- mh[, .(
  cluster_member_n=as.numeric(max(cluster_member_n, na.rm=TRUE)),
  host_family_n=as.numeric(uniqueN(host_family[!is.na(host_family) & nzchar(host_family)])),
  host_annotated=as.logical(any(!is.na(host_family) & nzchar(host_family)))
), by=representative_sequence_id]

if (fam == "Geminiviridae") {
  qcf <- file.path(
    kroot,
    "07_Geminiviridae_Rep_conserved",
    "AA95_alignment_QC",
    "Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv"
  )
} else {
  qcf <- file.path(
    kroot,
    "06_AA95_alignment_QC",
    fam,
    paste0(fam, ".alignment_sequence_QC.tsv")
  )
}
if (!file.exists(qcf)) stop("Missing alignment QC report: ", qcf)

aq <- fread(qcf, sep="\t")
if (!"tip" %in% names(aq)) stop("Alignment QC file lacks 'tip': ", qcf)
setnames(aq, "tip", "representative_sequence_id")
aq[, representative_sequence_id := as.character(representative_sequence_id)]

bias <- merge(
  fullrep,
  aq,
  by="representative_sequence_id",
  all.y=TRUE,
  sort=FALSE
)

bias[, drop_gap := as.logical(gap_fraction > 0.50)]
bias[, drop_effective := as.logical(effective_aa < 100)]
bias[, drop_reason := fifelse(
  action == "keep", "keep",
  fifelse(
    drop_gap & drop_effective, "gap_and_effective",
    fifelse(
      drop_gap, "gap",
      fifelse(drop_effective, "effective_aa", "other")
    )
  )
)]

fwrite(
  bias,
  file.path(outdir, "09_alignment_filtering_bias_detail.tsv"),
  sep="\t"
)

# Explicit coercion to numeric prevents integer/double inconsistency across groups.
bias_groups <- unique(bias[, .(action, drop_reason)])
bias_summary_list <- vector("list", nrow(bias_groups))

for (ii in seq_len(nrow(bias_groups))) {
  a <- bias_groups$action[ii]
  r <- bias_groups$drop_reason[ii]
  z <- bias[action == a & drop_reason == r]

  bias_summary_list[[ii]] <- data.table(
    virus_family=fam,
    action=as.character(a),
    drop_reason=as.character(r),
    representative_n=as.numeric(nrow(z)),
    host_annotated_n=as.numeric(sum(z$host_annotated %in% TRUE, na.rm=TRUE)),
    host_annotated_fraction=as.numeric(mean(z$host_annotated %in% TRUE, na.rm=TRUE)),
    cluster_member_median=as.numeric(median(z$cluster_member_n, na.rm=TRUE)),
    host_family_n_median=as.numeric(median(z$host_family_n, na.rm=TRUE)),
    effective_aa_median=as.numeric(median(z$effective_aa, na.rm=TRUE)),
    gap_fraction_median=as.numeric(median(z$gap_fraction, na.rm=TRUE))
  )
}

bias_summary <- rbindlist(bias_summary_list, fill=TRUE)
fwrite(
  bias_summary,
  file.path(outdir, "10_alignment_filtering_bias_summary.tsv"),
  sep="\t"
)

# ---------------------------
# Family QC
# ---------------------------
family_qc <- data.table(
  virus_family=fam,
  final_tree_tip_n=as.numeric(length(tr$tip.label)),
  host_annotated_tip_n=as.numeric(length(ordered_tips)),
  host_family_universe_n=as.numeric(length(host_universe)),
  distinct_host_spectrum_n=as.numeric(uniqueN(tip2$host_families)),
  total_possible_pair_n=as.numeric(total_pairs),
  tested_pair_n=as.numeric(chosen$tested_n),
  permutations=as.numeric(PERM_ITER),
  pair_sampling_replicates=as.numeric(if(total_pairs > MAX_PAIRS) PAIR_STABILITY_REPS else 1L),
  support_stratum_n=as.numeric(nlevels(strata)),
  min_support_thresholds="2;3",
  status="PASS"
)

fwrite(
  family_qc,
  file.path(outdir, "11_family_QC.tsv"),
  sep="\t"
)

cat("03L FINAL family PASS:", fam, "\n")
print(main_summary)
print(stab_summary)
print(support_summary)
print(rbindlist(threshold_summary, fill=TRUE))
print(bias_summary)
print(family_qc)
