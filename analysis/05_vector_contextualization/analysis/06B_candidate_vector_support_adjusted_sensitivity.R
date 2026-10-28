#!/usr/bin/env Rscript

# Fig. 4B - sampling/support-adjusted contextualization of the prioritized candidates
#
# Scientific question:
#   Do Fig.5 network-prioritized genera show more systematically documented VIP-DB vector
#   relationships than background genera after accounting for unequal public sequence support?
#
# This is a sensitivity/contextualization analysis, NOT a causal model and NOT a new
# candidate score. Host-family breadth and bridge score are deliberately NOT adjusted
# for because they define candidate status. The primary observation-support covariate
# is per-genus VC80 support (number of VC80 units with released sequence support);
# sequence_n is an alternative sensitivity covariate.
#
# Formal inputs:
#   Fig6/tables/01_network_genus_VIPDB_systematic_context.tsv
#   Fig5/tables/05_candidate_thresholds.tsv
#   Fig4/tables/01_VIPDB_relation_normalized.tsv
#
# Formal outputs:
#   Fig6/tables/11_candidate_VIPDB_support_eligible_enrichment.tsv
#   Fig6/tables/12_candidate_VIPDB_support_adjusted_models.tsv
#   Fig6/tables/13_network_genus_sampling_support.tsv
#   Fig6/tables/14_candidate_VIPDB_vector_ecology.tsv
#   Fig6/tables/15_candidate_VIPDB_adjusted_predictions.tsv
#   Fig6/qc/02_Fig6_support_adjusted_QC.tsv

suppressPackageStartupMessages({
  library(data.table)
})

default_root <- Sys.getenv("ANALYSIS_ROOT", unset = ".")
ANALYSIS_ROOT <- Sys.getenv("ANALYSIS_ROOT", unset = default_root)
ANALYSIS_ROOT <- normalizePath(ANALYSIS_ROOT, winslash = "/", mustWork = TRUE)

FIG6_BASE <- file.path(ANALYSIS_ROOT, "05_vector_contextualization")
FIG5_BASE <- file.path(ANALYSIS_ROOT, "04_genus_host_network")
FIG4_BASE <- file.path(ANALYSIS_ROOT, "03_arthropod_transmission_ecology")

TAB6 <- file.path(FIG6_BASE, "tables")
QC6  <- file.path(FIG6_BASE, "qc")
dir.create(TAB6, recursive = TRUE, showWarnings = FALSE)
dir.create(QC6, recursive = TRUE, showWarnings = FALSE)

P_CONTEXT <- file.path(TAB6, "01_network_genus_VIPDB_systematic_context.tsv")
P_THRESH  <- file.path(FIG5_BASE, "tables", "05_candidate_thresholds.tsv")
P_VIPREL  <- file.path(FIG4_BASE, "tables", "01_VIPDB_relation_normalized.tsv")

for (p in c(P_CONTEXT, P_THRESH, P_VIPREL)) {
  if (!file.exists(p)) stop("Required formal input not found: ", p)
}

ctx <- fread(P_CONTEXT, sep = "\t", showProgress = FALSE)
thr <- fread(P_THRESH,  sep = "\t", showProgress = FALSE)
vip <- fread(P_VIPREL, sep = "\t", showProgress = FALSE)

required_ctx <- c(
  "virus_node", "candidate", "virus_family",
  "host_family_n", "bridge_score",
  "qc_pass_VC80_support", "sequence_n",
  "VIPDB_documented_vector", "VIPDB_matched_taxid_n"
)
miss <- setdiff(required_ctx, names(ctx))
if (length(miss)) stop("Missing columns in Fig6 systematic context: ", paste(miss, collapse = ", "))

if (!"support_cutoff" %in% names(thr) || nrow(thr) != 1L) {
  stop("Fig5 candidate threshold table must contain exactly one support_cutoff row.")
}
support_cutoff <- as.numeric(thr$support_cutoff[1])
if (!is.finite(support_cutoff)) stop("Invalid support_cutoff in Fig5 threshold table.")

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  z <- toupper(trimws(as.character(x)))
  z %in% c("TRUE", "T", "1", "YES", "Y")
}
clean_text <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z)] <- ""
  z
}
collapse_unique <- function(x) {
  z <- sort(unique(clean_text(x)))
  z <- z[nzchar(z)]
  if (!length(z)) "" else paste(z, collapse = ";")
}

ctx[, candidate := as_flag(candidate)]
ctx[, VIPDB_documented_vector := as_flag(VIPDB_documented_vector)]

if (!is.logical(ctx$candidate)) stop("Internal type error: candidate must be logical.")
if (!is.logical(ctx$VIPDB_documented_vector)) stop("Internal type error: VIPDB_documented_vector must be logical.")
ctx[, qc_pass_VC80_support := as.numeric(qc_pass_VC80_support)]
ctx[, sequence_n := as.numeric(sequence_n)]

if (any(!is.finite(ctx$qc_pass_VC80_support))) stop("Non-finite qc_pass_VC80_support values.")
if (any(!is.finite(ctx$sequence_n))) stop("Non-finite sequence_n values.")

ctx[, `:=`(
  support_cutoff = support_cutoff,
  support_eligible = qc_pass_VC80_support >= support_cutoff,
  log1p_VC80_support = log1p(qc_pass_VC80_support),
  log1p_sequence_n = log1p(sequence_n),
  vector_binary = as.integer(VIPDB_documented_vector)
)]

# ------------------------------------------------------------------
# 11. Exact enrichment: all genera and support-eligible background
# ------------------------------------------------------------------
fisher_row <- function(d, analysis_set) {
  tab <- table(
    candidate = factor(d$candidate, levels = c(FALSE, TRUE)),
    vector = factor(d$VIPDB_documented_vector, levels = c(FALSE, TRUE))
  )
  ft <- fisher.test(tab)
  data.table(
    analysis_set = analysis_set,
    support_cutoff = support_cutoff,
    genus_n = nrow(d),
    candidate_n = sum(d$candidate),
    noncandidate_n = sum(!d$candidate),
    candidate_vector_n = sum(d$candidate & d$VIPDB_documented_vector),
    candidate_no_vector_n = sum(d$candidate & !d$VIPDB_documented_vector),
    noncandidate_vector_n = sum(!d$candidate & d$VIPDB_documented_vector),
    noncandidate_no_vector_n = sum(!d$candidate & !d$VIPDB_documented_vector),
    fisher_odds_ratio = unname(ft$estimate),
    fisher_ci_low = unname(ft$conf.int[1]),
    fisher_ci_high = unname(ft$conf.int[2]),
    fisher_p = ft$p.value
  )
}

enrich <- rbind(
  fisher_row(ctx, "all_network_genera"),
  fisher_row(ctx[support_eligible == TRUE], "support_eligible_genera")
)
fwrite(enrich, file.path(TAB6, "11_candidate_VIPDB_support_eligible_enrichment.tsv"), sep = "\t")

# ------------------------------------------------------------------
# 12. Logistic sensitivity models
#     These estimate the candidate term after controlling observation
#     support; they are not causal and may be imprecise because candidate
#     genera occupy the high-support tail by design.
# ------------------------------------------------------------------
fit_candidate_model <- function(d, formula, model_id, analysis_set, support_variable) {
  fit <- glm(formula, data = d, family = binomial())
  sm <- summary(fit)$coefficients

  candidate_term <- grep("^candidateTRUE$", rownames(sm), value = TRUE)
  if (!length(candidate_term)) {
    stop("Candidate coefficient not estimable for model: ", model_id)
  }
  candidate_term <- candidate_term[1]
  beta <- sm[candidate_term, "Estimate"]
  se   <- sm[candidate_term, "Std. Error"]
  pval <- sm[candidate_term, "Pr(>|z|)"]

  support_term <- setdiff(rownames(sm), c("(Intercept)", candidate_term))
  support_beta <- support_se <- support_p <- NA_real_
  if (length(support_term)) {
    st <- support_term[1]
    support_beta <- sm[st, "Estimate"]
    support_se   <- sm[st, "Std. Error"]
    support_p    <- sm[st, "Pr(>|z|)"]
  }

  data.table(
    model_id = model_id,
    analysis_set = analysis_set,
    support_variable = support_variable,
    formula = paste(deparse(formula), collapse = ""),
    genus_n = nrow(d),
    candidate_n = sum(d$candidate),
    vector_positive_n = sum(d$VIPDB_documented_vector),
    candidate_beta = beta,
    candidate_se = se,
    candidate_OR = exp(beta),
    candidate_OR_ci_low = exp(beta - 1.96 * se),
    candidate_OR_ci_high = exp(beta + 1.96 * se),
    candidate_p = pval,
    support_beta = support_beta,
    support_se = support_se,
    support_p = support_p,
    AIC = AIC(fit),
    converged = isTRUE(fit$converged)
  )
}

# Keep candidate as LOGICAL throughout.
# glm() handles a logical binary predictor directly and names the coefficient
# candidateTRUE. This avoids accidental factor arithmetic in summary counts.
eligible <- ctx[support_eligible == TRUE]

models <- rbindlist(list(
  fit_candidate_model(
    ctx,
    vector_binary ~ candidate,
    "M0_unadjusted_all",
    "all_network_genera",
    "none"
  ),
  fit_candidate_model(
    ctx,
    vector_binary ~ candidate + log1p_VC80_support,
    "M1_VC80_adjusted_all",
    "all_network_genera",
    "log1p(qc_pass_VC80_support)"
  ),
  fit_candidate_model(
    eligible,
    vector_binary ~ candidate + log1p_VC80_support,
    "M2_VC80_adjusted_support_eligible",
    "support_eligible_genera",
    "log1p(qc_pass_VC80_support)"
  ),
  fit_candidate_model(
    ctx,
    vector_binary ~ candidate + log1p_sequence_n,
    "M3_sequence_adjusted_all",
    "all_network_genera",
    "log1p(sequence_n)"
  ),
  fit_candidate_model(
    eligible,
    vector_binary ~ candidate + log1p_sequence_n,
    "M4_sequence_adjusted_support_eligible",
    "support_eligible_genera",
    "log1p(sequence_n)"
  )
))
fwrite(models, file.path(TAB6, "12_candidate_VIPDB_support_adjusted_models.tsv"), sep = "\t")

# ------------------------------------------------------------------
# 13. Genus-level support audit
# ------------------------------------------------------------------
support_audit <- ctx[, .(
  virus_node,
  virus_family,
  candidate = as.logical(candidate),
  host_family_n,
  bridge_score,
  qc_pass_VC80_support,
  sequence_n,
  log1p_VC80_support,
  log1p_sequence_n,
  support_cutoff,
  support_eligible,
  VIPDB_documented_vector,
  VIPDB_matched_taxid_n
)]
fwrite(support_audit, file.path(TAB6, "13_network_genus_sampling_support.tsv"), sep = "\t")

# ------------------------------------------------------------------
# 14. Candidate vector ecology from systematic normalized VIPDB table
# ------------------------------------------------------------------
cand_meta <- support_audit[candidate == TRUE]

vip[, virus_genus := clean_text(virus_genus)]
vip_cand <- vip[virus_genus %in% cand_meta$virus_node]

if (nrow(vip_cand)) {
  eco <- vip_cand[, .(
    VIPDB_relation_n = .N,
    VIPDB_virus_taxid_n = uniqueN(virus_taxid[!is.na(virus_taxid)]),
    vector_order_n = uniqueN(vector_order[nzchar(clean_text(vector_order))]),
    vector_family_n = uniqueN(vector_family[nzchar(clean_text(vector_family))]),
    vector_genus_n = uniqueN(vector_genus[nzchar(clean_text(vector_genus))]),
    vector_name_n = uniqueN(vector_name[nzchar(clean_text(vector_name))]),
    transmission_mode_n = uniqueN(fine4_mode[nzchar(clean_text(fine4_mode))]),
    fine4_modes = collapse_unique(fine4_mode),
    retention_groups = collapse_unique(retention_group),
    vector_orders = collapse_unique(vector_order),
    vector_families = collapse_unique(vector_family)
  ), by = .(virus_node = virus_genus)]
} else {
  eco <- data.table(virus_node = character())
}

candidate_ecology <- merge(cand_meta, eco, by = "virus_node", all.x = TRUE)
count_cols <- c(
  "VIPDB_relation_n", "VIPDB_virus_taxid_n", "vector_order_n", "vector_family_n",
  "vector_genus_n", "vector_name_n", "transmission_mode_n"
)
for (z in count_cols) {
  if (!z %in% names(candidate_ecology)) candidate_ecology[, (z) := 0L]
  candidate_ecology[is.na(get(z)), (z) := 0L]
}
text_cols <- c("fine4_modes", "retention_groups", "vector_orders", "vector_families")
for (z in text_cols) {
  if (!z %in% names(candidate_ecology)) candidate_ecology[, (z) := ""]
  candidate_ecology[is.na(get(z)), (z) := ""]
}
setorder(candidate_ecology, -host_family_n, -bridge_score, virus_node)
fwrite(candidate_ecology, file.path(TAB6, "14_candidate_VIPDB_vector_ecology.tsv"), sep = "\t")

# ------------------------------------------------------------------
# 15. Per-genus fitted probabilities for Fig. 4 sensitivity visualization
# ------------------------------------------------------------------
fit_vc_all <- glm(vector_binary ~ candidate + log1p_VC80_support, data = ctx, family = binomial())
fit_seq_all <- glm(vector_binary ~ candidate + log1p_sequence_n, data = ctx, family = binomial())
fit_vc_eligible <- glm(vector_binary ~ candidate + log1p_VC80_support, data = eligible, family = binomial())

pred <- support_audit[, .(
  virus_node, virus_family, candidate,
  qc_pass_VC80_support, sequence_n, support_eligible,
  VIPDB_documented_vector
)]
pred[, predicted_vector_prob_VC80_adjusted_all :=
       as.numeric(predict(fit_vc_all, newdata = ctx, type = "response"))]
pred[, predicted_vector_prob_sequence_adjusted_all :=
       as.numeric(predict(fit_seq_all, newdata = ctx, type = "response"))]

pred[, predicted_vector_prob_VC80_adjusted_eligible := NA_real_]
idx <- which(pred$support_eligible)
pred$predicted_vector_prob_VC80_adjusted_eligible[idx] <-
  as.numeric(predict(fit_vc_eligible, newdata = eligible, type = "response"))

fwrite(pred, file.path(TAB6, "15_candidate_VIPDB_adjusted_predictions.tsv"), sep = "\t")

# ------------------------------------------------------------------
# QC / run summary
# ------------------------------------------------------------------
qc <- data.table(
  metric = c(
    "status",
    "network_genus_n",
    "candidate_n",
    "candidate_vector_positive_n",
    "support_cutoff",
    "support_eligible_genus_n",
    "support_eligible_noncandidate_n",
    "support_eligible_candidate_n",
    "all_genera_fisher_p",
    "support_eligible_fisher_p",
    "VC80_adjusted_all_candidate_OR",
    "VC80_adjusted_all_candidate_p",
    "VC80_adjusted_eligible_candidate_OR",
    "VC80_adjusted_eligible_candidate_p",
    "sequence_adjusted_all_candidate_OR",
    "sequence_adjusted_all_candidate_p"
  ),
  value = as.character(c(
    "PASS",
    nrow(ctx),
    sum(ctx$candidate),
    sum(ctx$candidate & ctx$VIPDB_documented_vector),
    support_cutoff,
    nrow(eligible),
    sum(!eligible$candidate),
    sum(eligible$candidate),
    enrich[analysis_set == "all_network_genera", fisher_p],
    enrich[analysis_set == "support_eligible_genera", fisher_p],
    models[model_id == "M1_VC80_adjusted_all", candidate_OR],
    models[model_id == "M1_VC80_adjusted_all", candidate_p],
    models[model_id == "M2_VC80_adjusted_support_eligible", candidate_OR],
    models[model_id == "M2_VC80_adjusted_support_eligible", candidate_p],
    models[model_id == "M3_sequence_adjusted_all", candidate_OR],
    models[model_id == "M3_sequence_adjusted_all", candidate_p]
  ))
)
fwrite(qc, file.path(QC6, "02_Fig6_support_adjusted_QC.tsv"), sep = "\t")

cat("Fig6 06B support-adjusted vector contextualization PASS\n")
cat("ANALYSIS_ROOT:", ANALYSIS_ROOT, "\n")
cat("Support cutoff:", support_cutoff, "\n\n")
cat("Exact enrichment:\n")
print(enrich)
cat("\nCandidate-effect sensitivity models:\n")
print(models[, .(
  model_id, analysis_set, support_variable,
  candidate_OR, candidate_OR_ci_low, candidate_OR_ci_high, candidate_p,
  converged
)])
cat("\nInterpretation guardrail:\n")
cat(
  "These are sampling/support sensitivities, not causal estimates. ",
  "Attenuation after support adjustment means the raw candidate-vector enrichment ",
  "is partly associated with public evidence intensity and should not be presented ",
  "as independent biological validation of candidate status.\n",
  sep = ""
)
