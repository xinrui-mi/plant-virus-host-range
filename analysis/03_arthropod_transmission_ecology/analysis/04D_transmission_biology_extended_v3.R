#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# 04D_transmission_biology_extended_v3.R
#
# Purpose
# -------
# Extend, but do not replace, the frozen 04A-04C analyses.
# v3 separates global lineage-composition sensitivity from within-mode lineage influence.
#
# The analysis keeps S0-S2 strictly as arthropod-associated sequence evidence.
# VIP-DB is used for documented vector relationships and transmission mode.
#
# New analyses:
#  1) Five-group comparison:
#       no documented vector relationship in VIP-DB
#       non-persistent
#       semi-persistent
#       persistent non-propagative
#       persistent propagative
#     Only viruses with exactly one documented four-class transmission mode are
#     assigned to a documented-mode group. Viruses with ambiguous/multiple
#     modes are excluded from five-group inferential models and audited.
#  2) Host-family and host-order breadth:
#       quasi-Poisson, crude -> +VC80 support -> +virus family
#  3) Maximum host taxonomic distance:
#       linear model, crude -> +VC80 support -> +virus family
#  4) Lineage dependence:
#       transmission-mode composition by family/genus
#       systematic leave-one-family-out and leave-one-genus-out
#       within-mode lineage deletion, preserving the reference-group composition
#  5) Geographic breadth:
#       country_n, with association-record adjustment as primary sampling
#       adjustment; VC80 support as sensitivity.
#  6) S0-S2 validation:
#       association with VIP-DB documented-vector status
#       association with country breadth after record/family adjustment
#
# This script makes no causal claim and does not infer vector competence from
# sequence evidence.
# =============================================================================


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Cannot determine script path from --file.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
script_dir  <- dirname(script_path)
base        <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(base, ".."), mustWork = TRUE)

# Released layout: 04D needs only five columns from each master
# (sequence_id, virus_taxid, virus_family, virus_genus, arthropod_sequence_supported
# from the sequence master; sequence_id plus host fields from the association
# master), all of which are present in the released analysis inputs. Those are
# therefore the defaults. On the original cluster the full masters were used; set
# SEQ_FINAL / ASSOC_FINAL to use them instead.
SEQ_FINAL_PATH <- Sys.getenv(
  "SEQ_FINAL",
  unset = file.path(dirname(analysis_root), "data", "analysis_input",
                    "sequence_analysis_input.tsv.gz")
)
ASSOC_FINAL_PATH <- Sys.getenv(
  "ASSOC_FINAL",
  unset = file.path(dirname(analysis_root), "data", "analysis_input",
                    "association_analysis_input.tsv.gz")
)

tab_dir <- file.path(base, "tables")
qc_dir  <- file.path(base, "qc")
log_dir <- file.path(base, "logs")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

metric_file <- file.path(tab_dir, "04_formal_taxid_arthropod_vector_metrics.tsv")
vip_relation_file <- file.path(tab_dir, "01_VIPDB_relation_normalized.tsv")

seq_final <- SEQ_FINAL_PATH
assoc_final <- ASSOC_FINAL_PATH

for (p in c(metric_file, vip_relation_file, seq_final, assoc_final)) {
  if (!file.exists(p) || file.info(p)$size <= 0) {
    stop("Missing or empty required input: ", p)
  }
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

present <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) & nzchar(y) &
    !(tolower(y) %in% c("na", "n/a", "nan", "none", "null", "."))
}

std_text <- function(x) {
  y <- trimws(as.character(x))
  y[!present(y)] <- NA_character_
  y
}

as_bool <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

num <- function(x) suppressWarnings(as.numeric(x))

safe_unique_one <- function(x) {
  z <- unique(std_text(x))
  z <- z[present(z)]
  if (length(z) == 1L) z else NA_character_
}

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

canonical_mode <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- gsub("[_/-]+", " ", s)
  s <- gsub("\\s+", " ", s)

  out <- rep(NA_character_, length(s))

  # Order matters.
  out[
    grepl("persistent", s) &
      grepl("propagative|propagative", s) &
      !grepl("non propagative|nonpropagative", s)
  ] <- "Persistent propagative"

  out[
    grepl("persistent", s) &
      (grepl("non propagative|nonpropagative|circulative", s))
  ] <- "Persistent non-propagative"

  out[grepl("semi persistent|semipersistent", s)] <- "Semi-persistent"

  out[
    grepl("non persistent|nonpersistent", s) &
      !grepl("semi", s)
  ] <- "Non-persistent"

  out
}

MODE_LEVELS <- c(
  "No documented vector relationship",
  "Non-persistent",
  "Semi-persistent",
  "Persistent non-propagative",
  "Persistent propagative"
)

tidy_fit <- function(fit, model_id, outcome, scenario, model_type) {
  if (is.null(fit)) {
    return(data.table(
      model_id = model_id,
      outcome = outcome,
      scenario = scenario,
      model_type = model_type,
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      exponentiated_effect = NA_real_,
      exponentiated_ci_low = NA_real_,
      exponentiated_ci_high = NA_real_,
      n = NA_integer_,
      residual_df = NA_real_,
      converged = FALSE
    ))
  }

  sm <- summary(fit)$coefficients
  if (is.null(dim(sm))) {
    sm <- matrix(sm, nrow = 1L)
    rownames(sm) <- names(coef(fit))[1L]
  }

  out <- data.table(
    term = rownames(sm),
    estimate = as.numeric(sm[, 1]),
    std_error = as.numeric(sm[, 2]),
    statistic = as.numeric(sm[, 3]),
    p_value = as.numeric(sm[, 4])
  )

  out[, `:=`(
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error
  )]

  # Exponentiation is directly interpretable for log-link count models and
  # logistic models. It is retained as a convenience column for all models,
  # but effect_scale below defines the valid interpretation.
  out[, `:=`(
    exponentiated_effect = exp(estimate),
    exponentiated_ci_low = exp(ci_low),
    exponentiated_ci_high = exp(ci_high)
  )]

  out[, `:=`(
    model_id = model_id,
    outcome = outcome,
    scenario = scenario,
    model_type = model_type,
    n = nobs(fit),
    residual_df = df.residual(fit),
    converged = if (!is.null(fit$converged)) isTRUE(fit$converged) else TRUE
  )]

  setcolorder(
    out,
    c(
      "model_id", "outcome", "scenario", "model_type", "term",
      "estimate", "std_error", "statistic", "p_value",
      "ci_low", "ci_high",
      "exponentiated_effect", "exponentiated_ci_low",
      "exponentiated_ci_high",
      "n", "residual_df", "converged"
    )
  )
  out
}

fit_qp <- function(formula, d, model_id, outcome, scenario) {
  fit <- tryCatch(
    glm(
      formula,
      data = d,
      family = quasipoisson(link = "log"),
      control = glm.control(maxit = 200)
    ),
    error = function(e) NULL
  )
  tidy_fit(fit, model_id, outcome, scenario, "quasi-Poisson(log)")
}

fit_lm <- function(formula, d, model_id, outcome, scenario) {
  fit <- tryCatch(
    lm(formula, data = d),
    error = function(e) NULL
  )
  tidy_fit(fit, model_id, outcome, scenario, "linear")
}

fit_logistic <- function(formula, d, model_id, outcome, scenario) {
  fit <- tryCatch(
    glm(
      formula,
      data = d,
      family = binomial(),
      control = glm.control(maxit = 200)
    ),
    error = function(e) NULL
  )
  tidy_fit(fit, model_id, outcome, scenario, "binomial(logit)")
}

primary_terms <- function(x, exposure_prefix) {
  x[
    !is.na(term) &
      startsWith(term, exposure_prefix)
  ]
}

prep_model_data <- function(d) {
  z <- copy(d)
  z <- z[
    present(virus_family) &
      present(transmission_group)
  ]
  z[, transmission_group := factor(
    transmission_group,
    levels = MODE_LEVELS
  )]
  z <- z[!is.na(transmission_group)]
  z[, virus_family := factor(virus_family)]
  z[, virus_genus := factor(virus_genus)]
  z
}


# -----------------------------------------------------------------------------
# Load frozen tables
# -----------------------------------------------------------------------------

metric <- fread(metric_file, sep = "\t", showProgress = FALSE)
viprel <- fread(vip_relation_file, sep = "\t", showProgress = FALSE)
seq <- fread(seq_final, sep = "\t", showProgress = FALSE)
assoc <- fread(assoc_final, sep = "\t", showProgress = FALSE)

required_metric <- c(
  "virus_taxid",
  "host_family_n",
  "host_order_n",
  "max_host_taxonomic_distance",
  "qc_pass_VC80_support",
  "documented_vector",
  "fine4_status",
  "fine4_mode"
)
missing_metric <- setdiff(required_metric, names(metric))
if (length(missing_metric)) {
  stop(
    "04_formal_taxid_arthropod_vector_metrics.tsv missing required columns: ",
    paste(missing_metric, collapse = ", ")
  )
}

required_seq <- c(
  "sequence_id", "virus_taxid", "virus_family", "virus_genus",
  "arthropod_sequence_supported"
)
missing_seq <- setdiff(required_seq, names(seq))
if (length(missing_seq)) {
  stop(
    "Final sequence master missing required columns: ",
    paste(missing_seq, collapse = ", ")
  )
}

required_assoc <- c(
  "virus_taxid", "country"
)
missing_assoc <- setdiff(required_assoc, names(assoc))
if (length(missing_assoc)) {
  stop(
    "Final association master missing required columns: ",
    paste(missing_assoc, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# Taxonomy / evidence / geography additions from frozen Stage-20 masters
# -----------------------------------------------------------------------------

seq[, virus_taxid := std_text(virus_taxid)]
assoc[, virus_taxid := std_text(virus_taxid)]
metric[, virus_taxid := std_text(virus_taxid)]

tax_audit <- seq[present(virus_taxid), .(
  virus_family_final = safe_unique_one(virus_family),
  virus_genus_final = safe_unique_one(virus_genus),
  family_label_n = uniqueN(std_text(virus_family)[present(virus_family)]),
  genus_label_n = uniqueN(std_text(virus_genus)[present(virus_genus)]),
  sequence_n_stage20 = uniqueN(sequence_id),
  seq_supported = any(
    as_bool(arthropod_sequence_supported) %in% TRUE,
    na.rm = TRUE
  )
), by = virus_taxid]

geo <- assoc[present(virus_taxid), .(
  association_n_stage20 = .N,
  country_n_stage20 = uniqueN(std_text(country)[present(country)])
), by = virus_taxid]

metric <- merge(metric, tax_audit, by = "virus_taxid", all.x = TRUE)
metric <- merge(metric, geo, by = "virus_taxid", all.x = TRUE)

# Prefer taxonomy already present in the frozen Fig4 metric table when valid;
# otherwise use the final Stage-20 taxonomy derived above.
if (!"virus_family" %in% names(metric)) {
  metric[, virus_family := virus_family_final]
} else {
  metric[, virus_family := fifelse(
    present(virus_family),
    std_text(virus_family),
    virus_family_final
  )]
}

if (!"virus_genus" %in% names(metric)) {
  metric[, virus_genus := virus_genus_final]
} else {
  metric[, virus_genus := fifelse(
    present(virus_genus),
    std_text(virus_genus),
    virus_genus_final
  )]
}

# Use existing metric-table counts when present; otherwise Stage-20-derived ones.
if (!"association_n" %in% names(metric)) {
  metric[, association_n := association_n_stage20]
} else {
  metric[!is.finite(num(association_n)), association_n := association_n_stage20]
}

if (!"country_n" %in% names(metric)) {
  metric[, country_n := country_n_stage20]
} else {
  metric[!is.finite(num(country_n)), country_n := country_n_stage20]
}

if (!"sequence_n" %in% names(metric)) {
  metric[, sequence_n := sequence_n_stage20]
}

metric[, `:=`(
  host_family_n = num(host_family_n),
  host_order_n = num(host_order_n),
  max_host_taxonomic_distance = num(max_host_taxonomic_distance),
  qc_pass_VC80_support = num(qc_pass_VC80_support),
  association_n = num(association_n),
  country_n = num(country_n),
  sequence_n = num(sequence_n),
  documented_vector = as_bool(documented_vector),
  fine4_status = std_text(fine4_status),
  fine4_mode = std_text(fine4_mode),
  virus_family = std_text(virus_family),
  virus_genus = std_text(virus_genus)
)]

metric[is.na(seq_supported), seq_supported := FALSE]


# -----------------------------------------------------------------------------
# Build the five-group exposure
# -----------------------------------------------------------------------------

metric[, canonical_four_mode := canonical_mode(fine4_mode)]

metric[, transmission_group := NA_character_]

metric[
  documented_vector %in% FALSE,
  transmission_group := "No documented vector relationship"
]

metric[
  documented_vector %in% TRUE &
    fine4_status == "single_four_class_mode" &
    present(canonical_four_mode),
  transmission_group := canonical_four_mode
]

metric[, transmission_group := factor(
  transmission_group,
  levels = MODE_LEVELS
)]

metric[, five_group_included := !is.na(transmission_group)]

exclusion_audit <- metric[, .(
  virus_taxid_n = .N
), by = .(
  documented_vector,
  fine4_status,
  fine4_mode,
  canonical_four_mode,
  five_group_included
)]
setorder(exclusion_audit, -virus_taxid_n)

write_tsv(
  exclusion_audit,
  file.path(tab_dir, "19_04D_five_group_inclusion_audit.tsv")
)


# -----------------------------------------------------------------------------
# Descriptive five-group summary
# -----------------------------------------------------------------------------

five <- prep_model_data(metric[five_group_included == TRUE])

five_summary <- five[, .(
  virus_taxid_n = .N,
  virus_family_n = uniqueN(virus_family),
  virus_genus_n = uniqueN(virus_genus[present(virus_genus)]),

  host_family_median = median(host_family_n, na.rm = TRUE),
  host_family_IQR = IQR(host_family_n, na.rm = TRUE),

  host_order_median = median(host_order_n, na.rm = TRUE),
  host_order_IQR = IQR(host_order_n, na.rm = TRUE),

  max_distance_median = median(max_host_taxonomic_distance, na.rm = TRUE),
  max_distance_IQR = IQR(max_host_taxonomic_distance, na.rm = TRUE),

  country_median = median(country_n, na.rm = TRUE),
  country_IQR = IQR(country_n, na.rm = TRUE),

  association_n_median = median(association_n, na.rm = TRUE),
  VC80_support_median = median(qc_pass_VC80_support, na.rm = TRUE)
), by = transmission_group]

write_tsv(
  five_summary,
  file.path(tab_dir, "20_04D_five_group_descriptive_summary.tsv")
)


# -----------------------------------------------------------------------------
# Composition audits
# -----------------------------------------------------------------------------

family_comp <- five[, .(
  virus_taxid_n = .N
), by = .(
  transmission_group,
  virus_family
)]
family_comp[, group_fraction := virus_taxid_n / sum(virus_taxid_n), by = transmission_group]
setorder(family_comp, transmission_group, -virus_taxid_n)

write_tsv(
  family_comp,
  file.path(tab_dir, "21_04D_transmission_group_by_virus_family.tsv")
)

genus_comp <- five[present(virus_genus), .(
  virus_taxid_n = .N
), by = .(
  transmission_group,
  virus_family,
  virus_genus
)]
genus_comp[, group_fraction := virus_taxid_n / sum(virus_taxid_n), by = transmission_group]
setorder(genus_comp, transmission_group, -virus_taxid_n)

write_tsv(
  genus_comp,
  file.path(tab_dir, "22_04D_transmission_group_by_virus_genus.tsv")
)

family_overlap <- five[, .(
  virus_taxid_n = .N,
  transmission_group_n = uniqueN(transmission_group),
  transmission_groups = paste(
    levels(droplevels(transmission_group)),
    collapse = " | "
  )
), by = virus_family]
setorder(family_overlap, -transmission_group_n, -virus_taxid_n)

write_tsv(
  family_overlap,
  file.path(tab_dir, "23_04D_within_family_mode_overlap_audit.tsv")
)


# Vector-family x virus-genus composition from normalized relation rows.
required_viprel <- c(
  "virus_taxid", "virus_family", "virus_genus",
  "vector_family", "fine4_mode"
)
missing_viprel <- setdiff(required_viprel, names(viprel))

if (!length(missing_viprel)) {
  viprel[, `:=`(
    virus_taxid = std_text(virus_taxid),
    virus_family = std_text(virus_family),
    virus_genus = std_text(virus_genus),
    vector_family = std_text(vector_family),
    canonical_four_mode = canonical_mode(fine4_mode)
  )]

  vector_system <- unique(
    viprel[
      present(virus_taxid) &
        present(virus_genus) &
        present(vector_family),
      .(
        virus_taxid,
        virus_family,
        virus_genus,
        vector_family,
        canonical_four_mode
      )
    ]
  )[, .(
    virus_taxid_n = uniqueN(virus_taxid)
  ), by = .(
    vector_family,
    virus_family,
    virus_genus,
    canonical_four_mode
  )]

  setorder(vector_system, -virus_taxid_n)

} else {
  vector_system <- data.table(
    vector_family = character(),
    virus_family = character(),
    virus_genus = character(),
    canonical_four_mode = character(),
    virus_taxid_n = integer()
  )
}

write_tsv(
  vector_system,
  file.path(tab_dir, "24_04D_vector_family_by_virus_genus.tsv")
)


# -----------------------------------------------------------------------------
# Main five-group host-range models
# -----------------------------------------------------------------------------

all_model_terms <- list()
k <- 0L

for (outcome in c("host_family_n", "host_order_n")) {
  d <- five[
    is.finite(get(outcome)) &
      get(outcome) >= 0 &
      is.finite(qc_pass_VC80_support)
  ]

  if (nrow(d) > 0 && uniqueN(d$transmission_group) >= 2) {
    k <- k + 1L
    all_model_terms[[k]] <- fit_qp(
      as.formula(paste0(outcome, " ~ transmission_group")),
      d,
      paste0(outcome, "_M0_crude"),
      outcome,
      "all_data"
    )

    k <- k + 1L
    all_model_terms[[k]] <- fit_qp(
      as.formula(
        paste0(
          outcome,
          " ~ transmission_group + log1p(qc_pass_VC80_support)"
        )
      ),
      d,
      paste0(outcome, "_M1_VC80_adjusted"),
      outcome,
      "all_data"
    )

    d2 <- d[present(virus_family)]
    if (nrow(d2) > 0 && uniqueN(d2$transmission_group) >= 2) {
      k <- k + 1L
      all_model_terms[[k]] <- fit_qp(
        as.formula(
          paste0(
            outcome,
            " ~ transmission_group + log1p(qc_pass_VC80_support)",
            " + factor(virus_family)"
          )
        ),
        d2,
        paste0(outcome, "_M2_VC80_family_adjusted"),
        outcome,
        "all_data"
      )
    }
  }
}

# Maximum host-taxonomic distance: preserve the linear-model framework used
# in the frozen 04B adjusted sensitivity.
d <- five[
  is.finite(max_host_taxonomic_distance) &
    is.finite(qc_pass_VC80_support)
]

if (nrow(d) > 0 && uniqueN(d$transmission_group) >= 2) {
  k <- k + 1L
  all_model_terms[[k]] <- fit_lm(
    max_host_taxonomic_distance ~ transmission_group,
    d,
    "max_distance_M0_crude",
    "max_host_taxonomic_distance",
    "all_data"
  )

  k <- k + 1L
  all_model_terms[[k]] <- fit_lm(
    max_host_taxonomic_distance ~
      transmission_group +
      log1p(qc_pass_VC80_support),
    d,
    "max_distance_M1_VC80_adjusted",
    "max_host_taxonomic_distance",
    "all_data"
  )

  d2 <- d[present(virus_family)]
  if (nrow(d2) > 0 && uniqueN(d2$transmission_group) >= 2) {
    k <- k + 1L
    all_model_terms[[k]] <- fit_lm(
      max_host_taxonomic_distance ~
        transmission_group +
        log1p(qc_pass_VC80_support) +
        factor(virus_family),
      d2,
      "max_distance_M2_VC80_family_adjusted",
      "max_host_taxonomic_distance",
      "all_data"
    )
  }
}

mode_model_terms <- rbindlist(all_model_terms, fill = TRUE)

write_tsv(
  mode_model_terms,
  file.path(tab_dir, "25_04D_five_group_host_range_model_terms.tsv")
)

mode_primary <- primary_terms(mode_model_terms, "transmission_group")
mode_primary[, effect_scale := fifelse(
  model_type == "quasi-Poisson(log)",
  "rate ratio vs no documented vector relationship",
  "difference vs no documented vector relationship"
)]

write_tsv(
  mode_primary,
  file.path(tab_dir, "26_04D_five_group_host_range_primary_contrasts.tsv")
)


# -----------------------------------------------------------------------------
# Systematic leave-one-lineage-out sensitivity
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# No lineage is pre-selected on the basis of the observed host-range result.
# Every observed virus family and every observed virus genus is evaluated in
# turn, provided that removing that lineage leaves all five transmission groups
# represented. This avoids outcome-driven selection of Begomovirus,
# Geminiviridae, or any other lineage.
#
# The same fully adjusted models used in the main analysis are refitted after
# each single-lineage removal. Driver strength is summarized by the change in
# the transmission-group coefficient relative to the full-data model.
# -----------------------------------------------------------------------------

full_primary_m2 <- mode_primary[
  grepl("_M2_VC80_family_adjusted$", model_id),
  .(
    outcome,
    term,
    full_estimate = estimate,
    full_std_error = std_error,
    full_p_value = p_value,
    full_ci_low = ci_low,
    full_ci_high = ci_high,
    full_model_type = model_type
  )
]

if (!nrow(full_primary_m2)) {
  stop("No fully adjusted transmission-group coefficients available for leave-one-lineage-out analysis.")
}

mode_full_counts <- five[, .(
  full_mode_n = .N
), by = .(
  transmission_group = as.character(transmission_group)
)]

loo_fit_one <- function(d, lineage_level, lineage_name, lineage_col) {
  keep <- is.na(d[[lineage_col]]) | as.character(d[[lineage_col]]) != lineage_name
  dd <- prep_model_data(d[keep])

  removed <- d[!keep]
  removed_total_n <- nrow(removed)

  after_counts <- dd[, .N, by = transmission_group]
  observed_after <- as.character(after_counts$transmission_group)
  all_groups_retained <- all(MODE_LEVELS %in% observed_after)

  audit_row <- data.table(
    lineage_level = lineage_level,
    removed_lineage = lineage_name,
    removed_total_n = removed_total_n,
    remaining_n = nrow(dd),
    all_five_groups_retained = all_groups_retained,
    eligible = all_groups_retained,
    exclusion_reason = if (all_groups_retained) "" else "removal_eliminated_at_least_one_transmission_group"
  )

  for (g in MODE_LEVELS) {
    audit_row[[paste0("removed_n__", g)]] <- removed[
      as.character(transmission_group) == g,
      .N
    ]
    audit_row[[paste0("remaining_n__", g)]] <- dd[
      as.character(transmission_group) == g,
      .N
    ]
  }

  if (!all_groups_retained) {
    return(list(audit = audit_row, terms = data.table()))
  }

  out <- list()
  k_local <- 0L

  for (outcome in c("host_family_n", "host_order_n")) {
    z <- dd[
      is.finite(get(outcome)) &
        is.finite(qc_pass_VC80_support) &
        present(virus_family)
    ]

    if (nrow(z) >= 30 && all(MODE_LEVELS %in% as.character(unique(z$transmission_group)))) {
      k_local <- k_local + 1L
      out[[k_local]] <- fit_qp(
        as.formula(
          paste0(
            outcome,
            " ~ transmission_group + log1p(qc_pass_VC80_support)",
            " + factor(virus_family)"
          )
        ),
        z,
        paste0(outcome, "_M2_VC80_family_adjusted"),
        outcome,
        paste0("remove_", lineage_level, ":", lineage_name)
      )
    }
  }

  z <- dd[
    is.finite(max_host_taxonomic_distance) &
      is.finite(qc_pass_VC80_support) &
      present(virus_family)
  ]

  if (nrow(z) >= 30 && all(MODE_LEVELS %in% as.character(unique(z$transmission_group)))) {
    k_local <- k_local + 1L
    out[[k_local]] <- fit_lm(
      max_host_taxonomic_distance ~
        transmission_group +
        log1p(qc_pass_VC80_support) +
        factor(virus_family),
      z,
      "max_distance_M2_VC80_family_adjusted",
      "max_host_taxonomic_distance",
      paste0("remove_", lineage_level, ":", lineage_name)
    )
  }

  terms <- if (length(out)) {
    primary_terms(rbindlist(out, fill = TRUE), "transmission_group")
  } else {
    data.table()
  }

  if (nrow(terms)) {
    terms[, `:=`(
      lineage_level = lineage_level,
      removed_lineage = lineage_name,
      removed_total_n = removed_total_n,
      remaining_n_after_removal = nrow(dd)
    )]
  }

  list(audit = audit_row, terms = terms)
}

lineage_specs <- list(
  virus_family = "virus_family",
  virus_genus = "virus_genus"
)

loo_audit_rows <- list()
loo_term_rows <- list()
la <- 0L
lt <- 0L

for (lineage_level in names(lineage_specs)) {
  lineage_col <- lineage_specs[[lineage_level]]
  lineage_values <- sort(unique(std_text(five[[lineage_col]])))
  lineage_values <- lineage_values[present(lineage_values)]

  for (lineage_name in lineage_values) {
    res <- loo_fit_one(
      five,
      lineage_level = lineage_level,
      lineage_name = lineage_name,
      lineage_col = lineage_col
    )

    la <- la + 1L
    loo_audit_rows[[la]] <- res$audit

    if (nrow(res$terms)) {
      lt <- lt + 1L
      loo_term_rows[[lt]] <- res$terms
    }
  }
}

loo_audit <- rbindlist(loo_audit_rows, fill = TRUE)
loo_terms <- rbindlist(loo_term_rows, fill = TRUE)

write_tsv(
  loo_audit,
  file.path(tab_dir, "35_04D_leave_one_lineage_out_eligibility_audit.tsv")
)

if (nrow(loo_terms)) {
  loo_primary <- merge(
    loo_terms,
    full_primary_m2,
    by = c("outcome", "term"),
    all.x = TRUE,
    sort = FALSE
  )

  loo_primary[, transmission_group := sub("^transmission_group", "", term)]

  loo_primary <- merge(
    loo_primary,
    mode_full_counts,
    by = "transmission_group",
    all.x = TRUE,
    sort = FALSE
  )

  # Number of removed TaxIDs belonging to the transmission mode represented by
  # each coefficient. This is composition information, not an outcome filter.
  removed_mode_counts <- five[, .N, by = .(
    virus_family = as.character(virus_family),
    transmission_group = as.character(transmission_group)
  )]
  setnames(removed_mode_counts, "N", "removed_mode_n")
  removed_mode_counts[, lineage_level := "virus_family"]
  setnames(removed_mode_counts, "virus_family", "removed_lineage")

  genus_removed_counts <- five[present(virus_genus), .N, by = .(
    virus_genus = as.character(virus_genus),
    transmission_group = as.character(transmission_group)
  )]
  setnames(genus_removed_counts, "N", "removed_mode_n")
  genus_removed_counts[, lineage_level := "virus_genus"]
  setnames(genus_removed_counts, "virus_genus", "removed_lineage")

  removed_mode_counts <- rbindlist(
    list(removed_mode_counts, genus_removed_counts),
    fill = TRUE
  )

  loo_primary <- merge(
    loo_primary,
    removed_mode_counts,
    by = c("lineage_level", "removed_lineage", "transmission_group"),
    all.x = TRUE,
    sort = FALSE
  )
  loo_primary[is.na(removed_mode_n), removed_mode_n := 0L]

  loo_primary[, `:=`(
    removed_mode_fraction = removed_mode_n / pmax(full_mode_n, 1),
    delta_estimate = estimate - full_estimate,
    abs_delta_estimate = abs(estimate - full_estimate),
    relative_abs_change = fifelse(
      is.finite(full_estimate) & abs(full_estimate) > 1e-12,
      abs(estimate - full_estimate) / abs(full_estimate),
      NA_real_
    ),
    sign_changed = is.finite(full_estimate) & is.finite(estimate) &
      sign(full_estimate) != sign(estimate),
    significance_lost_at_0_05 = is.finite(full_p_value) & full_p_value < 0.05 &
      is.finite(p_value) & p_value >= 0.05,
    significance_gained_at_0_05 = is.finite(full_p_value) & full_p_value >= 0.05 &
      is.finite(p_value) & p_value < 0.05
  )]

  loo_primary[, effect_scale := fifelse(
    model_type == "quasi-Poisson(log)",
    "rate ratio vs no documented vector relationship",
    "difference vs no documented vector relationship"
  )]

  loo_primary[, `:=`(
    full_reported_effect = fifelse(
      model_type == "quasi-Poisson(log)",
      exp(full_estimate),
      full_estimate
    ),
    leaveout_reported_effect = fifelse(
      model_type == "quasi-Poisson(log)",
      exp(estimate),
      estimate
    )
  )]

  # Rank impact independently within each lineage level, outcome, and
  # transmission-mode coefficient. Ranking is based only on the absolute change
  # on the fitted model scale; no p-value is used to select lineages.
  loo_primary[, impact_rank := frank(
    -abs_delta_estimate,
    ties.method = "min",
    na.last = "keep"
  ), by = .(
    lineage_level,
    outcome,
    transmission_group
  )]

  setorder(
    loo_primary,
    lineage_level,
    outcome,
    transmission_group,
    impact_rank,
    removed_lineage
  )

  write_tsv(
    loo_primary,
    file.path(tab_dir, "27_04D_lineage_leave_one_out_primary_contrasts.tsv")
  )

  loo_top <- loo_primary[impact_rank <= 10]
  write_tsv(
    loo_top,
    file.path(tab_dir, "34_04D_lineage_driver_ranked_summary.tsv")
  )

} else {
  loo_primary <- data.table()
  loo_top <- data.table()

  write_tsv(
    loo_primary,
    file.path(tab_dir, "27_04D_lineage_leave_one_out_primary_contrasts.tsv")
  )
  write_tsv(
    loo_top,
    file.path(tab_dir, "34_04D_lineage_driver_ranked_summary.tsv")
  )
}


# -----------------------------------------------------------------------------
# Within-mode lineage influence sensitivity
# -----------------------------------------------------------------------------
#
# The global leave-one-lineage-out analysis above removes a lineage from every
# transmission group, including the reference group. That analysis is useful for
# overall lineage-composition sensitivity but a large coefficient change can be
# caused by changing the composition of the reference group even when the
# lineage is absent from the target documented mode.
#
# This complementary analysis therefore evaluates only lineages that are
# actually present in a documented transmission mode. For each non-reference
# mode and each family/genus represented in that mode, only the TaxIDs belonging
# to that lineage x mode combination are removed. All reference-group records
# and the same lineage in any other mode are retained. The same fully adjusted
# models are then refitted. This provides a direct sensitivity analysis of which
# lineages within a transmission mode most influence that mode's coefficient,
# without altering reference-group composition.
# -----------------------------------------------------------------------------

DOCUMENTED_MODE_LEVELS <- MODE_LEVELS[MODE_LEVELS != "No documented vector relationship"]

within_mode_fit_one <- function(
    d,
    lineage_level,
    lineage_name,
    lineage_col,
    target_mode
) {
  in_target <- present(d[[lineage_col]]) &
    as.character(d[[lineage_col]]) == lineage_name &
    as.character(d$transmission_group) == target_mode

  removed <- d[in_target]
  if (!nrow(removed)) {
    return(list(audit = data.table(), terms = data.table()))
  }

  dd <- prep_model_data(d[!in_target])
  observed_after <- as.character(unique(dd$transmission_group))
  all_groups_retained <- all(MODE_LEVELS %in% observed_after)

  audit_row <- data.table(
    lineage_level = lineage_level,
    target_transmission_group = target_mode,
    removed_lineage = lineage_name,
    removed_target_mode_n = nrow(removed),
    full_target_mode_n = d[as.character(transmission_group) == target_mode, .N],
    remaining_target_mode_n = dd[as.character(transmission_group) == target_mode, .N],
    remaining_n = nrow(dd),
    all_five_groups_retained = all_groups_retained,
    eligible = all_groups_retained,
    exclusion_reason = if (all_groups_retained) "" else "removal_eliminated_target_transmission_group"
  )

  if (!all_groups_retained) {
    return(list(audit = audit_row, terms = data.table()))
  }

  out <- list()
  k_local <- 0L

  for (outcome in c("host_family_n", "host_order_n")) {
    z <- dd[
      is.finite(get(outcome)) &
        is.finite(qc_pass_VC80_support) &
        present(virus_family)
    ]

    if (nrow(z) >= 30 && all(MODE_LEVELS %in% as.character(unique(z$transmission_group)))) {
      k_local <- k_local + 1L
      out[[k_local]] <- fit_qp(
        as.formula(
          paste0(
            outcome,
            " ~ transmission_group + log1p(qc_pass_VC80_support)",
            " + factor(virus_family)"
          )
        ),
        z,
        paste0(outcome, "_M2_VC80_family_adjusted"),
        outcome,
        paste0(
          "remove_within_mode_", lineage_level, ":", lineage_name,
          "|mode:", target_mode
        )
      )
    }
  }

  z <- dd[
    is.finite(max_host_taxonomic_distance) &
      is.finite(qc_pass_VC80_support) &
      present(virus_family)
  ]

  if (nrow(z) >= 30 && all(MODE_LEVELS %in% as.character(unique(z$transmission_group)))) {
    k_local <- k_local + 1L
    out[[k_local]] <- fit_lm(
      max_host_taxonomic_distance ~
        transmission_group +
        log1p(qc_pass_VC80_support) +
        factor(virus_family),
      z,
      "max_distance_M2_VC80_family_adjusted",
      "max_host_taxonomic_distance",
      paste0(
        "remove_within_mode_", lineage_level, ":", lineage_name,
        "|mode:", target_mode
      )
    )
  }

  terms <- if (length(out)) {
    rbindlist(out, fill = TRUE)[term == paste0("transmission_group", target_mode)]
  } else {
    data.table()
  }

  if (nrow(terms)) {
    terms[, `:=`(
      lineage_level = lineage_level,
      target_transmission_group = target_mode,
      removed_lineage = lineage_name,
      removed_target_mode_n = nrow(removed),
      full_target_mode_n = d[as.character(transmission_group) == target_mode, .N],
      remaining_target_mode_n = dd[as.character(transmission_group) == target_mode, .N]
    )]
  }

  list(audit = audit_row, terms = terms)
}

within_audit_rows <- list()
within_term_rows <- list()
wa <- 0L
wt <- 0L

for (lineage_level in names(lineage_specs)) {
  lineage_col <- lineage_specs[[lineage_level]]

  for (target_mode in DOCUMENTED_MODE_LEVELS) {
    lineage_values <- sort(unique(std_text(
      five[as.character(transmission_group) == target_mode][[lineage_col]]
    )))
    lineage_values <- lineage_values[present(lineage_values)]

    for (lineage_name in lineage_values) {
      res <- within_mode_fit_one(
        five,
        lineage_level = lineage_level,
        lineage_name = lineage_name,
        lineage_col = lineage_col,
        target_mode = target_mode
      )

      if (nrow(res$audit)) {
        wa <- wa + 1L
        within_audit_rows[[wa]] <- res$audit
      }

      if (nrow(res$terms)) {
        wt <- wt + 1L
        within_term_rows[[wt]] <- res$terms
      }
    }
  }
}

within_audit <- if (length(within_audit_rows)) {
  rbindlist(within_audit_rows, fill = TRUE)
} else data.table()

within_terms <- if (length(within_term_rows)) {
  rbindlist(within_term_rows, fill = TRUE)
} else data.table()

write_tsv(
  within_audit,
  file.path(tab_dir, "38_04D_within_mode_lineage_eligibility_audit.tsv")
)

if (nrow(within_terms)) {
  within_primary <- merge(
    within_terms,
    full_primary_m2,
    by = c("outcome", "term"),
    all.x = TRUE,
    sort = FALSE
  )

  within_primary[, `:=`(
    removed_mode_fraction = removed_target_mode_n / pmax(full_target_mode_n, 1),
    delta_estimate = estimate - full_estimate,
    abs_delta_estimate = abs(estimate - full_estimate),
    relative_abs_change = fifelse(
      is.finite(full_estimate) & abs(full_estimate) > 1e-12,
      abs(estimate - full_estimate) / abs(full_estimate),
      NA_real_
    ),
    sign_changed = is.finite(full_estimate) & is.finite(estimate) &
      sign(full_estimate) != sign(estimate),
    significance_lost_at_0_05 = is.finite(full_p_value) & full_p_value < 0.05 &
      is.finite(p_value) & p_value >= 0.05,
    significance_gained_at_0_05 = is.finite(full_p_value) & full_p_value >= 0.05 &
      is.finite(p_value) & p_value < 0.05
  )]

  # Descriptive sensitivity classification. This is not a causal decomposition.
  within_primary[, removal_effect_on_full_coefficient := fifelse(
    sign_changed,
    "reverses coefficient direction",
    fifelse(
      abs(estimate) + 1e-12 < abs(full_estimate),
      "reduces coefficient magnitude",
      fifelse(
        abs(estimate) > abs(full_estimate) + 1e-12,
        "increases coefficient magnitude",
        "minimal coefficient change"
      )
    )
  )]

  within_primary[, effect_scale := fifelse(
    model_type == "quasi-Poisson(log)",
    "rate ratio vs no documented vector relationship",
    "difference vs no documented vector relationship"
  )]

  within_primary[, `:=`(
    full_reported_effect = fifelse(
      model_type == "quasi-Poisson(log)",
      exp(full_estimate),
      full_estimate
    ),
    leaveout_reported_effect = fifelse(
      model_type == "quasi-Poisson(log)",
      exp(estimate),
      estimate
    )
  )]

  within_primary[, influence_rank := frank(
    -abs_delta_estimate,
    ties.method = "min",
    na.last = "keep"
  ), by = .(
    lineage_level,
    outcome,
    target_transmission_group
  )]

  setorder(
    within_primary,
    lineage_level,
    outcome,
    target_transmission_group,
    influence_rank,
    removed_lineage
  )

  write_tsv(
    within_primary,
    file.path(tab_dir, "36_04D_within_mode_lineage_leaveout_primary_contrasts.tsv")
  )

  within_top <- within_primary[influence_rank <= 10]
  write_tsv(
    within_top,
    file.path(tab_dir, "37_04D_within_mode_lineage_influence_ranked_summary.tsv")
  )
} else {
  within_primary <- data.table()
  within_top <- data.table()
  write_tsv(
    within_primary,
    file.path(tab_dir, "36_04D_within_mode_lineage_leaveout_primary_contrasts.tsv")
  )
  write_tsv(
    within_top,
    file.path(tab_dir, "37_04D_within_mode_lineage_influence_ranked_summary.tsv")
  )
}


# -----------------------------------------------------------------------------
# Geographic breadth: transmission group
# -----------------------------------------------------------------------------

geo_terms <- list()
kg <- 0L

g <- five[
  is.finite(country_n) &
    country_n >= 0 &
    is.finite(association_n) &
    association_n >= 0 &
    is.finite(qc_pass_VC80_support)
]

if (nrow(g) > 0 && uniqueN(g$transmission_group) >= 2) {
  kg <- kg + 1L
  geo_terms[[kg]] <- fit_qp(
    country_n ~ transmission_group,
    g,
    "country_M0_crude",
    "country_n",
    "all_data"
  )

  kg <- kg + 1L
  geo_terms[[kg]] <- fit_qp(
    country_n ~ transmission_group + log1p(association_n),
    g,
    "country_M1_association_adjusted",
    "country_n",
    "all_data"
  )

  g2 <- g[present(virus_family)]
  if (nrow(g2) > 0 && uniqueN(g2$transmission_group) >= 2) {
    kg <- kg + 1L
    geo_terms[[kg]] <- fit_qp(
      country_n ~
        transmission_group +
        log1p(association_n) +
        factor(virus_family),
      g2,
      "country_M2_association_family_adjusted",
      "country_n",
      "all_data"
    )

    kg <- kg + 1L
    geo_terms[[kg]] <- fit_qp(
      country_n ~
        transmission_group +
        log1p(qc_pass_VC80_support) +
        factor(virus_family),
      g2,
      "country_M3_VC80_family_sensitivity",
      "country_n",
      "all_data"
    )

    kg <- kg + 1L
    geo_terms[[kg]] <- fit_qp(
      country_n ~
        transmission_group +
        log1p(association_n) +
        log1p(qc_pass_VC80_support) +
        factor(virus_family),
      g2,
      "country_M4_association_VC80_family_sensitivity",
      "country_n",
      "all_data"
    )
  }
}

geo_model_terms <- rbindlist(geo_terms, fill = TRUE)
geo_primary <- primary_terms(geo_model_terms, "transmission_group")
geo_primary[, effect_scale := "rate ratio vs no documented vector relationship"]

write_tsv(
  geo_model_terms,
  file.path(tab_dir, "28_04D_country_breadth_model_terms.tsv")
)

write_tsv(
  geo_primary,
  file.path(tab_dir, "29_04D_country_breadth_primary_contrasts.tsv")
)


# -----------------------------------------------------------------------------
# S0-S2 arthropod-associated sequence evidence vs VIP-DB documented vector
# -----------------------------------------------------------------------------

sv <- metric[
  !is.na(seq_supported) &
    !is.na(documented_vector) &
    present(virus_family) &
    is.finite(qc_pass_VC80_support)
]

sv[, `:=`(
  seq_supported_factor = factor(
    seq_supported,
    levels = c(FALSE, TRUE),
    labels = c("No sequence evidence", "Arthropod-associated sequence evidence")
  ),
  documented_vector_binary = as.integer(documented_vector)
)]

fisher_tab <- table(
  seq_supported = factor(sv$seq_supported, levels = c(FALSE, TRUE)),
  documented_vector = factor(sv$documented_vector, levels = c(FALSE, TRUE))
)

fisher_fit <- fisher.test(fisher_tab)

fisher_out <- data.table(
  test = "Fisher exact",
  no_sequence_no_documented_vector = unname(fisher_tab[1, 1]),
  no_sequence_documented_vector = unname(fisher_tab[1, 2]),
  sequence_evidence_no_documented_vector = unname(fisher_tab[2, 1]),
  sequence_evidence_documented_vector = unname(fisher_tab[2, 2]),
  odds_ratio = unname(fisher_fit$estimate),
  ci_low = unname(fisher_fit$conf.int[1]),
  ci_high = unname(fisher_fit$conf.int[2]),
  p_value = fisher_fit$p.value
)

write_tsv(
  fisher_out,
  file.path(tab_dir, "30_04D_S0S2_vs_VIPDB_raw_fisher.tsv")
)

sv_logit <- fit_logistic(
  documented_vector_binary ~
    seq_supported_factor +
    log1p(qc_pass_VC80_support) +
    factor(virus_family),
  sv,
  "S0S2_vs_documented_vector_VC80_family_adjusted",
  "VIPDB_documented_vector",
  "all_data"
)

sv_logit_primary <- sv_logit[
  !is.na(term) &
    startsWith(term, "seq_supported_factor")
]
sv_logit_primary[, effect_scale := "odds ratio for documented vector relationship"]

write_tsv(
  sv_logit,
  file.path(tab_dir, "31_04D_S0S2_vs_VIPDB_adjusted_model_terms.tsv")
)

write_tsv(
  sv_logit_primary,
  file.path(tab_dir, "32_04D_S0S2_vs_VIPDB_adjusted_primary.tsv")
)


# -----------------------------------------------------------------------------
# S0-S2 and documented-vector status vs country breadth
# -----------------------------------------------------------------------------

geo_binary_rows <- list()
kb <- 0L

for (exposure in c("seq_supported_factor", "documented_vector")) {
  d <- copy(metric)[
    is.finite(country_n) &
      is.finite(association_n) &
      is.finite(qc_pass_VC80_support) &
      present(virus_family)
  ]

  if (exposure == "seq_supported_factor") {
    d[, seq_supported_factor := factor(
      seq_supported,
      levels = c(FALSE, TRUE),
      labels = c(
        "No sequence evidence",
        "Arthropod-associated sequence evidence"
      )
    )]
  } else {
    d[, documented_vector := factor(
      documented_vector,
      levels = c(FALSE, TRUE),
      labels = c(
        "No documented vector relationship",
        "Documented vector relationship"
      )
    )]
  }

  if (uniqueN(d[[exposure]]) >= 2) {
    kb <- kb + 1L
    geo_binary_rows[[kb]] <- fit_qp(
      as.formula(
        paste0(
          "country_n ~ ", exposure,
          " + log1p(association_n) + factor(virus_family)"
        )
      ),
      d,
      paste0("country_by_", exposure, "_association_family"),
      "country_n",
      "all_data"
    )

    kb <- kb + 1L
    geo_binary_rows[[kb]] <- fit_qp(
      as.formula(
        paste0(
          "country_n ~ ", exposure,
          " + log1p(association_n)",
          " + log1p(qc_pass_VC80_support)",
          " + factor(virus_family)"
        )
      ),
      d,
      paste0("country_by_", exposure, "_association_VC80_family"),
      "country_n",
      "all_data"
    )
  }
}

geo_binary_terms <- rbindlist(geo_binary_rows, fill = TRUE)
geo_binary_primary <- geo_binary_terms[
  !is.na(term) &
    (
      startsWith(term, "seq_supported_factor") |
        startsWith(term, "documented_vector")
    )
]
geo_binary_primary[, effect_scale := "country-count rate ratio"]

write_tsv(
  geo_binary_primary,
  file.path(tab_dir, "33_04D_binary_evidence_country_breadth_primary.tsv")
)


# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

mode_counts <- five[, .N, by = transmission_group]
mode_count_text <- paste(
  paste(mode_counts$transmission_group, mode_counts$N, sep = "="),
  collapse = "; "
)

multiple_mode_vector_n <- metric[
  documented_vector %in% TRUE &
    !(
      fine4_status == "single_four_class_mode" &
        present(canonical_four_mode)
    ),
  .N
]

families_multi_mode <- family_overlap[
  transmission_group_n >= 2,
  .N
]

main_mode_models <- mode_primary[
  grepl("_M2_VC80_family_adjusted$", model_id)
]

loo_family_tested_n <- if (nrow(loo_audit)) {
  loo_audit[lineage_level == "virus_family" & eligible %in% TRUE, .N]
} else 0L

loo_genus_tested_n <- if (nrow(loo_audit)) {
  loo_audit[lineage_level == "virus_genus" & eligible %in% TRUE, .N]
} else 0L

loo_primary_row_n <- nrow(loo_primary)
loo_top_row_n <- nrow(loo_top)
within_primary_row_n <- nrow(within_primary)
within_top_row_n <- nrow(within_top)
within_eligible_n <- if (nrow(within_audit)) within_audit[eligible %in% TRUE, .N] else 0L

qc <- data.table(
  metric = c(
    "formal_taxid_metric_rows",
    "unique_virus_taxid_n",
    "five_group_included_n",
    "documented_vector_ambiguous_or_unclassified_mode_n",
    "five_group_counts",
    "families_with_at_least_two_five_groups",
    "leave_one_family_out_tested_n",
    "leave_one_genus_out_tested_n",
    "leave_one_lineage_primary_contrast_rows",
    "ranked_top_impact_rows",
    "within_mode_lineage_eligible_tests_n",
    "within_mode_lineage_primary_contrast_rows",
    "within_mode_ranked_top_impact_rows",
    "main_family_adjusted_mode_contrast_rows",
    "S0S2_VIPDB_adjusted_primary_rows",
    "vector_system_rows",
    "status"
  ),
  value = as.character(c(
    nrow(metric),
    uniqueN(metric$virus_taxid),
    nrow(five),
    multiple_mode_vector_n,
    mode_count_text,
    families_multi_mode,
    loo_family_tested_n,
    loo_genus_tested_n,
    loo_primary_row_n,
    loo_top_row_n,
    within_eligible_n,
    within_primary_row_n,
    within_top_row_n,
    nrow(main_mode_models),
    nrow(sv_logit_primary),
    nrow(vector_system),
    "PASS"
  )),
  note = c(
    metric_file,
    "",
    "No documented vector relationship plus single documented four-class modes only.",
    "Audited but excluded from five-group inferential models.",
    "",
    "Shows how much mode information overlaps within virus families.",
    "All observed virus families were considered; a lineage was skipped only if its removal eliminated at least one of the five transmission groups.",
    "All observed virus genera were considered; a lineage was skipped only if its removal eliminated at least one of the five transmission groups.",
    "Systematic outcome-independent leave-one-lineage-out sensitivity results.",
    "Top 10 global lineage-composition sensitivities per lineage level x outcome x transmission mode; ranking does not use p-values.",
    "Within-mode tests remove only lineage members belonging to the target documented transmission mode; reference-group composition is unchanged.",
    "Complete within-mode lineage influence contrasts.",
    "Top 10 within-mode lineage influences per lineage level x outcome x target mode; ranking does not use p-values.",
    "Expected >0 if family-adjusted five-group models are estimable.",
    "Expected 1 primary sequence-evidence term.",
    "Unique virus TaxID x vector family x virus genus x canonical mode combinations.",
    "04D v3 finished without modifying frozen 04A-04C outputs."
  )
)

write_tsv(
  qc,
  file.path(qc_dir, "03_Fig4_04D_transmission_biology_QC.tsv")
)


# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

cat("\nFig4 04D transmission biology extension PASS\n")
cat("Base:", base, "\n")
cat("Five-group virus TaxIDs:", nrow(five), "\n\n")

cat("Five-group counts:\n")
print(mode_counts)

cat("\nFamilies containing >=2 five-group levels:",
    families_multi_mode, "\n")

cat("\nTop family composition:\n")
print(head(family_comp, 20))

cat("\nTop genus composition:\n")
print(head(genus_comp, 20))

cat("\nS0-S2 vs documented-vector raw Fisher:\n")
print(fisher_out)

cat("\nS0-S2 vs documented-vector adjusted primary term:\n")
print(sv_logit_primary)

cat("\nSystematic global leave-one-lineage-out top impacts (first 20 rows):\n")
if (nrow(loo_top)) print(head(loo_top, 20))

cat("\nWithin-mode lineage influence top impacts (first 20 rows):\n")
if (nrow(within_top)) print(head(within_top, 20))

cat("\nMain outputs written to Fig4 tables 19-38 and QC/03.\n")
