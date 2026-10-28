suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

# Released inputs are stored gzipped, and some are exposed under both their
# historical stage-number name and a plain name via symlinks. Resolve `X.tsv` to
# `X.tsv.gz` when only the gzipped file exists, so the published script keeps the
# exact file names it used on the cluster.
resolve_release_input <- function(path) {
  if (!file.exists(path) && file.exists(paste0(path, ".gz"))) return(paste0(path, ".gz"))
  path
}

read_current_tsv <- function(path) {
  path <- assert_under_final_root(resolve_release_input(path), mustWork = TRUE)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

read_fig3_table <- function(filename) {
  read_current_tsv(file.path(TABLE_DIR, filename))
}

require_columns <- function(dat, cols, label) {
  miss <- setdiff(cols, names(dat))
  if (length(miss)) {
    stop(
      label, " missing required columns: ", paste(miss, collapse = ", "),
      "\nAvailable columns: ", paste(names(dat), collapse = ", ")
    )
  }
}

present <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) & nzchar(y) &
    !(tolower(y) %in% c("na", "n/a", "nan", "none", "null", "."))
}

as_bool <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

canonical_mode <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- gsub("[_/-]+", " ", s)
  s <- gsub("\\s+", " ", s)
  out <- rep(NA_character_, length(s))

  out[
    grepl("persistent", s) &
      grepl("propagative", s) &
      !grepl("non propagative|nonpropagative", s)
  ] <- "Persistent propagative"

  out[
    grepl("persistent", s) &
      grepl("non propagative|nonpropagative|circulative", s)
  ] <- "Persistent non-propagative"

  out[grepl("semi persistent|semipersistent", s)] <- "Semi-persistent"
  out[
    grepl("non persistent|nonpersistent", s) & !grepl("semi", s)
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
DOCUMENTED_MODE_LEVELS <- MODE_LEVELS[-1]

std_seq_group <- function(x) {
  z <- tolower(trimws(as.character(x)))
  case_when(
    str_detect(z, "no.*sequence.*support|^false$|^0$") ~ "No sequence support",
    str_detect(z, "sequence.*support|^true$|^1$") ~ "Sequence-supported",
    TRUE ~ NA_character_
  )
}

std_outcome <- function(x) {
  z <- tolower(trimws(as.character(x)))
  case_when(
    str_detect(z, "host[_ -]?family") ~ "Host family breadth",
    str_detect(z, "host[_ -]?order") ~ "Host order breadth",
    str_detect(z, "max.*taxonomic.*distance|max[_ -]?distance") ~
      "Maximum host taxonomic distance",
    str_detect(z, "country") ~ "Country breadth",
    TRUE ~ NA_character_
  )
}

outcome_display <- function(x) {
  recode(
    x,
    `Host family breadth` = "Host family breadth\nRate ratio",
    `Host order breadth` = "Host order breadth\nRate ratio",
    `Maximum host taxonomic distance` = "Maximum host taxonomic\ndistance\nDifference",
    `Country breadth` = "Country breadth\nRate ratio"
  )
}

format_p_compact <- function(p) {
  ifelse(
    is.na(p), "P = NA",
    ifelse(
      p < 0.001,
      paste0("P = ", formatC(p, format = "e", digits = 1)),
      paste0("P = ", formatC(p, format = "f", digits = 3))
    )
  )
}

# -----------------------------------------------------------------------------
# S0/S1/S2: derive from frozen final Stage-20 evidence table and cross-check union.
# -----------------------------------------------------------------------------
build_final_s012 <- function() {
  final_path <- file.path(PIPELINE_RESULTS, "20_arthropod_sequence_evidence_FINAL.tsv")
  dat <- read_current_tsv(final_path)

  id_candidates <- c("sequence_id", "plant_accession", "plant_accession_id", "accession", "plant_id")
  id_col <- id_candidates[id_candidates %in% names(dat)][1]
  code_candidates <- c("arthropod_evidence_codes", "evidence_codes", "arthropod_sequence_evidence_codes")
  code_col <- code_candidates[code_candidates %in% names(dat)][1]
  support_candidates <- c("arthropod_sequence_supported", "sequence_supported", "any_arthropod_sequence_support")
  support_col <- support_candidates[support_candidates %in% names(dat)][1]

  if (any(is.na(c(id_col, code_col, support_col)))) {
    stop("Cannot identify Stage-20 S0/S1/S2 columns in ", final_path)
  }

  z <- dat |>
    transmute(
      sequence_id = trimws(as.character(.data[[id_col]])),
      evidence_codes = as.character(.data[[code_col]]),
      sequence_supported = as_bool(.data[[support_col]]) %in% TRUE
    ) |>
    filter(present(sequence_id))

  count_code <- function(code) {
    z |>
      filter(
        !is.na(evidence_codes),
        str_detect(
          str_to_upper(evidence_codes),
          paste0("(^|[^A-Z0-9])", code, "([^A-Z0-9]|$)")
        )
      ) |>
      summarise(n = n_distinct(sequence_id)) |>
      pull(n)
  }

  out <- tibble(
    evidence = c(
      "S0 exact accession match",
      "S1 strict sequence similarity",
      "S2 cross-source 80% clustering",
      "S0–S2 union"
    ),
    value = c(
      count_code("S0"), count_code("S1"), count_code("S2"),
      z |> filter(sequence_supported) |> summarise(n = n_distinct(sequence_id)) |> pull(n)
    )
  )

  union_path <- file.path(PIPELINE_RESULTS, "13_arthropod_sequence_supported_plant_ids.tsv")
  union_dat <- read_current_tsv(union_path)
  union_id_candidates <- c("sequence_id", "plant_accession", "plant_accession_id", "accession", "plant_id", "id")
  union_id_col <- union_id_candidates[union_id_candidates %in% names(union_dat)][1]
  if (is.na(union_id_col)) stop("Cannot identify ID column in frozen Stage-13 union table.")
  union_n <- union_dat |>
    transmute(id = trimws(as.character(.data[[union_id_col]]))) |>
    filter(present(id)) |>
    summarise(n = n_distinct(id)) |>
    pull(n)

  if (out$value[out$evidence == "S0–S2 union"] != union_n) {
    stop("S0–S2 union mismatch between frozen Stage20 and Stage13 inputs.")
  }
  out
}

# -----------------------------------------------------------------------------
# Formal TaxID metrics and exact 04D-v3 five-group reconstruction.
# -----------------------------------------------------------------------------
prepare_formal_metrics <- function(t04) {
  require_columns(
    t04,
    c(
      "virus_taxid", "host_family_n", "host_order_n",
      "max_host_taxonomic_distance", "country_n",
      "any_arthropod_sequence_support", "qc_pass_VC80_support",
      "virus_family", "virus_genus", "documented_vector",
      "fine4_status", "fine4_mode"
    ),
    "Table 04"
  )

  # IMPORTANT: final 04D v3 did NOT use raw Table 04 taxonomy/geography alone.
  # Before building the five groups it merged taxonomy and geography from the
  # frozen Stage-20 masters, then backfilled missing virus family/genus and
  # non-finite country counts. Fig.S3A must reproduce that preprocessing
  # exactly; otherwise four reference-group TaxIDs with a missing family in
  # raw Table 04 are incorrectly dropped.
  seq_file <- file.path(PIPELINE_RESULTS, "20_three_source_sequence_master_FINAL.tsv")
  assoc_file <- file.path(PIPELINE_RESULTS, "20_three_source_association_master_FINAL.tsv")
  seq_file <- assert_under_final_root(resolve_release_input(seq_file), mustWork = TRUE)
  assoc_file <- assert_under_final_root(resolve_release_input(assoc_file), mustWork = TRUE)

  seq_cols <- c("virus_taxid", "virus_family", "virus_genus")
  assoc_cols <- c("virus_taxid", "country")

  seq_names <- names(data.table::fread(seq_file, nrows = 0, showProgress = FALSE))
  assoc_names <- names(data.table::fread(assoc_file, nrows = 0, showProgress = FALSE))
  miss_seq <- setdiff(seq_cols, seq_names)
  miss_assoc <- setdiff(assoc_cols, assoc_names)
  if (length(miss_seq)) {
    stop("Final sequence master missing required columns for exact 04D backfill: ",
         paste(miss_seq, collapse = ", "))
  }
  if (length(miss_assoc)) {
    stop("Final association master missing required columns for exact 04D backfill: ",
         paste(miss_assoc, collapse = ", "))
  }

  seq <- data.table::fread(
    seq_file, select = seq_cols, sep = "\t", showProgress = FALSE
  )
  assoc <- data.table::fread(
    assoc_file, select = assoc_cols, sep = "\t", showProgress = FALSE
  )
  metric <- data.table::as.data.table(t04)

  std_text_exact <- function(x) {
    y <- trimws(as.character(x))
    y[!present(y)] <- NA_character_
    y
  }
  safe_unique_one_exact <- function(x) {
    z <- unique(std_text_exact(x))
    z <- z[present(z)]
    if (length(z) == 1L) z else NA_character_
  }

  metric[, virus_taxid := std_text_exact(virus_taxid)]
  metric[, virus_family_raw_plotQC := std_text_exact(virus_family)]
  metric[, virus_genus_raw_plotQC := std_text_exact(virus_genus)]
  metric[, country_n_raw_plotQC := suppressWarnings(as.numeric(country_n))]

  seq[, virus_taxid := std_text_exact(virus_taxid)]
  assoc[, virus_taxid := std_text_exact(virus_taxid)]

  # Exact taxonomy fallback used by final 04D v3 (scientifically relevant
  # fields only; audit-count columns from 04D are not needed for plotting).
  tax_audit <- seq[present(virus_taxid), .(
    virus_family_final = safe_unique_one_exact(virus_family),
    virus_genus_final = safe_unique_one_exact(virus_genus)
  ), by = virus_taxid]

  geo <- assoc[present(virus_taxid), .(
    country_n_stage20 = data.table::uniqueN(std_text_exact(country)[present(country)])
  ), by = virus_taxid]

  metric <- merge(metric, tax_audit, by = "virus_taxid", all.x = TRUE, sort = FALSE)
  metric <- merge(metric, geo, by = "virus_taxid", all.x = TRUE, sort = FALSE)

  metric[, virus_family := data.table::fifelse(
    present(virus_family), std_text_exact(virus_family), virus_family_final
  )]
  metric[, virus_genus := data.table::fifelse(
    present(virus_genus), std_text_exact(virus_genus), virus_genus_final
  )]

  metric[, country_n := suppressWarnings(as.numeric(country_n))]
  metric[!is.finite(country_n), country_n := country_n_stage20]

  backfill_audit <- metric[
    (!present(virus_family_raw_plotQC) & present(virus_family)) |
      (!present(virus_genus_raw_plotQC) & present(virus_genus)) |
      (!is.finite(country_n_raw_plotQC) & is.finite(country_n)),
    .(
      virus_taxid,
      virus_family_raw = virus_family_raw_plotQC,
      virus_family_after_04D_backfill = virus_family,
      virus_genus_raw = virus_genus_raw_plotQC,
      virus_genus_after_04D_backfill = virus_genus,
      country_n_raw = country_n_raw_plotQC,
      country_n_after_04D_backfill = country_n,
      documented_vector = as_bool(documented_vector),
      fine4_status = std_text_exact(fine4_status),
      fine4_mode = std_text_exact(fine4_mode)
    )
  ]

  formal <- metric[, .(
    virus_taxid = std_text_exact(virus_taxid),
    virus_family = std_text_exact(virus_family),
    virus_genus = std_text_exact(virus_genus),
    host_family_n = suppressWarnings(as.numeric(host_family_n)),
    host_order_n = suppressWarnings(as.numeric(host_order_n)),
    max_host_taxonomic_distance = suppressWarnings(as.numeric(max_host_taxonomic_distance)),
    country_n = suppressWarnings(as.numeric(country_n)),
    qc_pass_VC80_support = suppressWarnings(as.numeric(qc_pass_VC80_support)),
    seq_support_group = ifelse(
      as_bool(any_arthropod_sequence_support) %in% TRUE,
      "Sequence-supported", "No sequence support"
    ),
    documented_vector = as_bool(documented_vector),
    fine4_status = std_text_exact(fine4_status),
    fine4_mode = std_text_exact(fine4_mode)
  )]
  formal[, canonical_four_mode := canonical_mode(fine4_mode)]

  list(
    formal = tibble::as_tibble(formal),
    backfill_audit = tibble::as_tibble(backfill_audit)
  )
}

build_exact_five_group <- function(formal) {
  z <- formal |>
    mutate(
      transmission_group = case_when(
        documented_vector %in% FALSE ~ "No documented vector relationship",
        documented_vector %in% TRUE &
          fine4_status == "single_four_class_mode" &
          present(canonical_four_mode) ~ canonical_four_mode,
        TRUE ~ NA_character_
      )
    ) |>
    filter(present(virus_family), present(transmission_group)) |>
    mutate(
      transmission_group = factor(transmission_group, levels = MODE_LEVELS),
      virus_family = factor(virus_family),
      virus_genus = factor(virus_genus)
    ) |>
    filter(!is.na(transmission_group))
  z
}

prepare_five_group_descriptives <- function(formal, t20) {
  five <- build_exact_five_group(formal)
  require_columns(
    t20,
    c(
      "transmission_group", "virus_taxid_n",
      "host_family_median", "host_family_IQR",
      "host_order_median", "host_order_IQR",
      "max_distance_median", "max_distance_IQR",
      "country_median", "country_IQR"
    ),
    "Table 20"
  )

  long <- five |>
    select(
      virus_taxid, transmission_group,
      host_family_n, host_order_n, max_host_taxonomic_distance, country_n
    ) |>
    pivot_longer(
      cols = c(host_family_n, host_order_n, max_host_taxonomic_distance, country_n),
      names_to = "outcome_key", values_to = "value"
    ) |>
    mutate(
      outcome = recode(
        outcome_key,
        host_family_n = "Host family breadth",
        host_order_n = "Host order breadth",
        max_host_taxonomic_distance = "Maximum host taxonomic distance",
        country_n = "Country breadth"
      )
    )

  summary <- long |>
    group_by(transmission_group, outcome) |>
    summarise(
      nonmissing_n = sum(is.finite(value)),
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE, type = 7),
      q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE, type = 7),
      iqr_width = IQR(value, na.rm = TRUE, type = 7),
      .groups = "drop"
    ) |>
    left_join(
      five |> count(transmission_group, name = "virus_taxid_n"),
      by = "transmission_group"
    ) |>
    mutate(
      transmission_group = as.character(transmission_group),
      outcome = factor(
        outcome,
        levels = c(
          "Host family breadth", "Host order breadth",
          "Maximum host taxonomic distance", "Country breadth"
        )
      )
    )

  expected <- bind_rows(
    t20 |> transmute(
      transmission_group = as.character(transmission_group), outcome = "Host family breadth",
      expected_n = as.numeric(virus_taxid_n), expected_median = as.numeric(host_family_median),
      expected_iqr = as.numeric(host_family_IQR)
    ),
    t20 |> transmute(
      transmission_group = as.character(transmission_group), outcome = "Host order breadth",
      expected_n = as.numeric(virus_taxid_n), expected_median = as.numeric(host_order_median),
      expected_iqr = as.numeric(host_order_IQR)
    ),
    t20 |> transmute(
      transmission_group = as.character(transmission_group), outcome = "Maximum host taxonomic distance",
      expected_n = as.numeric(virus_taxid_n), expected_median = as.numeric(max_distance_median),
      expected_iqr = as.numeric(max_distance_IQR)
    ),
    t20 |> transmute(
      transmission_group = as.character(transmission_group), outcome = "Country breadth",
      expected_n = as.numeric(virus_taxid_n), expected_median = as.numeric(country_median),
      expected_iqr = as.numeric(country_IQR)
    )
  )

  qc <- summary |>
    transmute(
      transmission_group, outcome = as.character(outcome),
      observed_n = virus_taxid_n,
      observed_median = median,
      observed_iqr = iqr_width
    ) |>
    full_join(expected, by = c("transmission_group", "outcome")) |>
    mutate(
      n_match = observed_n == expected_n,
      median_match = abs(observed_median - expected_median) < 1e-8,
      iqr_match = abs(observed_iqr - expected_iqr) < 1e-8,
      all_match = n_match & median_match & iqr_match
    )

  list(summary = summary, qc = qc, five = five)
}

# -----------------------------------------------------------------------------
# Binary evidence associations. Exact frozen primary rows only.
# -----------------------------------------------------------------------------
prepare_binary_effects <- function(t07, t33) {
  require_columns(t07, c("model", "term", "estimate", "std_error", "p_value"), "Table 07")
  require_columns(
    t33,
    c(
      "model_id", "term", "estimate", "std_error", "p_value",
      "exponentiated_effect", "exponentiated_ci_low", "exponentiated_ci_high"
    ),
    "Table 33"
  )

  specs <- tribble(
    ~model, ~term, ~exposure, ~outcome, ~scale,
    "host_family_by_sequence_evidence_quasipoisson", "seq_support_groupsequence-supported",
      "Arthropod-associated sequence evidence", "Host family breadth", "RR",
    "host_order_by_sequence_evidence_quasipoisson", "seq_support_groupsequence-supported",
      "Arthropod-associated sequence evidence", "Host order breadth", "RR",
    "max_distance_by_seq_support_group", "seq_support_groupsequence-supported",
      "Arthropod-associated sequence evidence", "Maximum host taxonomic distance", "Difference",
    "host_family_by_documented_vector_quasipoisson", "documented_vectorTRUE",
      "Documented vector relationship", "Host family breadth", "RR",
    "host_order_by_documented_vector_quasipoisson", "documented_vectorTRUE",
      "Documented vector relationship", "Host order breadth", "RR",
    "max_distance_by_documented_vector", "documented_vectorTRUE",
      "Documented vector relationship", "Maximum host taxonomic distance", "Difference"
  )

  host <- specs |>
    left_join(t07, by = c("model", "term"))
  if (any(!is.finite(host$estimate) | !is.finite(host$std_error))) {
    stop("Table 07 exact primary binary-evidence rows could not all be recovered.")
  }
  host <- host |>
    mutate(
      raw_lo = estimate - 1.96 * std_error,
      raw_hi = estimate + 1.96 * std_error,
      effect = ifelse(scale == "RR", exp(estimate), estimate),
      conf_low = ifelse(scale == "RR", exp(raw_lo), raw_lo),
      conf_high = ifelse(scale == "RR", exp(raw_hi), raw_hi)
    ) |>
    select(exposure, outcome, scale, effect, conf_low, conf_high, p_value, model)

  country_specs <- tribble(
    ~model_id, ~term, ~exposure,
    "country_by_seq_supported_factor_association_family",
      "seq_supported_factorArthropod-associated sequence evidence",
      "Arthropod-associated sequence evidence",
    "country_by_documented_vector_association_family",
      "documented_vectorDocumented vector relationship",
      "Documented vector relationship"
  )

  country <- country_specs |>
    left_join(t33, by = c("model_id", "term"))
  if (any(!is.finite(country$exponentiated_effect))) {
    stop("Table 33 primary association-count + family country rows could not both be recovered.")
  }
  country <- country |>
    transmute(
      exposure, outcome = "Country breadth", scale = "RR",
      effect = as.numeric(exponentiated_effect),
      conf_low = as.numeric(exponentiated_ci_low),
      conf_high = as.numeric(exponentiated_ci_high),
      p_value = as.numeric(p_value), model = model_id
    )

  bind_rows(host, country) |>
    mutate(
      exposure = factor(
        exposure,
        levels = c("Documented vector relationship", "Arthropod-associated sequence evidence")
      ),
      outcome = factor(
        outcome,
        levels = c(
          "Host family breadth", "Host order breadth",
          "Maximum host taxonomic distance", "Country breadth"
        )
      ),
      outcome_display = outcome_display(as.character(outcome)),
      null_value = ifelse(scale == "RR", 1, 0)
    )
}

prepare_country_sensitivity <- function(t33) {
  require_columns(
    t33,
    c(
      "model_id", "term", "p_value",
      "exponentiated_effect", "exponentiated_ci_low", "exponentiated_ci_high"
    ),
    "Table 33"
  )

  specs <- tribble(
    ~model_id, ~term, ~exposure, ~adjustment,
    "country_by_seq_supported_factor_association_family",
      "seq_supported_factorArthropod-associated sequence evidence",
      "Arthropod-associated sequence evidence", "Association count + family",
    "country_by_seq_supported_factor_association_VC80_family",
      "seq_supported_factorArthropod-associated sequence evidence",
      "Arthropod-associated sequence evidence", "Association count + VC80 + family",
    "country_by_documented_vector_association_family",
      "documented_vectorDocumented vector relationship",
      "Documented vector relationship", "Association count + family",
    "country_by_documented_vector_association_VC80_family",
      "documented_vectorDocumented vector relationship",
      "Documented vector relationship", "Association count + VC80 + family"
  )

  x <- specs |>
    left_join(t33, by = c("model_id", "term"))
  if (any(!is.finite(x$exponentiated_effect))) {
    stop("Table 33 country sensitivity rows could not all be recovered.")
  }
  x |>
    transmute(
      exposure,
      adjustment = factor(
        adjustment,
        levels = c("Association count + family", "Association count + VC80 + family")
      ),
      effect = as.numeric(exponentiated_effect),
      conf_low = as.numeric(exponentiated_ci_low),
      conf_high = as.numeric(exponentiated_ci_high),
      p_value = as.numeric(p_value)
    ) |>
    mutate(
      exposure = factor(
        exposure,
        levels = c("Documented vector relationship", "Arthropod-associated sequence evidence")
      )
    )
}

# -----------------------------------------------------------------------------
# S0-S2 vs VIP-DB documented-vector consistency: raw Fisher and adjusted logistic.
# -----------------------------------------------------------------------------
prepare_external_consistency <- function(t30, t32) {
  require_columns(t30, c("odds_ratio", "ci_low", "ci_high", "p_value"), "Table 30")
  require_columns(
    t32,
    c(
      "model_id", "term", "p_value",
      "exponentiated_effect", "exponentiated_ci_low", "exponentiated_ci_high"
    ),
    "Table 32"
  )
  adj <- t32 |>
    filter(
      model_id == "S0S2_vs_documented_vector_VC80_family_adjusted",
      term == "seq_supported_factorArthropod-associated sequence evidence"
    )
  if (nrow(adj) != 1) stop("Table 32 adjusted S0-S2 vs VIP-DB primary row is not unique.")

  bind_rows(
    tibble(
      model = "Raw Fisher",
      odds_ratio = as.numeric(t30$odds_ratio[1]),
      ci_low = as.numeric(t30$ci_low[1]),
      ci_high = as.numeric(t30$ci_high[1]),
      p_value = as.numeric(t30$p_value[1])
    ),
    tibble(
      model = "VC80 + family adjusted",
      odds_ratio = as.numeric(adj$exponentiated_effect[1]),
      ci_low = as.numeric(adj$exponentiated_ci_low[1]),
      ci_high = as.numeric(adj$exponentiated_ci_high[1]),
      p_value = as.numeric(adj$p_value[1])
    )
  ) |>
    mutate(model = factor(model, levels = c("VC80 + family adjusted", "Raw Fisher")))
}

# -----------------------------------------------------------------------------
# Host-distance observations and prediction curves from frozen Table 14.
# -----------------------------------------------------------------------------
prepare_distance_observed <- function(t13) {
  require_columns(t13, c("group", "taxonomy_distance", "sharing_probability"), "Table 13")
  t13 |>
    transmute(
      group = std_seq_group(group),
      taxonomy_distance = as.numeric(taxonomy_distance),
      sharing_probability = as.numeric(sharing_probability),
      total_pairs = if ("total_pairs" %in% names(t13)) as.numeric(total_pairs) else NA_real_,
      positive_pairs = if ("positive_pairs" %in% names(t13)) as.numeric(positive_pairs) else NA_real_
    ) |>
    filter(!is.na(group), is.finite(taxonomy_distance), is.finite(sharing_probability))
}

build_distance_predictions <- function(t14, observed) {
  require_columns(t14, c("term", "estimate"), "Table 14")
  ctab <- t14 |> transmute(term = as.character(term), estimate = as.numeric(estimate))
  pull_unique <- function(pattern, exclude = NULL) {
    h <- ctab |> filter(str_detect(term, regex(pattern, ignore_case = TRUE)))
    if (!is.null(exclude)) h <- h |> filter(!str_detect(term, regex(exclude, ignore_case = TRUE)))
    if (nrow(h) != 1) stop("Could not uniquely identify Table 14 coefficient: ", pattern)
    h$estimate[1]
  }
  b0 <- pull_unique("intercept")
  bdist <- pull_unique("taxonomy.*distance|host.*distance", exclude = ":|sequence|support")
  bgroup <- pull_unique("sequence|support", exclude = ":")
  hint <- ctab |>
    filter(str_detect(term, ":")) |>
    filter(str_detect(term, regex("distance", ignore_case = TRUE))) |>
    filter(str_detect(term, regex("sequence|support", ignore_case = TRUE)))
  if (nrow(hint) != 1) stop("Could not uniquely identify Table 14 distance x sequence-support term.")
  bint <- hint$estimate[1]

  xx <- seq(min(observed$taxonomy_distance), max(observed$taxonomy_distance), length.out = 200)
  bind_rows(
    tibble(group = "No sequence support", taxonomy_distance = xx, eta = b0 + bdist * xx),
    tibble(group = "Sequence-supported", taxonomy_distance = xx, eta = b0 + bgroup + (bdist + bint) * xx)
  ) |>
    mutate(sharing_probability = plogis(eta))
}

# -----------------------------------------------------------------------------
# Fully adjusted transmission-mode effects. Exact M2 primary rows.
# -----------------------------------------------------------------------------
prepare_mode_effects <- function(t26, t29, t20) {
  require_columns(
    t26,
    c(
      "model_id", "outcome", "term", "estimate", "p_value",
      "ci_low", "ci_high", "exponentiated_effect",
      "exponentiated_ci_low", "exponentiated_ci_high"
    ),
    "Table 26"
  )
  require_columns(
    t29,
    c(
      "model_id", "outcome", "term", "estimate", "p_value",
      "exponentiated_effect", "exponentiated_ci_low", "exponentiated_ci_high"
    ),
    "Table 29"
  )

  host_ids <- c(
    "host_family_n_M2_VC80_family_adjusted",
    "host_order_n_M2_VC80_family_adjusted",
    "max_distance_M2_VC80_family_adjusted"
  )
  host <- t26 |>
    filter(model_id %in% host_ids, str_starts(term, "transmission_group")) |>
    mutate(
      mode = str_remove(term, "^transmission_group"),
      outcome_std = std_outcome(outcome),
      scale = ifelse(outcome_std == "Maximum host taxonomic distance", "Difference", "RR"),
      effect = ifelse(scale == "RR", as.numeric(exponentiated_effect), as.numeric(estimate)),
      conf_low = ifelse(scale == "RR", as.numeric(exponentiated_ci_low), as.numeric(ci_low)),
      conf_high = ifelse(scale == "RR", as.numeric(exponentiated_ci_high), as.numeric(ci_high))
    ) |>
    select(mode, outcome = outcome_std, scale, effect, conf_low, conf_high, p_value, model_id)

  country <- t29 |>
    filter(
      model_id == "country_M2_association_family_adjusted",
      str_starts(term, "transmission_group")
    ) |>
    transmute(
      mode = str_remove(term, "^transmission_group"),
      outcome = "Country breadth", scale = "RR",
      effect = as.numeric(exponentiated_effect),
      conf_low = as.numeric(exponentiated_ci_low),
      conf_high = as.numeric(exponentiated_ci_high),
      p_value = as.numeric(p_value), model_id
    )

  out <- bind_rows(host, country)
  expected <- tidyr::crossing(
    mode = DOCUMENTED_MODE_LEVELS,
    outcome = c(
      "Host family breadth", "Host order breadth",
      "Maximum host taxonomic distance", "Country breadth"
    )
  )
  if (nrow(anti_join(expected, out, by = c("mode", "outcome"))) > 0) {
    stop("Table 26/29 do not contain all 16 exact fully adjusted primary contrasts.")
  }

  counts <- t20 |>
    select(transmission_group, virus_taxid_n) |>
    mutate(transmission_group = as.character(transmission_group))

  out |>
    left_join(counts, by = c("mode" = "transmission_group")) |>
    mutate(
      mode = factor(mode, levels = rev(DOCUMENTED_MODE_LEVELS)),
      outcome = factor(
        outcome,
        levels = c(
          "Host family breadth", "Host order breadth",
          "Maximum host taxonomic distance", "Country breadth"
        )
      ),
      outcome_display = outcome_display(as.character(outcome)),
      null_value = ifelse(scale == "RR", 1, 0)
    )
}

# -----------------------------------------------------------------------------
# Frozen RdRP65 sensitivity summary.
# -----------------------------------------------------------------------------
prepare_rdrp_summary <- function(t15) {
  require_columns(
    t15,
    c("rdrp_operational_cluster_id", "host_family_n", "host_order_n", "any_seq_support", "dataset"),
    "Table 15"
  )
  t15 |>
    transmute(
      lineage = as.character(rdrp_operational_cluster_id),
      dataset = as.character(dataset),
      support = ifelse(as_bool(any_seq_support) %in% TRUE, "Sequence-supported", "No sequence support"),
      host_family_n = as.numeric(host_family_n),
      host_order_n = as.numeric(host_order_n)
    ) |>
    pivot_longer(c(host_family_n, host_order_n), names_to = "outcome_key", values_to = "value") |>
    mutate(
      outcome = recode(
        outcome_key,
        host_family_n = "Host family breadth",
        host_order_n = "Host order breadth"
      )
    ) |>
    group_by(dataset, support, outcome) |>
    summarise(
      n = n(), median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

validate_rdrp_medians <- function(summary_raw, t16) {
  require_columns(t16, c("dataset", "any_seq_support", "host_family_median", "host_order_median"), "Table 16")
  expected <- t16 |>
    transmute(
      dataset = as.character(dataset),
      support = ifelse(as_bool(any_seq_support) %in% TRUE, "Sequence-supported", "No sequence support"),
      `Host family breadth` = as.numeric(host_family_median),
      `Host order breadth` = as.numeric(host_order_median)
    ) |>
    pivot_longer(c(`Host family breadth`, `Host order breadth`), names_to = "outcome", values_to = "expected_median")

  summary_raw |>
    select(dataset, support, outcome, observed_median = median) |>
    full_join(expected, by = c("dataset", "support", "outcome")) |>
    mutate(match = abs(observed_median - expected_median) < 1e-8)
}

# -----------------------------------------------------------------------------
# Family composition: read exact Table 21 fractions; select union of top 3
# families per documented mode for a compact heatmap.
# -----------------------------------------------------------------------------
prepare_family_composition <- function(t21) {
  require_columns(t21, c("transmission_group", "virus_family", "virus_taxid_n", "group_fraction"), "Table 21")
  z <- t21 |>
    filter(transmission_group %in% DOCUMENTED_MODE_LEVELS) |>
    mutate(
      transmission_group = as.character(transmission_group),
      virus_family = as.character(virus_family),
      group_fraction = as.numeric(group_fraction),
      virus_taxid_n = as.numeric(virus_taxid_n)
    )

  keep_fam <- z |>
    group_by(transmission_group) |>
    arrange(desc(group_fraction), .by_group = TRUE) |>
    slice_head(n = 3) |>
    ungroup() |>
    distinct(virus_family) |>
    pull(virus_family)

  fam_order <- z |>
    filter(virus_family %in% keep_fam) |>
    group_by(virus_family) |>
    summarise(max_fraction = max(group_fraction, na.rm = TRUE), .groups = "drop") |>
    arrange(max_fraction, virus_family) |>
    pull(virus_family)

  z |>
    filter(virus_family %in% keep_fam) |>
    select(transmission_group, virus_family, virus_taxid_n, group_fraction) |>
    complete(
      transmission_group = DOCUMENTED_MODE_LEVELS,
      virus_family = keep_fam,
      fill = list(virus_taxid_n = 0, group_fraction = 0)
    ) |>
    mutate(
      transmission_group = factor(transmission_group, levels = DOCUMENTED_MODE_LEVELS),
      virus_family = factor(virus_family, levels = fam_order)
    )
}

# -----------------------------------------------------------------------------
# Within-mode lineage influence: exact rank-1 virus-genus row per mode x host
# outcome. Country breadth was not part of Tables 36-38 and is not fabricated.
# -----------------------------------------------------------------------------
prepare_within_mode_influence <- function(t37) {
  require_columns(
    t37,
    c(
      "outcome", "lineage_level", "target_transmission_group",
      "removed_lineage", "removed_target_mode_n", "full_target_mode_n",
      "full_reported_effect", "leaveout_reported_effect",
      "relative_abs_change", "influence_rank"
    ),
    "Table 37"
  )

  out <- t37 |>
    filter(lineage_level == "virus_genus", as.numeric(influence_rank) == 1) |>
    mutate(
      outcome_std = std_outcome(outcome),
      mode = as.character(target_transmission_group)
    ) |>
    filter(
      mode %in% DOCUMENTED_MODE_LEVELS,
      outcome_std %in% c(
        "Host family breadth", "Host order breadth", "Maximum host taxonomic distance"
      )
    ) |>
    transmute(
      mode,
      outcome = outcome_std,
      removed_lineage = as.character(removed_lineage),
      removed_n = as.numeric(removed_target_mode_n),
      full_n = as.numeric(full_target_mode_n),
      full_effect = as.numeric(full_reported_effect),
      leaveout_effect = as.numeric(leaveout_reported_effect),
      relative_abs_change = as.numeric(relative_abs_change),
      scale = ifelse(outcome == "Maximum host taxonomic distance", "Difference", "RR")
    )

  expected <- tidyr::crossing(
    mode = DOCUMENTED_MODE_LEVELS,
    outcome = c("Host family breadth", "Host order breadth", "Maximum host taxonomic distance")
  )
  if (nrow(anti_join(expected, out, by = c("mode", "outcome"))) > 0) {
    stop("Table 37 is missing at least one exact rank-1 genus-level within-mode sensitivity row.")
  }

  out |>
    mutate(
      mode = factor(mode, levels = rev(DOCUMENTED_MODE_LEVELS)),
      outcome = factor(
        outcome,
        levels = c("Host family breadth", "Host order breadth", "Maximum host taxonomic distance")
      ),
      outcome_display = outcome_display(as.character(outcome)),
      null_value = ifelse(scale == "RR", 1, 0),
      lineage_label = paste0(
        removed_lineage, "\n",
        removed_n, "/", full_n, " TaxIDs; ",
        scales::percent(relative_abs_change, accuracy = 1)
      )
    )
}

# =============================================================================
# FINAL v6 display helpers
# =============================================================================

# VIP-DB is a separate evidence source from S0-S2. For the compact source
# panel, show the four single-mode groups used in 04D plus the audited
# ambiguous/unclassified documented-vector records that were excluded from
# five-group inference. Do not infer a transmission mode for excluded records.
prepare_vipdb_source_summary <- function(t20, t04) {
  require_columns(t20, c("transmission_group", "virus_taxid_n"), "Table 20")
  require_columns(t04, c("documented_vector", "fine4_status", "fine4_mode"), "Table 04")

  modes <- t20 |>
    transmute(
      category = as.character(transmission_group),
      value = as.numeric(virus_taxid_n)
    ) |>
    filter(category %in% DOCUMENTED_MODE_LEVELS)

  documented <- as_bool(t04$documented_vector)
  status <- trimws(as.character(t04$fine4_status))
  can_mode <- canonical_mode(t04$fine4_mode)
  excluded_n <- sum(
    documented %in% TRUE &
      !(status == "single_four_class_mode" & present(can_mode)),
    na.rm = TRUE
  )

  bind_rows(
    modes,
    tibble(category = "Ambiguous/unclassified", value = excluded_n)
  ) |>
    mutate(
      category = factor(
        category,
        levels = rev(c(
          "Non-persistent",
          "Semi-persistent",
          "Persistent non-propagative",
          "Persistent propagative",
          "Ambiguous/unclassified"
        ))
      ),
      colour_key = if_else(
        as.character(category) == "Ambiguous/unclassified",
        "Excluded",
        as.character(category)
      )
    )
}

# Minimal user-facing labels for the evidence consistency estimate.
prepare_external_consistency_v6 <- function(t30, t32) {
  x <- prepare_external_consistency(t30, t32) |>
    mutate(
      model = case_when(
        as.character(model) == "Raw Fisher" ~ "Unadjusted",
        TRUE ~ "Adjusted"
      ),
      model = factor(model, levels = c("Adjusted", "Unadjusted"))
    )
  x
}

# Combine the binary evidence layer and the documented transmission-mode layer
# into one forest matrix. The blank row is display-only and separates the two
# inferential questions. Reference groups differ by layer and must be stated in
# the figure legend/caption, not implied to be identical.
prepare_combined_effect_matrix <- function(binary_effects, mode_effects) {
  evidence <- binary_effects |>
    transmute(
      row_key = as.character(exposure),
      section = "Evidence status",
      outcome = as.character(outcome),
      outcome_display = as.character(outcome_display),
      scale,
      effect,
      conf_low,
      conf_high,
      p_value,
      null_value,
      colour_key = as.character(exposure)
    )

  modes <- mode_effects |>
    transmute(
      row_key = as.character(mode),
      section = "Transmission mode",
      outcome = as.character(outcome),
      outcome_display = as.character(outcome_display),
      scale,
      effect,
      conf_low,
      conf_high,
      p_value,
      null_value,
      colour_key = as.character(mode)
    )

  row_order_top_to_bottom <- c(
    "Arthropod-associated sequence evidence",
    "Documented vector relationship",
    " ",
    "Non-persistent",
    "Semi-persistent",
    "Persistent non-propagative",
    "Persistent propagative"
  )

  out <- bind_rows(evidence, modes) |>
    mutate(
      row_key = factor(row_key, levels = rev(row_order_top_to_bottom)),
      outcome = factor(
        outcome,
        levels = c(
          "Host family breadth",
          "Host order breadth",
          "Maximum host taxonomic distance",
          "Country breadth"
        )
      ),
      outcome_display = factor(
        .env$outcome_display(as.character(outcome)),
        levels = .env$outcome_display(c(
          "Host family breadth",
          "Host order breadth",
          "Maximum host taxonomic distance",
          "Country breadth"
        ))
      )
    )

  attr(out, "row_limits") <- rev(row_order_top_to_bottom)
  out
}

# Main-figure family-composition heatmap: keep only the two most abundant virus
# families per documented transmission mode (union across modes). This retains
# the dominant lineage structure while avoiding a supplementary-style catalog.
prepare_family_composition_main <- function(t21, top_n = 2L) {
  require_columns(
    t21,
    c("transmission_group", "virus_family", "virus_taxid_n", "group_fraction"),
    "Table 21"
  )

  z <- t21 |>
    filter(transmission_group %in% DOCUMENTED_MODE_LEVELS) |>
    mutate(
      transmission_group = as.character(transmission_group),
      virus_family = as.character(virus_family),
      group_fraction = as.numeric(group_fraction),
      virus_taxid_n = as.numeric(virus_taxid_n)
    )

  keep_fam <- z |>
    group_by(transmission_group) |>
    arrange(desc(group_fraction), .by_group = TRUE) |>
    slice_head(n = top_n) |>
    ungroup() |>
    distinct(virus_family) |>
    pull(virus_family)

  fam_order <- z |>
    filter(virus_family %in% keep_fam) |>
    group_by(virus_family) |>
    summarise(max_fraction = max(group_fraction, na.rm = TRUE), .groups = "drop") |>
    arrange(max_fraction, virus_family) |>
    pull(virus_family)

  z |>
    filter(virus_family %in% keep_fam) |>
    select(transmission_group, virus_family, virus_taxid_n, group_fraction) |>
    complete(
      transmission_group = DOCUMENTED_MODE_LEVELS,
      virus_family = keep_fam,
      fill = list(virus_taxid_n = 0, group_fraction = 0)
    ) |>
    mutate(
      transmission_group = factor(transmission_group, levels = DOCUMENTED_MODE_LEVELS),
      virus_family = factor(virus_family, levels = fam_order)
    )
}

# Same frozen rank-1 within-mode sensitivity as v5, but with a deliberately
# minimal label. Relative coefficient-change percentages remain in the TSV
# outputs and are not printed on the figure because they can be confused with
# within-mode proportions.
prepare_within_mode_influence_v6 <- function(t37) {
  x <- prepare_within_mode_influence(t37)
  x |>
    mutate(
      lineage_label = paste0(removed_lineage, " (", removed_n, "/", full_n, ")")
    )
}
