rm(list = ls())
options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_script <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
} else {
  normalizePath("03_Fig3_publication_plot_FINAL_v22.R", mustWork = TRUE)
}
SCRIPT_DIR <- dirname(this_script)

source(file.path(SCRIPT_DIR, "00_check_packages.R"))
source(file.path(SCRIPT_DIR, "00_display_rules.R"))
source(file.path(SCRIPT_DIR, "00_paths.R"))
source(file.path(SCRIPT_DIR, "00_plot_theme.R"))
source(file.path(SCRIPT_DIR, "00_helpers.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tibble)
})

log_file <- file.path(LOG_DIR, "Fig3_publication_plot_FINAL_v22.log")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
on.exit({ sink(); close(zz) }, add = TRUE)

cat("NEW Fig.3 publication plotting FINAL v22 started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("FINAL_THREE_SOURCE_ROOT: ", FINAL_THREE_SOURCE_ROOT, "\n", sep = "")
cat("PIPELINE_RESULTS: ", PIPELINE_RESULTS, "\n", sep = "")
cat("FIG3_ANALYSIS_ROOT: ", FIG3_ANALYSIS_ROOT, "\n", sep = "")
cat("Resolved figure font: ", FIG_FONT, "\n", sep = "")
cat("Typography: EXACTLY inherited from Fig4_publication_FINAL_v25/v26.\n")
cat("Scientific source: final 04D_transmission_biology_extended_v3 definitions/results.\n")
cat("v13 layout: A grouped horizontal lollipop summary; B OR + C 0-6 distance on same row; D documented-vector ecology full width; S3A bold x labels; S3B two outcome panels with offset genus labels.\n")
cat("R version: ", R.version.string, "\n\n", sep = "")

# =============================================================================
# Utility: derive virus-TaxID counts for S0/S1/S2/union from frozen inputs.
# =============================================================================
build_s012_taxid_counts <- function() {
  ev_path <- file.path(PIPELINE_RESULTS, "20_arthropod_sequence_evidence_FINAL.tsv")
  seq_path <- file.path(PIPELINE_RESULTS, "20_three_source_sequence_master_FINAL.tsv")
  ev <- read_current_tsv(ev_path)
  seqm <- read_current_tsv(seq_path)

  ev_id_candidates <- c("sequence_id", "plant_accession", "plant_accession_id", "accession", "plant_id")
  seq_id_candidates <- c("sequence_id", "plant_accession", "plant_accession_id", "accession", "plant_id")
  ev_id_col <- ev_id_candidates[ev_id_candidates %in% names(ev)][1]
  seq_id_col <- seq_id_candidates[seq_id_candidates %in% names(seqm)][1]
  code_candidates <- c("arthropod_evidence_codes", "evidence_codes", "arthropod_sequence_evidence_codes")
  support_candidates <- c("arthropod_sequence_supported", "sequence_supported", "any_arthropod_sequence_support")
  code_col <- code_candidates[code_candidates %in% names(ev)][1]
  support_col <- support_candidates[support_candidates %in% names(ev)][1]
  if (any(is.na(c(ev_id_col, seq_id_col, code_col, support_col))) || !"virus_taxid" %in% names(seqm)) {
    stop("Cannot derive S0/S1/S2 virus-TaxID counts from frozen Stage-20 inputs.")
  }

  ev2 <- ev |>
    transmute(
      sequence_id = trimws(as.character(.data[[ev_id_col]])),
      evidence_codes = as.character(.data[[code_col]]),
      sequence_supported = as_bool(.data[[support_col]]) %in% TRUE
    ) |>
    filter(present(sequence_id))

  seq2 <- seqm |>
    transmute(
      sequence_id = trimws(as.character(.data[[seq_id_col]])),
      virus_taxid = suppressWarnings(as.numeric(as.character(virus_taxid)))
    ) |>
    filter(present(sequence_id), !is.na(virus_taxid)) |>
    distinct(sequence_id, virus_taxid)

  j <- ev2 |>
    inner_join(seq2, by = "sequence_id")

  count_code_taxid <- function(code) {
    j |>
      filter(
        !is.na(evidence_codes),
        str_detect(str_to_upper(evidence_codes), paste0("(^|[^A-Z0-9])", code, "([^A-Z0-9]|$)"))
      ) |>
      summarise(n = n_distinct(virus_taxid)) |>
      pull(n)
  }

  tibble(
    row = c("S0", "S1", "S2", "S0–S2 union"),
    count = c(
      count_code_taxid("S0"),
      count_code_taxid("S1"),
      count_code_taxid("S2"),
      j |>
        filter(sequence_supported) |>
        summarise(n = n_distinct(virus_taxid)) |>
        pull(n)
    )
  )
}

# =============================================================================
# 1. Frozen inputs only
# =============================================================================
s012 <- build_final_s012()
s012_taxids <- build_s012_taxid_counts()
t04 <- read_fig3_table("04_formal_taxid_arthropod_vector_metrics.tsv")
t07 <- read_fig3_table("07_support_and_family_adjusted_models.tsv")
t13 <- read_fig3_table("13_host_distance_by_sequence_evidence_group.tsv")
t14 <- read_fig3_table("14_host_distance_interaction_model.tsv")
t20 <- read_fig3_table("20_04D_five_group_descriptive_summary.tsv")
t21 <- read_fig3_table("21_04D_transmission_group_by_virus_family.tsv")
t26 <- read_fig3_table("26_04D_five_group_host_range_primary_contrasts.tsv")
t29 <- read_fig3_table("29_04D_country_breadth_primary_contrasts.tsv")
t30 <- read_fig3_table("30_04D_S0S2_vs_VIPDB_raw_fisher.tsv")
t32 <- read_fig3_table("32_04D_S0S2_vs_VIPDB_adjusted_primary.tsv")
t33 <- read_fig3_table("33_04D_binary_evidence_country_breadth_primary.tsv")
t37 <- read_fig3_table("37_04D_within_mode_lineage_influence_ranked_summary.tsv")

formal_obj <- prepare_formal_metrics(t04)
formal <- formal_obj$formal
backfill_audit <- formal_obj$backfill_audit
external_consistency <- prepare_external_consistency_v6(t30, t32)
binary_effects <- prepare_binary_effects(t07, t33)
mode_effects <- prepare_mode_effects(t26, t29, t20)
distance_obs <- prepare_distance_observed(t13)
distance_pred <- build_distance_predictions(t14, distance_obs)
family_comp <- prepare_family_composition(t21)
lineage_influence <- prepare_within_mode_influence(t37)

write_tsv(backfill_audit, file.path(LOG_DIR, "Fig3_FINAL_v22_QC_04D_taxonomy_geography_backfill.tsv"))
write_tsv(s012, file.path(LOG_DIR, "Fig3_FINAL_v22_input_S012_sequence_counts.tsv"))
write_tsv(s012_taxids, file.path(LOG_DIR, "Fig3_FINAL_v22_input_S012_taxid_counts.tsv"))
write_tsv(external_consistency, file.path(LOG_DIR, "Fig3_FINAL_v22_input_evidence_consistency.tsv"))
write_tsv(binary_effects, file.path(LOG_DIR, "Fig3_FINAL_v22_input_binary_effects.tsv"))
write_tsv(mode_effects, file.path(LOG_DIR, "Fig3_FINAL_v22_input_mode_effects.tsv"))
write_tsv(distance_obs, file.path(LOG_DIR, "Fig3_FINAL_v22_input_distance_observed.tsv"))
write_tsv(family_comp, file.path(LOG_DIR, "FigS3_FINAL_v22_input_family_composition.tsv"))
write_tsv(lineage_influence, file.path(LOG_DIR, "FigS3_FINAL_v22_input_within_mode_influence.tsv"))

# =============================================================================
# 2. Derived summary counts for Fig.3A
# =============================================================================
seq_taxid_n <- formal |>
  filter(seq_support_group == "Sequence-supported") |>
  summarise(n = n_distinct(virus_taxid)) |>
  pull(n)

docvec_taxid_n <- formal |>
  filter(documented_vector %in% TRUE) |>
  summarise(n = n_distinct(virus_taxid)) |>
  pull(n)

s012_union_n <- s012 |>
  filter(evidence == "S0–S2 union") |>
  pull(value)

mode_count_map <- t20 |>
  transmute(mode = as.character(transmission_group), n = as.numeric(virus_taxid_n)) |>
  filter(mode %in% DOCUMENTED_MODE_LEVELS)
mode_n <- setNames(mode_count_map$n, mode_count_map$mode)

A_seq <- s012_taxids |>
  mutate(
    group = "Arthropod-associated\nsequence evidence",
    colour_key = if_else(row == "S0–S2 union", "S0–S2 union", "Arthropod-associated sequence evidence"),
    y = c(8, 7, 6, 5)
  )
A_doc <- tibble(
  row = c("Overall", DOCUMENTED_MODE_LEVELS),
  count = c(docvec_taxid_n, unname(mode_n[DOCUMENTED_MODE_LEVELS])),
  group = "Documented vector\nrelationship",
  colour_key = c("Documented vector relationship", DOCUMENTED_MODE_LEVELS),
  y = c(3.5, 2.5, 1.8, 1.1, 0.4)
)
A_dat <- bind_rows(A_seq, A_doc)

A_lab <- tibble(
  y = c(mean(A_seq$y), mean(A_doc$y)),
  label = c(
    "Arthropod-associated\nsequence evidence",
    "VIP-DB documented\nvector evidence"
  )
)

A_row_labels <- A_dat |>
  transmute(
    row, y,
    label_colour = case_when(
      row %in% c("S0", "S1", "S2") ~ "grey18",
      row == "S0–S2 union" ~ COL_BLUE_DARK,
      row == "Overall" ~ COL_ORANGE,
      TRUE ~ "grey18"
    ),
    label_face = if_else(row %in% c("S0–S2 union", "Overall"), "bold", "plain")
  )

# =============================================================================
# 3. Main Fig.3A: grouped labels with brackets and leader lines (union/overall excluded)
# =============================================================================
# Render the left information area entirely within its own panel so long group
# titles are never clipped. Use brackets to link each group title to the
# relevant detailed rows while leaving S0-S2 union and Overall outside the
# brackets, as requested.
A_brackets <- tibble(
  x = c(0.60, 0.60),
  x_tip = c(0.67, 0.67),
  y_top = c(8.25, 2.75),
  y_bot = c(5.75, 0.15),
  y_mid = c((8.25 + 5.75) / 2, (2.75 + 0.15) / 2)
)

A_group_text <- tibble(
  x = c(0.42, 0.42),
  y = c(7.0, 1.45),
  label = c(
    "Arthropod-associated\nsequence evidence",
    "VIP-DB documented\nvector evidence"
  )
)

A_leaders <- tibble(
  x = c(0.45, 0.45),
  xend = c(0.57, 0.57),
  y = A_brackets$y_mid,
  yend = A_brackets$y_mid
)

pA_info <- ggplot() +
  geom_text(
    data = A_group_text,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, family = FIG_FONT, size = TEXT_NORMAL / .pt,
    fontface = "bold", lineheight = 0.92, colour = "grey18"
  ) +
  geom_segment(
    data = A_leaders,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = LINE_THIN, colour = "grey55", lineend = "round"
  ) +
  geom_segment(
    data = A_brackets,
    aes(x = x, xend = x, y = y_bot, yend = y_top),
    inherit.aes = FALSE,
    linewidth = LINE_THIN, colour = "grey55", lineend = "round"
  ) +
  geom_segment(
    data = A_brackets,
    aes(x = x, xend = x_tip, y = y_top, yend = y_top),
    inherit.aes = FALSE,
    linewidth = LINE_THIN, colour = "grey55", lineend = "round"
  ) +
  geom_segment(
    data = A_brackets,
    aes(x = x, xend = x_tip, y = y_bot, yend = y_bot),
    inherit.aes = FALSE,
    linewidth = LINE_THIN, colour = "grey55", lineend = "round"
  ) +
  geom_text(
    data = A_row_labels,
    aes(x = 0.985, y = y, label = row, colour = label_colour, fontface = label_face),
    inherit.aes = FALSE,
    hjust = 1, family = FIG_FONT, size = TEXT_SMALL / .pt,
    lineheight = 0.92, show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(-0.35, 8.75), breaks = NULL,
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = FIG_FONT) +
  theme(plot.margin = margin(6, 6, 4, 6))

pA_data <- ggplot(A_dat, aes(x = count, y = y, colour = colour_key)) +
  geom_segment(aes(x = 0, xend = count, yend = y), linewidth = 0.80, lineend = "round") +
  geom_point(size = 2.45) +
  geom_text(
    aes(x = count + 14, label = scales::comma(count)),
    hjust = 0, family = FIG_FONT, size = TEXT_SMALL / .pt,
    colour = "grey18", show.legend = FALSE
  ) +
  scale_colour_manual(values = c(ASSOC_COLS, "S0–S2 union" = COL_BLUE_DARK), guide = "none") +
  scale_x_continuous(
    limits = c(0, 450),
    breaks = c(0, 100, 200, 300, 400),
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_y_continuous(limits = c(-0.35, 8.75), breaks = NULL,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Virus TaxIDs", y = NULL) +
  theme_publication() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 0.40),
    axis.line.y = element_line(colour = "black", linewidth = 0.40),
    plot.margin = margin(6, 8, 4, 0)
  )

pA <- pA_info + pA_data +
  plot_layout(widths = c(94, 86))

# =============================================================================
# 4. Main Fig.3B: evidence correspondence (left side of B+C row)
# =============================================================================
ext_plot <- external_consistency

pB <- ggplot(ext_plot, aes(x = odds_ratio, y = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.34) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.14, linewidth = 0.58, lineend = "round", colour = COL_GREY_DARK) +
  geom_point(size = 2.10, colour = COL_GREY_DARK) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4, 8),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.10, 0.12))
  ) +
  labs(x = "Odds ratio (95% CI)", y = NULL) +
  theme_publication() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = TEXT_SMALL),
    plot.margin = margin(6, 10, 4, 4)
  )

# =============================================================================
# 5. Main Fig.3C: 0-6 host-distance decay (right side of B+C row)
# =============================================================================
distance_obs <- distance_obs |>
  mutate(group = recode(group, "No sequence support" = "Absent", "Sequence-supported" = "Present"))
distance_pred <- distance_pred |>
  mutate(group = recode(group, "No sequence support" = "Absent", "Sequence-supported" = "Present"))

distance_levels <- c("Absent", "Present")
distance_obs$group <- factor(distance_obs$group, levels = distance_levels)
distance_pred$group <- factor(distance_pred$group, levels = distance_levels)

pC <- ggplot() +
  geom_line(data = distance_pred,
            aes(x = taxonomy_distance, y = sharing_probability, colour = group), linewidth = 0.92) +
  geom_point(data = distance_obs |> filter(sharing_probability > 0),
             aes(x = taxonomy_distance, y = sharing_probability, colour = group), size = 2.05) +
  scale_colour_manual(values = SEQ_COLS, name = "Arthropod-associated sequence evidence") +
  guides(colour = guide_legend(title.position = "left", nrow = 1, byrow = TRUE)) +
  scale_x_continuous(
    breaks = sort(unique(distance_obs$taxonomy_distance)),
    expand = expansion(mult = c(0.015, 0.02))
  ) +
  scale_y_log10(
    labels = function(x) {
      p <- x * 100
      ifelse(
        p >= 1,
        paste0(formatC(p, format = "f", digits = 0), "%"),
        ifelse(
          p >= 0.1,
          paste0(formatC(p, format = "f", digits = 1), "%"),
          paste0(formatC(p, format = "f", digits = 2), "%")
        )
      )
    },
    breaks = scales::breaks_log(n = 5),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(x = "Host taxonomic distance", y = "Virus sharing probability (%)") +
  theme_publication() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.height = grid::unit(8, "pt"),
    legend.key.width = grid::unit(11, "pt"),
    legend.spacing.x = grid::unit(2, "pt"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.20),
    plot.margin = margin(6, 4, 4, 4)
  )

# =============================================================================
# 6. Main Fig.3D: documented-vector ecology (full width, final row)
# =============================================================================
outcomes_keep <- c("Host family breadth", "Maximum host taxonomic distance", "Country breadth")
main_outcome_labels <- c(
  "Host family breadth Rate ratio",
  "Maximum host distance Difference",
  "Country breadth\nRate ratio"
)

docrow <- binary_effects |>
  filter(exposure == "Documented vector relationship", outcome %in% outcomes_keep) |>
  mutate(
    outcome_key = as.character(outcome),
    outcome_short = factor(outcome_key, levels = outcomes_keep, labels = main_outcome_labels)
  ) |>
  transmute(
    row_name = "Documented vector relationship",
    y = 6.0,
    outcome_short,
    null_value, effect, conf_low, conf_high,
    colour_key = "Documented vector relationship"
  )

moderows <- mode_effects |>
  filter(outcome %in% outcomes_keep) |>
  mutate(
    mode_chr = as.character(mode),
    outcome_key = as.character(outcome),
    outcome_short = factor(outcome_key, levels = outcomes_keep, labels = main_outcome_labels),
    y = case_when(
      mode_chr == "Non-persistent" ~ 4.0,
      mode_chr == "Semi-persistent" ~ 3.0,
      mode_chr == "Persistent non-propagative" ~ 2.0,
      mode_chr == "Persistent propagative" ~ 1.0,
      TRUE ~ NA_real_
    )
  ) |>
  transmute(
    row_name = mode_chr,
    y = as.numeric(y),
    outcome_short,
    null_value, effect, conf_low, conf_high,
    colour_key = mode_chr
  )

D_dat <- bind_rows(docrow, moderows)
null_df_D <- D_dat |> distinct(outcome_short, null_value)

y_breaks_D <- c(6, 4, 3, 2, 1)
y_labels_D <- c(
  "Documented vector\nrelationship",
  "Non-persistent",
  "Semi-persistent",
  "Persistent non-propagative",
  "Persistent propagative"
)

pD <- ggplot(D_dat, aes(x = effect, y = y, colour = colour_key)) +
  geom_vline(
    data = null_df_D,
    aes(xintercept = null_value),
    inherit.aes = FALSE,
    linetype = "dashed", colour = "grey72", linewidth = 0.34
  ) +
  geom_hline(yintercept = 5, colour = "grey84", linewidth = 0.36) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    orientation = "y", width = 0.14, linewidth = 0.58, lineend = "round"
  ) +
  geom_point(size = 2.10) +
  facet_wrap(~ outcome_short, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = ASSOC_COLS, guide = "none") +
  scale_y_continuous(
    breaks = y_breaks_D,
    labels = y_labels_D,
    limits = c(0.55, 6.35),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.10, 0.12))
  ) +
  labs(x = "Adjusted association vs no documented vector relationship (95% CI)", y = NULL) +
  theme_publication() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = TEXT_SMALL, lineheight = 0.92),
    axis.text.y = element_text(size = TEXT_SMALL, lineheight = 0.92),
    plot.margin = margin(6, 4, 4, 8)
  )

# =============================================================================
# 7. Save main panels and assemble: A / (B + C) / D
# =============================================================================
save_plot_all_formats(pA, "Fig3_panelA_evidence_summary_FINAL_v22", A_WIDTH_MM, A_HEIGHT_MM, PANEL_DIR)
save_plot_all_formats(pB, "Fig3_panelB_evidence_correspondence_FINAL_v22", BC_WIDTHS_MM[1], BC_HEIGHT_MM, PANEL_DIR)
save_plot_all_formats(pC, "Fig3_panelC_host_distance_decay_FINAL_v22", BC_WIDTHS_MM[2], BC_HEIGHT_MM, PANEL_DIR)
save_plot_all_formats(pD, "Fig3_panelD_documented_vector_ecology_FINAL_v22", D_WIDTH_MM, D_HEIGHT_MM, PANEL_DIR)

BC_row <- wrap_elements(full = pB) + wrap_elements(full = pC) + plot_layout(widths = BC_WIDTHS_MM)
main_plot <- wrap_elements(full = pA) / wrap_elements(full = BC_row) / wrap_elements(full = pD) +
  plot_layout(heights = c(A_HEIGHT_MM, BC_HEIGHT_MM, D_HEIGHT_MM))

save_plot_all_formats(main_plot, "Fig3_main_publication_FINAL_v22", MAIN_WIDTH_MM, MAIN_HEIGHT_MM, MAIN_DIR)
cat("Saved main Fig.3 and all standalone main panels.\n")

# =============================================================================
# 8. Fig.S3A: virus-family composition; no n labels; bold x-axis labels
# =============================================================================
pSA <- ggplot(family_comp, aes(x = transmission_group, y = virus_family, fill = group_fraction)) +
  geom_tile(colour = "white", linewidth = 0.32) +
  geom_text(
    data = family_comp |> filter(group_fraction >= 0.07),
    aes(label = label_percent(accuracy = 1)(group_fraction)),
    family = FIG_FONT, size = TEXT_SMALL / .pt, colour = "grey15"
  ) +
  scale_fill_gradient(
    low = "white", high = COL_BLUE_DARK,
    labels = label_percent(accuracy = 1), name = "Within mode\nproportion"
  ) +
  scale_x_discrete(labels = c(
    "Non-persistent" = "Non-persistent",
    "Semi-persistent" = "Semi-persistent",
    "Persistent non-propagative" = "Persistent non-propagative",
    "Persistent propagative" = "Persistent propagative"
  )) +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = TEXT_SMALL, lineheight = 0.92, face = "bold"),
    axis.text.y = element_text(size = TEXT_SMALL, face = "italic"),
    legend.position = "right",
    plot.margin = margin(4, 8, 4, 4)
  )

# =============================================================================
# 9. Fig.S3B: within-mode lineage influence; retain two outcomes only.
# =============================================================================
S3_outcome_levels <- c(
  "Host family breadth",
  "Maximum host taxonomic distance"
)
S3_outcome_labels <- c(
  "Host family breadth Rate ratio",
  "Maximum host distance Difference"
)

lineage_plot <- lineage_influence |>
  filter(outcome %in% S3_outcome_levels) |>
  mutate(
    outcome_key = as.character(outcome),
    outcome_short = factor(outcome_key, levels = S3_outcome_levels, labels = S3_outcome_labels)
  ) |>
  group_by(outcome_short) |>
  mutate(
    range_min = min(c(full_effect, leaveout_effect), na.rm = TRUE),
    range_max = max(c(full_effect, leaveout_effect), na.rm = TRUE),
    range_span = pmax(range_max - range_min, 0.25),
    endpoint = pmax(full_effect, leaveout_effect),
    label_x = endpoint + 0.12 * range_span,
    leader_x0 = endpoint + 0.03 * range_span,
    leader_x1 = endpoint + 0.09 * range_span
  ) |>
  ungroup()

null_SB <- lineage_plot |> distinct(outcome_short, null_value)
point_SB <- bind_rows(
  lineage_plot |> transmute(mode, outcome_short, state = "Full model", effect = full_effect),
  lineage_plot |> transmute(mode, outcome_short, state = "After lineage removal", effect = leaveout_effect)
)

pSB <- ggplot(lineage_plot, aes(y = mode, colour = mode)) +
  geom_vline(
    data = null_SB,
    aes(xintercept = null_value),
    inherit.aes = FALSE,
    linetype = "dashed", colour = "grey75", linewidth = 0.32
  ) +
  geom_segment(
    aes(x = full_effect, xend = leaveout_effect, yend = mode),
    linewidth = 0.70, alpha = 0.80, lineend = "round"
  ) +
  geom_point(
    data = point_SB,
    aes(x = effect, y = mode, shape = state, colour = mode),
    size = 2.15, stroke = 0.55
  ) +
  geom_segment(
    aes(x = leader_x0, xend = leader_x1, y = mode, yend = mode),
    inherit.aes = FALSE,
    colour = "grey55", linewidth = 0.28
  ) +
  geom_text(
    aes(x = label_x, y = mode, label = removed_lineage),
    inherit.aes = FALSE,
    hjust = 0,
    family = FIG_FONT, fontface = "italic",
    size = TEXT_SMALL / .pt, colour = "grey20", lineheight = 0.90
  ) +
  facet_wrap(~ outcome_short, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = MODE_COLS, guide = "none") +
  scale_shape_manual(values = c("Full model" = 1, "After lineage removal" = 16), name = NULL) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.08, 0.50))
  ) +
  labs(x = "Adjusted transmission mode association", y = NULL) +
  theme_publication() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = TEXT_SMALL, lineheight = 0.92),
    axis.text.y = element_text(size = TEXT_SMALL),
    plot.margin = margin(4, 24, 4, 4)
  )

# =============================================================================
# 10. Save supplementary panels and assemble
# =============================================================================
save_plot_all_formats(pSA, "FigS3_panelA_family_composition_FINAL_v22", SUPP_WIDTH_MM, SA_HEIGHT_MM, SUPP_DIR)
save_plot_all_formats(pSB, "FigS3_panelB_within_mode_lineage_influence_FINAL_v22", SUPP_WIDTH_MM, SB_HEIGHT_MM, SUPP_DIR)

supp_plot <- wrap_elements(full = pSA) / wrap_elements(full = pSB) +
  plot_layout(heights = c(SA_HEIGHT_MM, SB_HEIGHT_MM))
save_plot_all_formats(
  supp_plot, "FigS3_transmission_ecology_FINAL_v22",
  SUPP_WIDTH_MM, SUPP_HEIGHT_MM, SUPP_DIR
)
cat("Saved Fig.S3A-B and assembled supplementary figure.\n\n")

# =============================================================================
# 11. Final audit
# =============================================================================
cat("Fig.3A evidence summary (all values read from frozen result files):\n")
print(A_dat |> select(group, row, count))
cat("  S0-S2 union supporting sequences: ", scales::comma(s012_union_n), "\n\n", sep = "")
cat("Documented-vector associations used in Fig.3D:\n")
print(D_dat |> select(row_name, outcome_short, effect, conf_low, conf_high))
cat("\nOutputs written to:\n", OUT_DIR, "\n", sep = "")
cat("Main + individual panels: PDF / PNG / TIFF\n")
cat("Supplement + individual panels: PDF / PNG / TIFF\n")
cat("NEW Fig.3 publication plotting FINAL v22 completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
