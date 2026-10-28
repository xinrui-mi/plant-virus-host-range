#!/usr/bin/env Rscript

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
source(file.path(SCRIPT_DIR, "00_display_rules.R"))
source(file.path(SCRIPT_DIR, "00_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

log_file <- file.path(OUT_LOG, "Fig4_publication_FINAL_v1_1.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Fig.4 publication plotting FINAL v1 started: ", format(Sys.time()), "\n", sep = "")
cat("NETWORK_ROOT: ", NETWORK_ROOT, "\n", sep = "")
cat("ASSOC_MASTER: ", ASSOC_MASTER, "\n", sep = "")
cat("Resolved figure font: ", FONT_FAMILY, "\n", sep = "")
cat("R version: ", R.version.string, "\n\n", sep = "")

# ============================================================
# Frozen inputs. Plotting only; no network reconstruction.
# ============================================================
f02 <- locate_table("02_network_edge_weights.tsv", "02_.*network.*edge.*weight.*\\.tsv$")
f03 <- locate_table("03_network_module_membership.tsv", "03_.*module.*membership.*\\.tsv$")
f04 <- locate_table("04_current_ICTV_genus_network_metrics.tsv", "04_.*genus.*network.*metric.*\\.tsv$")
f05 <- locate_table("05_candidate_thresholds.tsv", "05_.*candidate.*threshold.*\\.tsv$")
f08 <- locate_table("08_network_80pct_VC80_resampling_long.tsv", "08_.*80.*VC80.*resampling.*\\.tsv$", required = FALSE)
f11 <- locate_table("11_candidate_threshold_perturbation_long.tsv", "11_.*threshold.*perturbation.*\\.tsv$", required = FALSE)
f17 <- locate_table("17_projection_degree_preserving_randomization.tsv", "17_.*projection.*random.*\\.tsv$", required = FALSE)
f18 <- locate_table("18_projection_small_world_summary.tsv", "18_.*small.*world.*summary.*\\.tsv$", required = FALSE)
f19 <- locate_table("19_genus_pure_VC80_network_metrics_sensitivity.tsv", "19_.*genus.*pure.*VC80.*metric.*\\.tsv$", required = FALSE)
f_vector_models <- VECTOR_MODEL_TABLE

edges_raw <- read_tsv(f02)
modules_raw <- read_tsv(f03)
metrics_raw <- read_tsv(f04)
threshold_raw <- read_tsv(f05)
assoc_master <- read_tsv(ASSOC_MASTER)

# Frozen-partition guard (U26, 2026-10-27): this figure must render the PUBLISHED
# 12-module partition, not a re-run partition. `cluster_louvain()` is randomised,
# so a re-run of 05A cannot reproduce it; re-run output is written to tables/rerun/
# and is never read here. Verify the module-membership table is the frozen one.
.frozen_md5 <- "171e8d46ce9604ea7726958c78185dc6"
.this_md5 <- unname(tools::md5sum(f03))
if (!identical(.this_md5, .frozen_md5)) {
  stop("Module-membership table is not the frozen published partition ",
       "(MD5 ", .this_md5, " != frozen ", .frozen_md5, ").\n",
       "Fig. 4 only renders the frozen tables. If you re-ran the network ",
       "analysis, its output is in analysis/04_genus_host_network/tables/rerun/ ",
       "and must not be plotted as the published figure.")
}
message("Frozen published partition verified (MD5 match).")

thresholds <- extract_official_thresholds(threshold_raw)
metrics <- standardize_metrics(metrics_raw, thresholds)
edges <- standardize_edge_table(edges_raw)
modules <- standardize_modules(modules_raw)
tax_maps <- make_taxonomy_maps(assoc_master)

breadth_cutoff <- thresholds[threshold == "Host-family breadth", value]
bridge_cutoff <- thresholds[threshold == "Network bridge score", value]
support_cutoff <- thresholds[threshold == "VC80 support", value]

# Display-only landmarks derived from the frozen result values.
# They are not candidate thresholds and are never used to select genera.

# Fill formal family from the frozen three-source master only if needed.
metrics <- merge(metrics, tax_maps$genus_family, by = "virus_genus", all.x = TRUE, suffixes = c("", "_master"))
metrics[is.na(virus_family) | !nzchar(virus_family), virus_family := virus_family_master]
metrics[, virus_family_master := NULL]

edges <- merge(edges, tax_maps$genus_family, by = "virus_genus", all.x = TRUE)
edges <- merge(edges, tax_maps$host_order, by = "host_family", all.x = TRUE)

cat("Official thresholds extracted from result file:\n")
print(thresholds)
cat("\nVirus genera in metrics: ", nrow(metrics), "\n", sep = "")
cat("Official candidate genera: ", metrics[candidate == TRUE, uniqueN(virus_genus)], "\n", sep = "")
cat("Weighted genus-host-family edges: ", nrow(edges), "\n", sep = "")
cat("Network modules: ", uniqueN(modules$module_id), "\n", sep = "")
cat("Genus -> family mapping coverage: ", percent(mean(!is.na(edges$virus_family))), "\n", sep = "")
cat("Host family -> order mapping coverage: ", percent(mean(!is.na(edges$host_order))), "\n\n", sep = "")

fwrite(thresholds, file.path(OUT_LOG, "Fig4_FINAL_v1_1_extracted_official_thresholds.tsv"), sep = "\t")

# ============================================================
# A. Legacy design restored exactly in structure/colors:
#    Virus family x Host order aggregated overview.
#    Display-only Top-N/edge thinning are inherited from old Fig.5.
# ============================================================
a_raw <- edges[
  !is.na(virus_family) & nzchar(virus_family) &
  !is.na(host_order) & nzchar(host_order),
  .(weight = sum(edge_weight, na.rm = TRUE)),
  by = .(virus_family, host_order)
]
if (!nrow(a_raw)) stop("Panel A cannot be drawn: virus-family / host-order mapping is empty.")

keep_families_A <- a_raw[, .(total_weight = sum(weight)), by = virus_family][
  order(-total_weight, virus_family)
][seq_len(min(.N, TOP_FAMILY_N_A)), virus_family]

keep_orders_A <- a_raw[, .(total_weight = sum(weight)), by = host_order][
  order(-total_weight, host_order)
][seq_len(min(.N, TOP_HOST_ORDER_N_A)), host_order]

a_full <- a_raw[virus_family %in% keep_families_A & host_order %in% keep_orders_A][order(-weight)]
edge_cutoff_A <- if (nrow(a_full)) as.numeric(quantile(a_full$weight, A_EDGE_KEEP_QUANTILE, na.rm = TRUE)) else 0

a_plot <- a_full[weight >= edge_cutoff_A][order(-weight, virus_family, host_order)]
a_plot <- head(a_plot, MAX_EDGES_A)

virus_nodes_A <- a_plot[, .(total_weight = sum(weight)), by = virus_family][order(-total_weight, virus_family)]
# Wider display spacing than v4: this is visualization-only and does not alter the network.
virus_nodes_A[, `:=`(x = 0.93, y = rev(seq_len(.N)))]
host_nodes_A <- a_plot[, .(total_weight = sum(weight)), by = host_order][order(-total_weight, host_order)]
host_nodes_A[, `:=`(x = 2.70, y = rev(seq_len(.N)))]

edge_plot_A <- merge(a_plot, virus_nodes_A[, .(virus_family, x_from = x, y_from = y)], by = "virus_family")
edge_plot_A <- merge(edge_plot_A, host_nodes_A[, .(host_order, x_to = x, y_to = y)], by = "host_order")

pA <- ggplot() +
  geom_curve(
    data = edge_plot_A,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = weight, alpha = weight),
    curvature = 0.15, colour = "grey35"
  ) +
  geom_point(
    data = virus_nodes_A, aes(x = x, y = y, size = total_weight),
    shape = 21, fill = COL_A_VIRUS, colour = "black", stroke = 0.35
  ) +
  geom_point(
    data = host_nodes_A, aes(x = x, y = y, size = total_weight),
    shape = 21, fill = COL_A_HOST, colour = "black", stroke = 0.35
  ) +
  geom_text(
    data = virus_nodes_A, aes(x = 0.71, y = y, label = virus_family),
    hjust = 1, family = FONT_FAMILY, fontface = "italic",
    size = TEXT_SMALL / ggplot2::.pt
  ) +
  geom_text(
    data = host_nodes_A, aes(x = 2.91, y = y, label = host_order),
    hjust = 0, family = FONT_FAMILY, size = TEXT_SMALL / ggplot2::.pt
  ) +
  scale_size_continuous(
    name = "Node strength", range = c(2.10, 5.35),
    breaks = auto_three_legend_breaks(c(virus_nodes_A$total_weight, host_nodes_A$total_weight))
  ) +
  scale_linewidth_continuous(name = "Link support", range = c(0.22, 1.55), breaks = pretty_breaks(n = 3)) +
  scale_alpha_continuous(range = c(0.22, 0.78), guide = "none") +
  scale_x_continuous(
    limits = c(-0.58, 4.12),
    breaks = c(0.84, 2.78),
    labels = c("Virus families", "Host orders"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.07))
  ) +
  guides(
    linewidth = guide_legend(
      order = 1, override.aes = list(alpha = 0.72),
      keyheight = grid::unit(3.6, "mm"), keywidth = grid::unit(8.0, "mm")
    ),
    size = guide_legend(
      order = 2, override.aes = list(fill = "white", colour = "black", stroke = 0.25),
      keyheight = grid::unit(4.0, "mm"), keywidth = grid::unit(7.8, "mm"),
      label.position = "right", label.hjust = 0
    )
  ) +
  labs(x = NULL, y = NULL) +
  LEGACY_FIG5_THEME() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = TEXT_NORMAL, face = "bold"),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.38),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.12, "cm"),
    legend.key.height = grid::unit(0.27, "cm"),
    legend.key.width = grid::unit(0.46, "cm"),
    legend.title = element_text(size = TEXT_SMALL, face = "bold"),
    legend.text = element_text(size = TEXT_SMALL, vjust = 0.5),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0.8, 0),
    plot.margin = margin(4, 10, 4, 10)
  ) +
  coord_cartesian(clip = "off")

# ============================================================
# B. Legacy module-composition bubble design restored:
#    Virus family + Host order, top eight modules.
# ============================================================
module_summary <- modules[, .(node_n = uniqueN(node)), by = module_id][order(-node_n, module_id)]
module_summary[, `:=`(module_rank = seq_len(.N), module_display = paste0("M", seq_len(.N)))]
module_top_B <- head(module_summary$module_id, TOP_MODULE_N_B)
module_levels_B <- module_summary[module_id %in% module_top_B][order(module_rank), module_display]

modules_B <- merge(modules, module_summary[, .(module_id, module_display)], by = "module_id", all.x = TRUE)
modules_B <- modules_B[module_id %in% module_top_B]

virus_mod_B <- modules_B[node_type == "virus_genus", .(virus_genus = node, module_display)]
virus_map_B <- merge(virus_mod_B, tax_maps$genus_family, by = "virus_genus", all.x = TRUE)
b_virus_raw <- virus_map_B[!is.na(virus_family) & nzchar(virus_family), .(n = .N), by = .(module_display, virus_family)]

keep_fam_B <- b_virus_raw[, .(total_n = sum(n)), by = virus_family][order(-total_n, virus_family)]
keep_fam_B <- head(keep_fam_B$virus_family, TOP_VIRUS_FAMILY_N_B)
b_virus_all <- b_virus_raw[, component := fifelse(virus_family %in% keep_fam_B, virus_family, "Other")][
  , .(n = sum(n)), by = .(module_display, component)
][, prop := n / sum(n), by = module_display][, component_type := "Virus family"]
b_virus <- if (DROP_OTHER_IN_PANEL_B) b_virus_all[component != "Other"] else copy(b_virus_all)

host_mod_B <- modules_B[node_type == "host_family", .(host_family = node, module_display)]
host_map_B <- merge(host_mod_B, tax_maps$host_order, by = "host_family", all.x = TRUE)
b_host_raw <- host_map_B[!is.na(host_order) & nzchar(host_order), .(n = .N), by = .(module_display, host_order)]

keep_order_B <- b_host_raw[, .(total_n = sum(n)), by = host_order][order(-total_n, host_order)]
keep_order_B <- head(keep_order_B$host_order, TOP_HOST_ORDER_N_B)
b_host_all <- b_host_raw[, component := fifelse(host_order %in% keep_order_B, host_order, "Other")][
  , .(n = sum(n)), by = .(module_display, component)
][, prop := n / sum(n), by = module_display][, component_type := "Host order"]
b_host <- if (DROP_OTHER_IN_PANEL_B) b_host_all[component != "Other"] else copy(b_host_all)

b_comp <- rbindlist(list(b_virus, b_host), fill = TRUE)
b_comp[, module_display := factor(module_display, levels = module_levels_B)]
b_comp[, component_type := factor(component_type, levels = c("Virus family", "Host order"))]
component_order <- b_comp[, .(total = sum(n)), by = .(component_type, component)][order(component_type, total)]
# Use per-facet factor order through a stable combined factor.
b_comp[, component_plot := factor(component, levels = unique(component_order$component))]

virus_family_components_B <- unique(as.character(b_virus$component))
taxon_axis_labels_B <- function(x) {
  lab <- ifelse(
    x %in% virus_family_components_B,
    paste0("italic('", x, "')"),
    paste0("plain('", x, "')")
  )
  parse(text = lab)
}

pB <- ggplot(b_comp, aes(x = module_display, y = component_plot)) +
  geom_point(aes(size = n, fill = prop), shape = 21, colour = "grey25", stroke = 0.25) +
  facet_grid(component_type ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(name = "Number of nodes", range = c(0.68, 2.90), breaks = pretty_breaks(n = 3)) +
  guides(
    size = guide_legend(
      title.position = "top",
      keyheight = grid::unit(2.9, "mm"),
      keywidth = grid::unit(3.1, "mm")
    )
  ) +
  scale_fill_gradient(
    name = "Proportion\nwithin module",
    low = COL_B_LOW, high = COL_B_HIGH,
    labels = percent_format(accuracy = 1),
    guide = guide_colorbar(
      title.position = "top",
      barheight = grid::unit(9.0, "mm"),
      barwidth = grid::unit(1.8, "mm")
    )
  ) +
  scale_y_discrete(
    labels = taxon_axis_labels_B,
    expand = expansion(add = c(0.85, 1.15))
  ) +
  labs(x = "Network module", y = NULL) +
  LEGACY_FIG5_THEME() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = TEXT_SMALL),
    axis.text.y = element_text(size = TEXT_SMALL, face = "plain", margin = margin(r = 2)),
    legend.position = "right",
    legend.title = element_text(size = TEXT_SMALL, face = "bold", lineheight = 0.92),
    legend.text = element_text(size = TEXT_TINY),
    legend.key.height = grid::unit(0.10, "cm"),
    legend.key.width = grid::unit(0.15, "cm"),
    legend.spacing.y = grid::unit(0.020, "cm"),
    legend.spacing.x = grid::unit(0.015, "cm"),
    legend.box.spacing = grid::unit(2.8, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 0, 5.5, 5.5, unit = "pt"),
    axis.title.x = element_text(size = TEXT_NORMAL, face = "bold")
  )

# ============================================================
# C. Legacy breadth x bridge design restored.
#    Virus family = fill color; VC80 support = circle size.
# ============================================================
family_cols_C <- candidate_family_palette(metrics)
candidate_families <- setdiff(names(family_cols_C), "Other")
metrics[, virus_family_plot := fifelse(virus_family %in% candidate_families, virus_family, "Other")]
metrics[is.na(virus_family_plot) | !nzchar(virus_family_plot), virus_family_plot := "Other"]
metrics[, virus_family_plot := factor(virus_family_plot, levels = names(family_cols_C))]

label_tbl_C <- metrics[candidate == TRUE][order(-host_family_breadth, -bridge_score, -vc80_support)]
size_breaks_C <- auto_three_legend_breaks(label_tbl_C$vc80_support)
x_breaks_C <- unique(c(1, 20, breadth_cutoff, 50, max(metrics[candidate == TRUE, host_family_breadth], na.rm = TRUE)))

cat("Panel C VC80 legend breaks: ", paste(size_breaks_C, collapse = ", "), "\n", sep = "")
cat("Panel C focused breadth-axis breaks: ", paste(x_breaks_C, collapse = ", "), "\n", sep = "")

pC <- ggplot(metrics, aes(x = host_family_breadth, y = bridge_score)) +
  geom_vline(xintercept = breadth_cutoff, linetype = "dashed", colour = "grey45", linewidth = 0.40) +
  geom_hline(yintercept = bridge_cutoff, linetype = "dashed", colour = "grey45", linewidth = 0.40) +
  geom_point(
    data = metrics[candidate == FALSE],
    aes(size = vc80_support, fill = virus_family_plot),
    shape = 21, colour = "black", alpha = 0.28, stroke = 0.14
  ) +
  geom_point(
    data = metrics[candidate == TRUE],
    aes(size = vc80_support, fill = virus_family_plot),
    shape = 21, colour = "black", alpha = 0.96, stroke = 0.38
  ) +
  ggrepel::geom_text_repel(
    data = label_tbl_C, aes(label = virus_genus),
    family = FONT_FAMILY, fontface = "italic", size = TEXT_TINY / ggplot2::.pt,
    box.padding = 0.38, point.padding = 0.22,
    segment.size = 0.18, segment.color = "grey58",
    min.segment.length = 0, force = 2.6, force_pull = 0.55,
    max.time = 5, max.iter = 30000, seed = DISPLAY_SEED,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(
    values = family_cols_C,
    name = "Virus family",
    breaks = candidate_families,
    labels = function(x) parse(text = paste0("italic('", x, "')")),
    drop = FALSE
  ) +
  scale_size_continuous(
    name = "VC80 support", range = c(1.0, 4.4), trans = "sqrt",
    breaks = size_breaks_C, labels = label_comma()
  ) +
  scale_x_continuous(
    trans = focus_above_cutoff_trans(breadth_cutoff, lower_slope = 0.50),
    breaks = x_breaks_C,
    labels = label_comma(),
    expand = expansion(mult = c(0.025, 0.48))
  ) +
  guides(
    fill = guide_legend(
      order = 1, ncol = 1, byrow = TRUE,
      override.aes = list(size = 2.05, alpha = 1),
      keyheight = grid::unit(2.15, "mm"), keywidth = grid::unit(2.65, "mm"),
      title.position = "top", label.position = "right", label.hjust = 0
    ),
    size = guide_legend(
      order = 2, ncol = 1, byrow = TRUE,
      override.aes = list(fill = "white", alpha = 1, colour = "black", stroke = 0.25),
      keyheight = grid::unit(4.0, "mm"), keywidth = grid::unit(5.4, "mm"),
      title.position = "top", label.position = "right", label.hjust = 0
    )
  ) +
  labs(x = "Host family breadth", y = "Network bridge score") +
  coord_cartesian(clip = "off") +
  LEGACY_FIG5_THEME() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.spacing.y = grid::unit(0.10, "cm"),
    legend.spacing.x = grid::unit(0.02, "cm"),
    legend.key.height = grid::unit(0.16, "cm"),
    legend.key.width = grid::unit(0.26, "cm"),
    legend.title = element_text(size = TEXT_SMALL, face = "bold"),
    legend.text = element_text(size = TEXT_SMALL, margin = margin(l = 0.2)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 1.2, 2),
    axis.text.x = element_text(size = TEXT_TINY, margin = margin(t = 1))
  )

# Frozen prioritized-genus set reused by Supplementary Fig. S4 robustness panels.
# This is the exact candidate definition from the frozen network plotting package.
# It is NOT a new selection step and is not used to calculate Panel D.
cand <- metrics[candidate == TRUE][order(-host_family_breadth, -bridge_score, -vc80_support, virus_genus)]
if (!nrow(cand)) stop("No official candidate virus genera detected.")
if (uniqueN(cand$virus_genus) != 12L) {
  warning("Expected 12 official prioritized genera; found ", uniqueN(cand$virus_genus), ". Using frozen candidate flag as stored in the result table.")
}

# ============================================================
# D. Documented vector relationship: unadjusted vs evidence-adjusted
#    candidate effect. Plotting only: M0/M1/M3 are read from the frozen
#    candidate-contextualization result table; no logistic model is re-fit.
# ============================================================
vector_models <- read_tsv(f_vector_models)
required_vector_cols <- c(
  "model_id", "candidate_or", "candidate_or_ci_low",
  "candidate_or_ci_high", "candidate_p", "converged"
)
missing_vector_cols <- setdiff(required_vector_cols, names(vector_models))
if (length(missing_vector_cols)) {
  stop("Vector model table lacks required columns: ", paste(missing_vector_cols, collapse = ", "))
}

find_vector_model <- function(id) {
  z <- vector_models[model_id == id]
  if (nrow(z) != 1L) stop("Expected exactly one frozen row for model_id = ", id)
  if (!flag_from_any(z$converged)) warning("Frozen model did not report converged=TRUE: ", id)
  z
}

m0 <- find_vector_model("M0_unadjusted_all")
m1 <- find_vector_model("M1_VC80_adjusted_all")
m3 <- find_vector_model("M3_sequence_adjusted_all")

vector_D <- rbindlist(list(
  data.table(
    adjustment = "VC80 support", model = "Unadjusted",
    OR = as.numeric(m0$candidate_or), lo = as.numeric(m0$candidate_or_ci_low),
    hi = as.numeric(m0$candidate_or_ci_high), p = as.numeric(m0$candidate_p)
  ),
  data.table(
    adjustment = "VC80 support", model = "Adjusted",
    OR = as.numeric(m1$candidate_or), lo = as.numeric(m1$candidate_or_ci_low),
    hi = as.numeric(m1$candidate_or_ci_high), p = as.numeric(m1$candidate_p)
  ),
  data.table(
    adjustment = "Sequence count", model = "Unadjusted",
    OR = as.numeric(m0$candidate_or), lo = as.numeric(m0$candidate_or_ci_low),
    hi = as.numeric(m0$candidate_or_ci_high), p = as.numeric(m0$candidate_p)
  ),
  data.table(
    adjustment = "Sequence count", model = "Adjusted",
    OR = as.numeric(m3$candidate_or), lo = as.numeric(m3$candidate_or_ci_low),
    hi = as.numeric(m3$candidate_or_ci_high), p = as.numeric(m3$candidate_p)
  )
))
vector_D[, adjustment := factor(adjustment, levels = c("VC80 support", "Sequence count"))]
vector_D[, y := ifelse(model == "Unadjusted", 2, 1)]

if (any(!is.finite(vector_D$OR) | !is.finite(vector_D$lo) | !is.finite(vector_D$hi))) {
  stop("Non-finite OR/CI detected in frozen vector-adjustment model table.")
}

D_COL_UNADJUSTED <- "#2F8FB3"
D_COL_ADJUSTED <- "#E64B35"
d_x_min <- min(c(1, vector_D$lo), na.rm = TRUE)
d_x_max <- max(c(1, vector_D$hi), na.rm = TRUE)
d_x_limits <- c(max(0.2, d_x_min * 0.75), d_x_max * 1.20)

pD <- ggplot(vector_D, aes(x = OR, y = y, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.40, colour = "grey45") +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y), linewidth = 0.65) +
  geom_segment(aes(x = lo, xend = lo, y = y - 0.055, yend = y + 0.055), linewidth = 0.65) +
  geom_segment(aes(x = hi, xend = hi, y = y - 0.055, yend = y + 0.055), linewidth = 0.65) +
  geom_point(size = 3.0) +
  facet_wrap(~adjustment, nrow = 1) +
  scale_colour_manual(values = c("Unadjusted" = D_COL_UNADJUSTED, "Adjusted" = D_COL_ADJUSTED), guide = "none") +
  scale_y_continuous(
    limits = c(0.65, 2.35), breaks = c(1, 2), labels = c("Adjusted", "Unadjusted"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_x_log10(breaks = c(0.5, 1, 10, 100, 1000), labels = log_tick_labels) +
  coord_cartesian(xlim = d_x_limits, clip = "on") +
  labs(x = "Odds ratio for documented vector relationships", y = NULL) +
  PUB_THEME() +
  theme(
    strip.background = element_rect(fill = "grey94", colour = "grey65", linewidth = 0.25),
    strip.text = element_text(size = TEXT_SMALL, face = "bold"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = TEXT_SMALL),
    axis.title.x = element_text(size = TEXT_NORMAL, face = "bold"),
    panel.spacing.x = grid::unit(5.0, "mm"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
  )

fwrite(vector_D, file.path(OUT_LOG, "Fig4_FINAL_v1_1_documented_vector_relationship_models_plotted.tsv"), sep = "\t")

# ============================================================
# Main assembly: fixed physical panel geometry, matching stand-alone exports.
# Each plot is isolated with wrap_elements(full=...) exactly as in frozen Fig.3,
# preventing patchwork from re-aligning internal axes/legends in the composite.
# ============================================================
MAIN_WIDTH_MM <- 190
TOP_HEIGHT_MM <- 64
BOTTOM_HEIGHT_MM <- 68
MAIN_HEIGHT_MM <- TOP_HEIGHT_MM + BOTTOM_HEIGHT_MM
TOP_WIDTHS_MM <- c(100, 90)      # A, B
BOTTOM_WIDTHS_MM <- c(100, 90)   # C, D

cat(
  "Main/panel geometry (mm): main=", sprintf("%.2f x %.2f", MAIN_WIDTH_MM, MAIN_HEIGHT_MM),
  "; A/B=", sprintf("%.2f x %.2f", TOP_WIDTHS_MM[1], TOP_HEIGHT_MM), " / ", sprintf("%.2f x %.2f", TOP_WIDTHS_MM[2], TOP_HEIGHT_MM),
  "; C/D=", sprintf("%.2f x %.2f", BOTTOM_WIDTHS_MM[1], BOTTOM_HEIGHT_MM), " / ", sprintf("%.2f x %.2f", BOTTOM_WIDTHS_MM[2], BOTTOM_HEIGHT_MM), "\n",
  sep = ""
)

top_row <- patchwork::wrap_plots(
  list(
    patchwork::wrap_elements(full = pA),
    patchwork::wrap_elements(full = pB)
  ),
  ncol = 2,
  widths = TOP_WIDTHS_MM
)

bottom_row <- patchwork::wrap_plots(
  list(
    patchwork::wrap_elements(full = pC),
    patchwork::wrap_elements(full = pD)
  ),
  ncol = 2,
  widths = BOTTOM_WIDTHS_MM
)

main_fig <- patchwork::wrap_plots(
  list(
    patchwork::wrap_elements(full = top_row),
    patchwork::wrap_elements(full = bottom_row)
  ),
  ncol = 1,
  heights = c(TOP_HEIGHT_MM, BOTTOM_HEIGHT_MM)
)

save_triplet(
  main_fig, "Fig4_main_publication_FINAL_v1_1", OUT_MAIN,
  width_mm = MAIN_WIDTH_MM, height_mm = MAIN_HEIGHT_MM
)

# Stand-alone panels use exactly the same physical width/height as their main-figure cells.
save_triplet(pA, "Fig4A_network_overview", OUT_PANEL,
             width_mm = TOP_WIDTHS_MM[1], height_mm = TOP_HEIGHT_MM)
save_triplet(pB, "Fig4B_network_module_composition", OUT_PANEL,
             width_mm = TOP_WIDTHS_MM[2], height_mm = TOP_HEIGHT_MM)
save_triplet(pC, "Fig4C_host_breadth_bridge_candidates", OUT_PANEL,
             width_mm = BOTTOM_WIDTHS_MM[1], height_mm = BOTTOM_HEIGHT_MM)
save_triplet(pD, "Fig4D_documented_vector_relationship_adjustment", OUT_PANEL,
             width_mm = BOTTOM_WIDTHS_MM[2], height_mm = BOTTOM_HEIGHT_MM)

# ============================================================
# Supplementary robustness figure.
# Story: data robustness -> threshold robustness -> taxonomic sensitivity.
# The former small-world/null visualization is intentionally not included here;
# its frozen analysis tables remain untouched and available for Methods/tables.
# ============================================================
supp_plots <- list()

# S1: candidate recovery under repeated 80% VC80 subsampling.
# Formal virus genus names on the y axis are intentionally italicized (ICTV taxon typography).
if (!is.na(f08) && file.exists(f08)) {
  r8 <- read_tsv(f08)
  gcol <- find_col(names(r8), exact = c("virus_node", "virus_genus"), regex = c("virus.*node", "virus.*genus"), required = FALSE)
  icol <- find_col(names(r8), exact = c("iteration", "iter"), regex = c("iter"), required = FALSE)
  ccol <- find_col(names(r8), exact = c("candidate", "is_candidate"), regex = c("candidate"), required = FALSE)

  if (!is.null(gcol) && !is.null(icol) && !is.null(ccol)) {
    rr <- r8[get(gcol) %in% cand$virus_genus,
             .(recovery = mean(truthy(get(ccol)), na.rm = TRUE)), by = .(virus_genus = get(gcol))]
    rr <- merge(data.table(virus_genus = cand$virus_genus), rr, by = "virus_genus", all.x = TRUE)
    rr[is.na(recovery), recovery := 0]
    rr <- rr[order(recovery, virus_genus)]
    rr[, virus_genus := factor(virus_genus, levels = virus_genus)]

    pS1 <- ggplot(rr, aes(x = recovery, y = virus_genus)) +
      geom_segment(
        aes(x = 0, xend = recovery, yend = virus_genus),
        colour = "grey84", linewidth = 0.48
      ) +
      geom_point(
        shape = 21, size = 2.20,
        fill = COL_BLUE, colour = "white", stroke = 0.24
      ) +
      scale_x_continuous(
        labels = percent_format(accuracy = 1),
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        expand = expansion(mult = c(0, 0.01))
      ) +
      labs(x = "Recovery frequency of prioritized genera (%)", y = NULL) +
      PUB_THEME() +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = TEXT_SMALL, face = "italic"),
        axis.title.x = element_text(size = TEXT_NORMAL, face = "bold"),
        panel.grid.major.x = element_line(linewidth = 0.22, colour = "grey92"),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(3, 5, 3, 5)
      )

    supp_plots[["resampling"]] <- pS1
    save_triplet(pS1, "FigS4_candidate_recovery", OUT_SUPP, width_mm = 165, height_mm = 66)
  }
}

# S2: threshold robustness profile.
# This replaces the v26 heatmap. The three threshold dimensions are represented as:
# x = host-family-breadth percentile, facet = network-bridge-score percentile,
# colour = VC80 support threshold, y = candidate-set Jaccard similarity.
if (!is.na(f11) && file.exists(f11)) {
  r11 <- read_tsv(f11)
  gcol <- find_col(names(r11), exact = c("virus_node", "virus_genus"), regex = c("virus.*node", "virus.*genus"), required = FALSE)
  ccol <- find_col(names(r11), exact = c("candidate", "is_candidate"), regex = c("candidate"), required = FALSE)
  hcol <- find_col(names(r11), exact = c("host_cutoff", "host_breadth_cutoff"), regex = c("host.*cut", "breadth.*cut"), required = FALSE)
  bcol <- find_col(names(r11), exact = c("bridge_cutoff"), regex = c("bridge.*cut"), required = FALSE)
  scol <- find_col(names(r11), exact = c("support_cutoff", "vc80_support_cutoff"), regex = c("support.*cut"), required = FALSE)
  hqcol <- find_col(names(r11), exact = c("host_q"), regex = c("^host.*q$"), required = FALSE)
  bqcol <- find_col(names(r11), exact = c("bridge_q"), regex = c("^bridge.*q$"), required = FALSE)

  if (!is.null(gcol) && !is.null(ccol) && !is.null(hcol) && !is.null(bcol) &&
      !is.null(scol) && !is.null(hqcol) && !is.null(bqcol)) {

    main_set <- sort(cand$virus_genus)

    stab <- r11[, .(
      candidate_set = list(sort(unique(get(gcol)[truthy(get(ccol))])))
    ), by = .(
      host_q = suppressWarnings(as.numeric(get(hqcol))),
      bridge_q = suppressWarnings(as.numeric(get(bqcol))),
      breadth = suppressWarnings(as.numeric(get(hcol))),
      bridge = suppressWarnings(as.numeric(get(bcol))),
      support = suppressWarnings(as.numeric(get(scol)))
    )]

    stab[, jaccard := vapply(candidate_set, function(x) {
      u <- union(main_set, x)
      if (!length(u)) 1 else length(intersect(main_set, x)) / length(u)
    }, numeric(1))]

    ordinal_pct <- function(q) {
      z <- 100 * q
      ztxt <- ifelse(
        abs(z - round(z)) < 1e-8,
        sprintf("%.0f", z),
        sub("\\.0$", "", sprintf("%.1f", z))
      )
      suffix <- ifelse(
        abs(z - 90) < 1e-8, "th",
        ifelse(abs(z - 95) < 1e-8, "th",
        ifelse(abs(z - 97.5) < 1e-8, "th",
        ifelse(abs(z - 50) < 1e-8, "th",
        ifelse(abs(z - 75) < 1e-8, "th", "th"))))
      )
      paste0(ztxt, suffix)
    }

    host_levels <- sort(unique(stab$host_q))
    bridge_levels <- sort(unique(stab$bridge_q))
    support_levels <- sort(unique(stab$support))

    stab[, host_percentile := factor(
      host_q,
      levels = host_levels,
      labels = ordinal_pct(host_levels)
    )]
    stab[, bridge_percentile := factor(
      bridge_q,
      levels = bridge_levels,
      labels = paste0("Network bridge score\n", ordinal_pct(bridge_levels), " percentile")
    )]
    stab[, support_f := factor(
      support,
      levels = support_levels,
      labels = as.character(support_levels)
    )]

    # The formal main-analysis threshold combination is determined dynamically
    # from the frozen official cutoffs; it is marked with a thin black ring.
    official_point <- stab[
      abs(breadth - breadth_cutoff) < 1e-10 &
      abs(bridge - bridge_cutoff) < 1e-10 &
      abs(support - support_cutoff) < 1e-10
    ]

    support_cols_S2 <- setNames(
      c("#56B4E9", "#E69F00", "#009E73")[seq_along(support_levels)],
      as.character(support_levels)
    )

    support_shapes_S2 <- setNames(c(16, 17, 15)[seq_along(support_levels)], as.character(support_levels))

    pS2 <- ggplot(
      stab,
      aes(
        x = host_percentile,
        y = jaccard,
        colour = support_f,
        shape = support_f,
        group = support_f
      )
    ) +
      geom_point(
        position = position_dodge(width = 0.34),
        size = 2.05, alpha = 0.96
      ) +
      geom_point(
        data = official_point,
        aes(x = host_percentile, y = jaccard),
        inherit.aes = FALSE,
        position = position_dodge(width = 0.34),
        shape = 21, fill = NA, colour = "black",
        size = 3.20, stroke = 0.65
      ) +
      facet_wrap(~ bridge_percentile, nrow = 1) +
      scale_colour_manual(
        name = "VC80 support threshold",
        values = support_cols_S2,
        breaks = as.character(support_levels),
        drop = FALSE
      ) +
      scale_shape_manual(
        name = "VC80 support threshold",
        values = support_shapes_S2,
        breaks = as.character(support_levels),
        drop = FALSE
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        labels = percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.09))
      ) +
      labs(
        x = "Host family breadth percentile",
        y = "Jaccard similarity of prioritized genus sets"
      ) +
      PUB_THEME() +
      theme(
        panel.grid.major.y = element_line(linewidth = 0.22, colour = "grey91"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(
          fill = "grey94", colour = "grey65", linewidth = 0.25
        ),
        strip.text = element_text(size = TEXT_SMALL, face = "bold"),
        axis.text.x = element_text(size = TEXT_SMALL),
        axis.text.y = element_text(size = TEXT_SMALL),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = TEXT_SMALL, face = "bold"),
        legend.text = element_text(size = TEXT_SMALL),
        legend.key.width = grid::unit(7.0, "mm"),
        legend.key.height = grid::unit(2.8, "mm"),
        legend.spacing.x = grid::unit(0.8, "mm"),
        legend.margin = margin(0, 0, 1.5, 0),
        panel.spacing.x = grid::unit(3.5, "mm"),
        plot.margin = margin(4, 5, 5, 5)
      ) +
      guides(
        colour = guide_legend(order = 1, override.aes = list(size = 2.6)),
        shape = guide_legend(order = 1)
      )

    # Export the actual cutoff <-> percentile mapping for the figure legend.
    cutoff_map <- unique(stab[, .(
      host_q, breadth, bridge_q, bridge, support
    )])[order(host_q, bridge_q, support)]
    fwrite(
      cutoff_map,
      file.path(OUT_LOG, "Fig4_FINAL_v1_1_threshold_perturbation_cutoff_map.tsv"),
      sep = "\t"
    )

    supp_plots[["perturbation"]] <- pS2
    save_triplet(
      pS2,
      "FigS4_threshold_robustness",
      OUT_SUPP,
      width_mm = 165,
      height_mm = 68
    )
  }
}

# S3: sensitivity to restricting VC80s to a single current virus genus.
# Compare the three quantities that define the formal candidate rule.
if (!is.na(f19) && file.exists(f19)) {
  sensitivity_plot <- tryCatch({
    pure_raw <- read_tsv(f19)
    pure_metrics <- standardize_metrics(pure_raw, thresholds)

    cmp <- merge(
      metrics[, .(
        virus_genus,
        main_host_family_breadth = host_family_breadth,
        main_bridge_score = bridge_score,
        main_vc80_support = vc80_support,
        main_candidate = candidate
      )],
      pure_metrics[, .(
        virus_genus,
        sensitivity_host_family_breadth = host_family_breadth,
        sensitivity_bridge_score = bridge_score,
        sensitivity_vc80_support = vc80_support,
        sensitivity_candidate = candidate
      )],
      by = "virus_genus",
      all = FALSE
    )

    main_set <- sort(cmp[main_candidate == TRUE, virus_genus])
    sensitivity_set <- sort(cmp[sensitivity_candidate == TRUE, virus_genus])
    union_set <- union(main_set, sensitivity_set)
    jac <- if (!length(union_set)) {
      1
    } else {
      length(intersect(main_set, sensitivity_set)) / length(union_set)
    }

    sens_long <- rbindlist(list(
      cmp[, .(
        virus_genus,
        candidate = main_candidate,
        metric = "Host family breadth",
        main_value = main_host_family_breadth,
        sensitivity_value = sensitivity_host_family_breadth
      )],
      cmp[, .(
        virus_genus,
        candidate = main_candidate,
        metric = "Network bridge score",
        main_value = main_bridge_score,
        sensitivity_value = sensitivity_bridge_score
      )],
      cmp[, .(
        virus_genus,
        candidate = main_candidate,
        metric = "VC80 support",
        main_value = main_vc80_support,
        sensitivity_value = sensitivity_vc80_support
      )]
    ))

    sens_long <- sens_long[
      is.finite(main_value) & is.finite(sensitivity_value)
    ]
    sens_long[, candidate_status := fifelse(candidate, "Prioritized genera", "Other genera")]
    sens_long[, metric := factor(
      metric,
      levels = c("Host family breadth", "Network bridge score", "VC80 support")
    )]

    # Force matching x/y ranges within each facet so the y=x line is visually meaningful.
    sens_limits <- sens_long[, {
      lo <- min(c(main_value, sensitivity_value), na.rm = TRUE)
      hi <- max(c(main_value, sensitivity_value), na.rm = TRUE)
      pad <- if (is.finite(hi - lo) && hi > lo) 0.035 * (hi - lo) else 0.5
      .(
        x = c(lo - pad, hi + pad),
        y = c(lo - pad, hi + pad)
      )
    }, by = metric]

    p <- ggplot(sens_long, aes(x = main_value, y = sensitivity_value)) +
      geom_blank(data = sens_limits, aes(x = x, y = y), inherit.aes = FALSE) +
      geom_abline(
        intercept = 0, slope = 1,
        linewidth = 0.42, linetype = "dashed", colour = "grey58"
      ) +
      geom_point(
        data = sens_long[candidate == FALSE],
        aes(colour = "Other genera"),
        size = 0.62, alpha = 0.24
      ) +
      geom_point(
        data = sens_long[candidate == TRUE],
        aes(colour = "Prioritized genera"),
        size = 1.72, alpha = 0.98
      ) +
      scale_colour_manual(
        name = NULL,
        values = c(
          "Other genera" = "grey68",
          "Prioritized genera" = COL_BLUE
        ),
        breaks = c("Other genera", "Prioritized genera")
      ) +
      guides(
        colour = guide_legend(
          override.aes = list(
            size = c(1.15, 1.80),
            alpha = c(0.55, 1.00)
          )
        )
      ) +
      facet_wrap(~ metric, nrow = 1, scales = "free") +
      labs(
        x = "Main analysis",
        y = "Genus-pure VC80 analysis"
      ) +
      PUB_THEME() +
      theme(
        aspect.ratio = 1,
        panel.grid.major = element_line(linewidth = 0.20, colour = "grey92"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(
          fill = "grey94", colour = "grey65", linewidth = 0.25
        ),
        strip.text = element_text(size = TEXT_SMALL, face = "bold"),
        axis.text = element_text(size = TEXT_TINY),
        axis.title = element_text(size = TEXT_NORMAL, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = TEXT_SMALL),
        legend.key.width = grid::unit(3.5, "mm"),
        legend.key.height = grid::unit(2.6, "mm"),
        legend.spacing.x = grid::unit(1.0, "mm"),
        legend.margin = margin(2.0, 0, 0, 0),
        panel.spacing.x = grid::unit(3.5, "mm"),
        plot.margin = margin(7, 4, 3, 4)
      )

    fwrite(
      data.table(
        main_candidate_n = length(main_set),
        single_genus_candidate_n = length(sensitivity_set),
        shared_candidate_n = length(intersect(main_set, sensitivity_set)),
        candidate_set_jaccard = jac
      ),
      file.path(OUT_LOG, "Fig4_FINAL_v1_1_single_genus_VC80_candidate_sensitivity.tsv"),
      sep = "\t"
    )

    p
  }, error = function(e) {
    warning("Single genus VC80 sensitivity panel skipped: ", conditionMessage(e))
    NULL
  })

  if (!is.null(sensitivity_plot)) {
    supp_plots[["single_genus"]] <- sensitivity_plot
    save_triplet(
      sensitivity_plot,
      "FigS4_single_genus_VC80_sensitivity",
      OUT_SUPP,
      width_mm = 165,
      height_mm = 70
    )
  }
}

if (length(supp_plots)) {
  height_map <- c(
    resampling = 0.84,
    perturbation = 0.98,
    single_genus = 1.08
  )
  supp_heights <- unname(height_map[names(supp_plots)])
  supp_heights[!is.finite(supp_heights)] <- 0.90

  supp_fig <- wrap_plots(
    supp_plots,
    ncol = 1,
    heights = supp_heights
  )

  supp_height <- min(225, 52 + 54 * sum(supp_heights))
  save_triplet(
    supp_fig,
    "FigS4_network_robustness_FINAL_v1_1",
    OUT_SUPP,
    width_mm = 180,
    height_mm = supp_height
  )
}

# Reproducibility logs.
fwrite(metrics[candidate == TRUE][order(-host_family_breadth, -bridge_score, -vc80_support)],
       file.path(OUT_LOG, "Fig4_FINAL_v1_1_plotted_candidate_genera.tsv"), sep = "\t")
fwrite(data.table(
  display_rule = c(
    "Panel A top virus families", "Panel A top host orders", "Panel A edge keep quantile", "Panel A maximum displayed edges",
    "Panel B displayed modules", "Panel B top virus families", "Panel B top host orders", "Panel B drop Other"
  ),
  value = c(
    TOP_FAMILY_N_A, TOP_HOST_ORDER_N_A, A_EDGE_KEEP_QUANTILE, MAX_EDGES_A,
    TOP_MODULE_N_B, TOP_VIRUS_FAMILY_N_B, TOP_HOST_ORDER_N_B, DROP_OTHER_IN_PANEL_B
  )
), file.path(OUT_LOG, "Fig4_FINAL_v1_1_display_only_rules.tsv"), sep = "\t")

cat("\nOutputs written to: ", OUT_ROOT, "\n", sep = "")
cat("Main figure: PDF + PNG + TIFF\n")
cat("Individual panels: PDF + PNG + TIFF\n")
cat("Supplementary robustness: PDF + PNG + TIFF\n")
cat("Fig.4 publication plotting FINAL v1 completed: ", format(Sys.time()), "\n", sep = "")
