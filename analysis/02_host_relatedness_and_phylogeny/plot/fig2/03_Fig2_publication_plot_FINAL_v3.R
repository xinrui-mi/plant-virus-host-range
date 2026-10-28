#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")) else normalizePath(getwd(), winslash = "/")

source(file.path(SCRIPT_DIR, "00_paths.R"))
source(file.path(SCRIPT_DIR, "00_check_packages.R"))
source(file.path(SCRIPT_DIR, "00_plot_theme.R"))
source(file.path(SCRIPT_DIR, "00_helpers.R"))
source(file.path(SCRIPT_DIR, "00_helpers_Fig2_v3.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ape)
  library(ggtree)
})

log_file <- file.path(OUT_LOG, "Fig2_publication_plot_FINAL_v14.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Fig.2 publication plotting FINAL v14 started: ", format(Sys.time()), "\n", sep = "")
cat("FIG3_ROOT: ", FIG3_ROOT, "\n", sep = "")
cat("Resolved figure font: ", FONT_FAMILY, "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")

f01 <- file.path(TABLE_DIR, "01_formal_taxid_host_sharing_by_distance.tsv")
f02 <- file.path(TABLE_DIR, "02_formal_taxid_global_decay_model.tsv")
f03 <- file.path(TABLE_DIR, "03_family_specific_decay_models.tsv")
f06 <- file.path(TABLE_DIR, "06_formal_taxid_degree_preserving_null_beta_iterations.tsv")
f12 <- file.path(FIG3_ROOT, "pipeline_results", "03L_FINAL_12family_host_spectrum", "07_FINAL_plot_ready_family_summary.tsv")
required <- c(f01, f02, f03, f06, f12)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) stop("Required frozen Fig.2 input(s) are missing:\n", paste(missing_required, collapse = "\n"))

input_manifest <- data.table(role = c("Panel A","Panel B observed","Panel B null","Panel C beta forest","Panel D rho dot plot","Supplementary Fig.2 atlas"), file = c(f01,f02,f06,f03,f12,f12), exists = TRUE)

# Panel A
A_dt <- read_tsv(f01)
pA <- make_distance_decay_plot(A_dt, line_colour = COL_BLUE, y_title = "Virus sharing probability (%)", show_pair_size = TRUE, legend_inside = TRUE) +
  theme(axis.text.x = element_text(angle = 23, hjust = 1, vjust = 1, face = "plain"))
save_triplet(pA, "Fig2A_host_taxonomic_distance_decay", OUT_PANEL, width_mm = 100, height_mm = 76)

# Panel B
obs_beta <- extract_observed_beta(f02)
B_dt <- read_tsv(f06)
bcol <- find_col(names(B_dt), exact = c("beta", "estimate", "coefficient", "null_beta"), regex = c("null.*beta", "^beta$", "estimate", "coef"))
B_dt[, .beta := as.numeric(get(bcol))]
B_dt <- B_dt[is.finite(.beta)]
if (!nrow(B_dt)) stop("Panel B null table contains no finite beta values: ", f06)
pB <- make_clean_hist_density_plot(B_dt, obs_beta)
save_triplet(pB, "Fig2B_degree_preserving_null", OUT_PANEL, width_mm = 84, height_mm = 76)

# Family-level data
c <- read_tsv(f03)
family_col <- find_col(names(c), exact = c("virus_family", "family"), regex = c("virus.*family", "^family$"))
beta_col <- find_col(names(c), exact = c("beta", "estimate", "coefficient"), regex = c("^beta$", "estimate", "coef"))
low_col <- find_col(names(c), exact = c("beta_ci_low", "conf_low", "conf.low", "lower_ci", "ci_lower", "beta_lower"), regex = c("beta.*ci.*low", "conf.*low", "lower.*ci", "ci.*lower", "q025"))
high_col <- find_col(names(c), exact = c("beta_ci_high", "conf_high", "conf.high", "upper_ci", "ci_upper", "beta_upper"), regex = c("beta.*ci.*high", "conf.*high", "upper.*ci", "ci.*upper", "q975"))
Cbeta <- data.table(virus_family = as.character(c[[family_col]]), beta = as.numeric(c[[beta_col]]), beta_low = as.numeric(c[[low_col]]), beta_high = as.numeric(c[[high_col]]))
Cbeta <- Cbeta[is.finite(beta)]
setorder(Cbeta, beta)
beta_family_order <- Cbeta$virus_family

d <- read_tsv(f12)
dfam <- find_col(names(d), exact = c("virus_family", "family"), regex = c("virus.*family", "^family$"))
rho_col <- find_col(names(d), exact = c("observed_rho", "rho"), regex = c("observed.*rho", "^rho$"))
sig_col <- find_signif_col(names(d))
qsup_col <- find_support_q_col(names(d))
Crho <- data.table(virus_family = as.character(d[[dfam]]), observed_rho = as.numeric(d[[rho_col]]))
if (!is.null(sig_col)) Crho[, support_sig := logicalize_signif(d[[sig_col]])] else if (!is.null(qsup_col)) Crho[, support_sig := as.numeric(d[[qsup_col]]) <= 0.05] else Crho[, support_sig := FALSE]
Crho <- Crho[virus_family %in% FINAL12_FAMILIES & is.finite(observed_rho)]
if (nrow(Crho) != 12L) stop("Final rho summary must contain 12 eligible families; got ", nrow(Crho))

extra_phylogeny_only <- setdiff(as.character(Crho$virus_family), beta_family_order)
family_order <- c(extra_phylogeny_only, beta_family_order)

missing_in_beta <- setdiff(family_order, as.character(Cbeta$virus_family))
if (length(missing_in_beta)) Cbeta <- rbind(Cbeta, data.table(virus_family = missing_in_beta, beta = NA_real_, beta_low = NA_real_, beta_high = NA_real_), fill = TRUE)
missing_in_rho <- setdiff(family_order, as.character(Crho$virus_family))
if (length(missing_in_rho)) Crho <- rbind(Crho, data.table(virus_family = missing_in_rho, observed_rho = NA_real_, support_sig = NA), fill = TRUE)

Cbeta[, in_rho_set := virus_family %in% as.character(Crho[is.finite(observed_rho), virus_family])]
Cbeta[, family_row := factor(virus_family, levels = rev(family_order))]
Crho[, family_row := factor(virus_family, levels = rev(family_order))]

cat("Family-set overlap audit:\n")
cat("  beta families: ", length(beta_family_order), "\n", sep = "")
cat("  rho families: ", nrow(Crho[is.finite(observed_rho)]), "\n", sep = "")
cat("  overlap: ", length(intersect(beta_family_order, as.character(Crho[is.finite(observed_rho), virus_family]))), "\n", sep = "")
cat("  union: ", length(family_order), "\n", sep = "")
cat("  rho-only family(ies): ", paste(extra_phylogeny_only, collapse = ";"), "\n", sep = "")

# Panel C: beta forest (24 families; shared y-axis with Panel D)
beta_xlim <- range(c(Cbeta$beta_low, Cbeta$beta_high, obs_beta), na.rm = TRUE)

pC <- ggplot(Cbeta, aes(y = family_row, x = beta)) +
  geom_hline(yintercept = seq_along(levels(Cbeta$family_row)),
             colour = "#EFEFEF", linewidth = 0.35) +
  geom_vline(xintercept = obs_beta, colour = COL_ORANGE,
             linewidth = 0.80, linetype = 2) +
  geom_segment(data = Cbeta[is.finite(beta_low) & is.finite(beta_high)],
               aes(x = beta_low, xend = beta_high, yend = family_row,
                   colour = in_rho_set),
               linewidth = 0.84, lineend = "round") +
  geom_point(data = Cbeta[is.finite(beta) & in_rho_set],
             shape = 21, size = 2.45, stroke = 0.30,
             colour = "black", fill = COL_BLUE) +
  geom_point(data = Cbeta[is.finite(beta) & !in_rho_set],
             shape = 21, size = 2.45, stroke = 0.30,
             colour = "black", fill = "white") +
  scale_colour_manual(values = c(`TRUE` = COL_BLUE,
                                  `FALSE` = "#7F7F7F"),
                      guide = "none") +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(limits = beta_xlim + c(-0.05, 0.03)) +
  labs(x = "Family-specific decay coefficient (β)", y = NULL) +
  PUB_THEME() +
  theme(axis.text.y = element_text(face = "italic", size = 6.2),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1.5, 1.0, 1.5, 3.0))


# Panel D: rho dot plot (same height and complete family rows as Panel C)
rho_xmax <- max(0.31, max(Crho$observed_rho, na.rm = TRUE) + 0.03)

# Keep the full family rows. Only draw points for families with rho values.
pD <- ggplot(Cbeta, aes(y = family_row, x = NA_real_)) +
  geom_hline(yintercept = seq_along(levels(Cbeta$family_row)),
             colour = "#EFEFEF", linewidth = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.50,
             colour = "#9E9E9E", linetype = 2) +
  geom_point(data = Crho[support_sig == TRUE & is.finite(observed_rho)],
             aes(x = observed_rho),
             shape = 21, size = 2.40, stroke = 0.42,
             colour = "black", fill = COL_BLUE) +
  geom_point(data = Crho[!is.na(support_sig) & !support_sig &
                           is.finite(observed_rho)],
             aes(x = observed_rho),
             shape = 21, size = 2.40, stroke = 0.42,
             colour = "black", fill = "white") +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(limits = c(-0.01, rho_xmax),
                     breaks = c(0, 0.10, 0.20, 0.30),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Spearman correlation (ρ)", y = NULL) +
  PUB_THEME() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1.5, 2.2, 1.5, 0.6))

# Embedded legend inside the rho panel
pD <- pD +
  annotate("point", x = 0.17, y = 24.5,
           colour = "black", fill = COL_BLUE,
           shape = 21, size = 2.2, stroke = 0.35) +
  annotate("point", x = 0.17, y = 22.8,
           colour = "black", fill = "white",
           shape = 21, size = 2.2, stroke = 0.35) +
  annotate("text", x = 0.22, y = 24.5,
           label = "Significant",
           hjust = 0, family = FONT_FAMILY, size = 2.4) +
  annotate("text", x = 0.22, y = 22.8,
           label = "Not significant",
           hjust = 0, family = FONT_FAMILY, size = 2.4)

cd_panel <- wrap_plots(wrap_elements(full = pC),
                       wrap_elements(full = pD),
                       ncol = 2,
                       widths = c(1.15, 0.85))

save_triplet(cd_panel,
             "Fig2CD_family_decay_and_phylogeny_v14",
             OUT_PANEL, width_mm = 180, height_mm = 136)

# Main figure
row_top <- wrap_plots(list(wrap_elements(full = pA),
                           wrap_elements(full = pB)),
                      ncol = 2,
                      widths = c(1.10, 0.90))

main_fig <- wrap_plots(list(wrap_elements(full = row_top),
                            wrap_elements(full = cd_panel)),
                       ncol = 1,
                       heights = c(0.34, 0.66))

save_triplet(main_fig,
             "Fig2_main_publication_FINAL_v14",
             OUT_MAIN, width_mm = 180, height_mm = 190)

FINAL_PANEL_DIR <- file.path(OUT_ROOT, "FINAL_individual_panels")
save_triplet(pC, "Fig2C_family_decay_beta_v14",
             FINAL_PANEL_DIR, 112, 136)
save_triplet(pD, "Fig2D_phylogenetic_association_dotplot_v14",
             FINAL_PANEL_DIR, 68, 136)
# Supplementary Fig.2 carried forward from v8
summary12 <- copy(Crho)
summary12[, virus_family_chr := as.character(virus_family)]
trait_list <- list(); trait_audit <- list(); tree_manifest <- list()
STAGE20_ASSOC <- resolve_stage20_association_master_v32(FIG3_ROOT)
HOST_ORDER_LOOKUP_OBJ <- read_stage20_host_taxid_order_lookup_v32(STAGE20_ASSOC)
HOST_ORDER_LOOKUP <- HOST_ORDER_LOOKUP_OBJ$lookup
if (nrow(HOST_ORDER_LOOKUP_OBJ$conflict)) fwrite(HOST_ORDER_LOOKUP_OBJ$conflict, file.path(OUT_LOG, "FigS2_FINAL_v14_host_taxid_order_conflicts.tsv"), sep = "\t")
for (fam in FINAL12_FAMILIES) {
  tf <- resolve_final12_tree_file(FIG3_ROOT, fam)
  tr <- ape::read.tree(tf)
  rr <- prepare_final12_family_traits_v32(FIG3_ROOT, fam, host_taxid_order_lookup = HOST_ORDER_LOOKUP, tree_file = tf)
  trait_list[[fam]] <- rr$trait
  trait_audit[[fam]] <- rr$audit
  tree_manifest[[fam]] <- data.table(virus_family = fam, tree_file = tf, tree_tip_n = length(tr$tip.label))
}
trait_audit_dt <- rbindlist(trait_audit, fill = TRUE)
tree_manifest_dt <- rbindlist(tree_manifest, fill = TRUE)
fwrite(trait_audit_dt, file.path(OUT_LOG, "FigS2_FINAL_v14_trait_source_audit.tsv"), sep = "\t")
fwrite(tree_manifest_dt, file.path(OUT_LOG, "FigS2_FINAL_v14_tree_manifest.tsv"), sep = "\t")
order_count <- rbindlist(lapply(trait_list, function(x) x[!is.na(host_order) & nzchar(host_order) & host_order != "Multiple host orders", .N, by = host_order]), fill = TRUE)[, .(N = sum(N)), by = host_order]
setorder(order_count, -N, host_order)
GLOBAL_TOP_ORDERS <- head(order_count$host_order, 8)
strong_cols <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2", "#EDC948", "#9C755F")
GLOBAL_ORDER_COLS <- setNames(strong_cols[seq_along(GLOBAL_TOP_ORDERS)], GLOBAL_TOP_ORDERS)
all_breadth <- unlist(lapply(trait_list, function(x) x$host_family_breadth[is.finite(x$host_family_breadth) & x$host_family_breadth > 0]))
if (!length(all_breadth)) stop("No host-family breadth values available for the final 12-family atlas.")
SUPP_BREADTH_MAX <- max(all_breadth, na.rm = TRUE)
SUPP_BREADTH_BREAKS <- unique(round(c(1, 5, 10, 20, SUPP_BREADTH_MAX)))
SUPP_BREADTH_BREAKS <- SUPP_BREADTH_BREAKS[SUPP_BREADTH_BREAKS > 0]
atlas_panels <- list()
for (fam in FINAL12_FAMILIES) {
  tf <- resolve_final12_tree_file(FIG3_ROOT, fam)
  tr <- ape::read.tree(tf)
  tr <- midpoint_for_display(tr)
  tr <- ape::ladderize(tr, right = FALSE)
  trait <- copy(trait_list[[fam]])
  pp <- make_endpoint_phylogram_v210(family = fam, tr = tr, trait = trait, top_host_orders = GLOBAL_TOP_ORDERS, order_palette = GLOBAL_ORDER_COLS, show_legend = FALSE, folded = FALSE, breadth_breaks = SUPP_BREADTH_BREAKS, breadth_limit_max = SUPP_BREADTH_MAX, inside_legend = FALSE, add_scale_bar = FALSE, scale_y_offset = 7, legend_position = "bottom") +
        labs(title = bquote(atop(italic(.(fam)), "n = " ~ .(length(tr$tip.label))))) +
    theme(plot.title = element_text(family = FONT_FAMILY, size = TEXT_SMALL, hjust = 0.5, lineheight = 0.95, margin = margin(b = 1.8)),
          plot.margin = margin(3.2, 4.2, 3.2, 4.2))
  atlas_panels[[fam]] <- pp
  save_triplet(pp, paste0("FigS2_", fam, "_full_phylogeny_v14"), OUT_SUPP, width_mm = 180, height_mm = 142)
}


pal_full <- c(GLOBAL_ORDER_COLS, Other = "#B8B8B8")
legend_df <- data.table(x = seq_along(c(GLOBAL_TOP_ORDERS, "Other")), y = 1, host_order_display = factor(c(GLOBAL_TOP_ORDERS, "Other"), levels = c(GLOBAL_TOP_ORDERS, "Other")), host_family_breadth = c(SUPP_BREADTH_BREAKS[1:min(length(SUPP_BREADTH_BREAKS), length(c(GLOBAL_TOP_ORDERS, "Other")))], rep(max(SUPP_BREADTH_BREAKS), max(0, length(c(GLOBAL_TOP_ORDERS, "Other")) - length(SUPP_BREADTH_BREAKS)))))
legend_plot <- ggplot(legend_df, aes(x = x, y = y, colour = host_order_display, size = log1p(host_family_breadth))) +
  geom_point(shape = 21, fill = "white", stroke = 0.38) +
  scale_colour_manual(values = pal_full, breaks = c(GLOBAL_TOP_ORDERS, "Other"), name = "Host order") +
  scale_size_continuous(range = c(0.40, 1.75), breaks = log1p(SUPP_BREADTH_BREAKS), labels = as.character(SUPP_BREADTH_BREAKS), name = "Host family breadth") +
  theme_void(base_family = FONT_FAMILY) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical", legend.title = element_text(size = TEXT_SMALL, face = "bold", colour = "black"), legend.text = element_text(size = TEXT_SMALL, face = "plain", colour = "black"), legend.key.height = grid::unit(3.0, "mm"), legend.key.width = grid::unit(2.8, "mm"), legend.box.spacing = grid::unit(0.5, "mm"), plot.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(order = 1, title.position = "top", nrow = 2, byrow = TRUE, override.aes = list(shape = 16, size = 3.2)), size = guide_legend(order = 2, title.position = "top", nrow = 1, override.aes = list(colour = "black", fill = "white", shape = 21, stroke = 0.38)))


atlas_3x4_grid <- wrap_plots(atlas_panels[FINAL12_FAMILIES], ncol = 4)
final_atlas <- wrap_plots(list(wrap_elements(full = atlas_3x4_grid), wrap_elements(full = legend_plot)), ncol = 1, heights = c(0.90, 0.10))
save_triplet(final_atlas, "FigS2_12family_phylogenetic_atlas_3x4_v14", OUT_SUPP, width_mm = 300, height_mm = 226)

# ========== 新增 4×3 ==========
atlas_4x3_grid <- wrap_plots(atlas_panels[FINAL12_FAMILIES], ncol = 3)
atlas_4x3 <- wrap_plots(list(wrap_elements(full = atlas_4x3_grid), wrap_elements(full = legend_plot)), ncol = 1, heights = c(0.90, 0.10))
save_triplet(atlas_4x3, "FigS2_12family_phylogenetic_atlas_4x3_v14", OUT_SUPP, width_mm = 260, height_mm = 260)
# ==============================

atlas_2x6_grid <- wrap_plots(atlas_panels[FINAL12_FAMILIES], ncol = 6)
atlas_2x6 <- wrap_plots(list(wrap_elements(full = atlas_2x6_grid), wrap_elements(full = legend_plot)), ncol = 1, heights = c(0.90, 0.10))
save_triplet(atlas_2x6, "FigS2_12family_phylogenetic_atlas_2x6_v14", OUT_SUPP, width_mm = 330, height_mm = 182)

fwrite(input_manifest, file.path(OUT_LOG, "Fig2_FINAL_v14_input_manifest.tsv"), sep = "\t")
writeLines(c(paste0("Resolved_font\t", FONT_FAMILY), paste0("R_version\t", R.version.string), paste0("FIG3_ROOT\t", FIG3_ROOT), paste0("Stage20_association_master_for_tree_colour_lookup\t", STAGE20_ASSOC), paste0("Global_top_host_orders\t", paste(GLOBAL_TOP_ORDERS, collapse = ";")), paste0("Shared_breadth_max\t", SUPP_BREADTH_MAX)), file.path(OUT_LOG, "Fig2_FINAL_v14_environment.txt"))
cat("Fig.2 publication plotting FINAL v14 completed: ", format(Sys.time()), "\n", sep = "")
cat("Main output: ", file.path(OUT_MAIN, "Fig2_main_publication_FINAL_v14.pdf"), "\n", sep = "")
cat("Supp 3x4: ", file.path(OUT_SUPP, "FigS2_12family_phylogenetic_atlas_3x4_v14.pdf"), "\n", sep = "")
cat("Supp 2x6: ", file.path(OUT_SUPP, "FigS2_12family_phylogenetic_atlas_2x6_v14.pdf"), "\n", sep = "")
