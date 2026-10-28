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
source(file.path(SCRIPT_DIR, "00_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ape)
  library(ggtree)
})

log_file <- file.path(OUT_LOG, "Fig3_publication_plot_v2_11_1.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Fig.3 publication plotting v2.11.1 started: ", format(Sys.time()), "\n", sep = "")
cat("FIG3_ROOT: ", FIG3_ROOT, "\n", sep = "")
cat("Resolved figure font: ", FONT_FAMILY, "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")

# ============================================================
# Frozen formal inputs
# ============================================================
f01 <- file.path(TABLE_DIR, "01_formal_taxid_host_sharing_by_distance.tsv")
f02 <- file.path(TABLE_DIR, "02_formal_taxid_global_decay_model.tsv")
f04 <- file.path(TABLE_DIR, "04_RdRP65_host_sharing_by_distance.tsv")
f06 <- file.path(TABLE_DIR, "06_formal_taxid_degree_preserving_null_beta_iterations.tsv")
f26 <- file.path(TABLE_DIR, "26_03G_pair_sampling_stability_summary.tsv")
f27 <- file.path(TABLE_DIR, "27_03G_support_stratified_permutation_summary.tsv")

required <- c(f01, f02, f06, f26, f27)
missing_required <- required[!file.exists(required)]

if (length(missing_required)) {
  stop(
    "Required frozen Fig.3 input(s) are missing:\n",
    paste(missing_required, collapse = "\n")
  )
}

input_manifest <- data.table(
  role = c(
    "Panel A",
    "Panel B observed",
    "Panel B null",
    "Panel D pair-sampling",
    "Panel D support-stratified permutation",
    "Supplementary RdRP65"
  ),
  file = c(f01, f02, f06, f26, f27, f04),
  exists = file.exists(c(f01, f02, f06, f26, f27, f04))
)

# ============================================================
# Panel A
# Global formal-virus-taxid sharing decay
# ============================================================
a <- read_tsv(f01)

cat("Panel A columns: ", paste(names(a), collapse = ", "), "\n", sep = "")

pA <- make_distance_decay_plot(
  a,
  line_colour = COL_BLUE,
  y_title = "Virus-sharing probability (%)",
  show_pair_size = TRUE
)

save_triplet(
  pA,
  "Fig3A_formal_taxid_host_taxonomic_decay",
  OUT_PANEL,
  width_mm = 98,
  height_mm = 72
)

# ============================================================
# Panel B
# Degree-preserving null model
# ============================================================
obs_beta <- extract_observed_beta(f02)
formal_p <- extract_empirical_p(f02)

b <- read_tsv(f06)
bcol <- find_col(
  names(b),
  exact = c("beta", "estimate", "coefficient", "null_beta"),
  regex = c("null.*beta", "^beta$", "estimate", "coef")
)
b[, .beta := as.numeric(get(bcol))]
b <- b[is.finite(.beta)]

if (!nrow(b)) stop("Panel B null table contains no finite beta values: ", f06)

if (!is.finite(formal_p)) {
  formal_p <- (sum(abs(b$.beta) >= abs(obs_beta)) + 1) / (nrow(b) + 1)
}

xrange <- range(c(b$.beta, obs_beta), finite = TRUE)
xspan <- diff(xrange)

dens_obj <- count_scaled_density_explicit(
  b$.beta,
  bins = 30
)

pB <- ggplot(
  b,
  aes(x = .beta)
)

if (is.finite(dens_obj$binwidth)) {
  pB <- pB +
    geom_histogram(
      binwidth = dens_obj$binwidth,
      boundary = dens_obj$boundary,
      closed = "left",
      fill = COL_GREY,
      colour = "white",
      linewidth = 0.24
    )
} else {
  pB <- pB +
    geom_histogram(
      bins = 30,
      fill = COL_GREY,
      colour = "white",
      linewidth = 0.24
    )
}

if (nrow(dens_obj$df)) {
  pB <- pB +
    geom_line(
      data = dens_obj$df,
      aes(
        x = x,
        y = y_count_scaled
      ),
      inherit.aes = FALSE,
      linewidth = 0.72,
      colour = "#666666"
    )
}

pB <- pB +
  geom_vline(
    xintercept = obs_beta,
    linewidth = 0.90,
    colour = COL_ORANGE
  ) +
  annotate(
    "text",
    x = obs_beta + 0.018 * xspan,
    y = Inf,
    label = paste0(
      "Observed β = ",
      formatC(
        obs_beta,
        digits = 3,
        format = "f"
      )
    ),
    hjust = 0,
    vjust = 1.35,
    family = FONT_FAMILY,
    size = 2.45,
    colour = COL_DARK
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.09))
  ) +
  labs(
    x = "Logistic decay coefficient (β)",
    y = "Permutation count"
  ) +
  PUB_THEME()

save_triplet(
  pB,
  "Fig3B_degree_preserving_null",
  OUT_PANEL,
  width_mm = 82,
  height_mm = 72
)

# ============================================================
# Panel C
# Representative Geminiviridae branch-length phylogram.
# Display compression remains topology + clade-size only.
# Formal 03F/03G statistics use the complete original tree.
# ============================================================
tree_manifest <- data.table(
  virus_family = c(
    "Geminiviridae",
    "Caulimoviridae",
    "Potyviridae",
    "Tombusviridae"
  ),
  tree_file =
    NA_character_,
  tree_tip_n =
    NA_integer_
)

for (
  i in seq_len(
    nrow(tree_manifest)
  )
) {
  fam <-
    tree_manifest$virus_family[i]

  tf <- find_tree_for_family(
    FIG3_ROOT,
    fam
  )

  tree_manifest[
    i,
    tree_file := tf
  ]

  if (
    !is.na(tf) &&
      file.exists(tf)
  ) {
    tr0 <- ape::read.tree(
      tf
    )

    tree_manifest[
      i,
      tree_tip_n :=
        length(
          tr0$tip.label
        )
    ]
  }
}

host_long_file <-
  find_tip_host_family_file(
    FIG3_ROOT
  )

host_long <-
  read_tip_host_family_long(
    host_long_file
  )

ann_master_file <-
  find_tip_annotation_master(
    FIG3_ROOT
  )

ann_master <-
  read_tip_annotation_master(
    ann_master_file
  )

cat(
  "Tip × host-family annotation file: ",
  host_long_file,
  "\n",
  sep = ""
)

cat(
  "Tip × host-family rows recovered: ",
  nrow(host_long),
  "\n",
  sep = ""
)

cat(
  "Tree-tip annotation master: ",
  ann_master_file,
  "\n",
  sep = ""
)

cat(
  "Tree-tip annotation columns: ",
  paste(
    names(
      ann_master
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

input_manifest <- rbind(
  input_manifest,
  data.table(
    role = c(
      "Panel C tip × host-family long",
      "Panel C tree-tip annotation master"
    ),
    file = c(
      host_long_file,
      ann_master_file
    ),
    exists = c(
      !is.na(
        host_long_file
      ) &&
        file.exists(
          host_long_file
        ),
      !is.na(
        ann_master_file
      ) &&
        file.exists(
          ann_master_file
        )
    )
  ),
  fill = TRUE
)

get_tree_file <- function(
    fam
) {
  z <- tree_manifest[
    virus_family == fam,
    tree_file
  ]

  if (
    !length(z) ||
      is.na(z)
  ) {
    stop(
      "Formal tree not detected for ",
      fam
    )
  }

  z
}

# Build traits for ALL four formal families here.
# This allows one shared host-order colour dictionary while preserving
# family-specific Top-order legends in Supplementary.
all_family_traits <- lapply(
  tree_manifest$virus_family,
  function(fam) {
    prepare_family_tip_traits(
      fam,
      get_tree_file(
        fam
      ),
      ann_master,
      host_long
    )
  }
)

names(
  all_family_traits
) <- tree_manifest$virus_family

trait_gem_raw <-
  all_family_traits[["Geminiviridae"]]

TOP6_ORDERS <-
  family_top_orders_v210(
    trait_gem_raw,
    top_n = 6
  )

# Family-specific Top 8 sets used later in Supplementary.
FAMILY_TOP8 <- lapply(
  all_family_traits,
  family_top_orders_v210,
  top_n = 8
)

# Build a deterministic global dictionary:
# main-panel Top 6 first, then additional Supplementary orders.
GLOBAL_ORDER_LEVELS <- unique(
  c(
    TOP6_ORDERS,
    unlist(
      FAMILY_TOP8,
      use.names = FALSE
    )
  )
)

GLOBAL_ORDER_COLS <-
  shared_order_palette_v210(
    GLOBAL_ORDER_LEVELS
  )

ORDER_COLS_MAIN <-
  GLOBAL_ORDER_COLS[
    TOP6_ORDERS
  ]

cat(
  "Panel C Top 6 dominant host orders: ",
  paste(
    TOP6_ORDERS,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

MAIN_BREADTH_MAX <-
  max(
    trait_gem_raw$host_family_breadth,
    na.rm = TRUE
  )

# Publication-facing representative legend values.
# The underlying point-size mapping remains continuous over the full observed
# range (1–17); the legend uses evenly interpretable landmarks.
MAIN_BREADTH_BREAKS <- c(
  1,
  5,
  10,
  15
)

cat(
  "Panel C host-family breadth breaks: ",
  paste(
    MAIN_BREADTH_BREAKS,
    collapse = ", "
  ),
  "; max = ",
  MAIN_BREADTH_MAX,
  "\n",
  sep = ""
)

tr_gem0 <- ape::read.tree(
  get_tree_file(
    "Geminiviridae"
  )
)

tr_gem_disp <-
  midpoint_for_display(
    tr_gem0
  )

fold_gem <-
  make_topology_folded_display_v27(
    tr_gem_disp,
    trait_gem_raw,
    target_display_n = 330,
    min_clade_n = 4,
    max_clade_fraction = 0.30
  )

cat(
  "Panel C topology folds selected: Geminiviridae=",
  nrow(
    fold_gem$selected
  ),
  "\n",
  sep = ""
)

cat(
  "Panel C display compression: Geminiviridae ",
  length(
    tr_gem0$tip.label
  ),
  " -> ",
  fold_gem$display_tip_n,
  " display units\n",
  sep = ""
)

if (
  nrow(
    fold_gem$selected
  )
) {
  fold_summary <-
    fold_gem$selected[
      ,
      .(
        fold_node =
          node,
        collapsed_tip_n =
          clade_n
      )
    ]

  fwrite(
    fold_summary,
    file.path(
      OUT_LOG,
      "Fig3C_v2_11_1_folded_clades.tsv"
    ),
    sep = "\t"
  )
}

pC <-
  make_endpoint_phylogram_v210(
    family =
      "Geminiviridae",
    tr =
      fold_gem$tree,
    trait =
      fold_gem$trait,
    top_host_orders =
      TOP6_ORDERS,
    order_palette =
      ORDER_COLS_MAIN,
    show_legend =
      TRUE,
    folded =
      TRUE,
    breadth_breaks =
      MAIN_BREADTH_BREAKS,
    breadth_limit_max =
      MAIN_BREADTH_MAX,
    inside_legend =
      TRUE,
    add_scale_bar =
      TRUE,
    scale_y_offset =
      10
  ) +
  theme(
    # Reserve a clean left gutter for panel label C.
    plot.title = element_text(
      family = FONT_FAMILY,
      size = TEXT_SMALL,
      face = "italic",
      colour = "black",
      hjust = 0,
      margin = margin(
        b = 2,
        l = 34
      )
    )
  )

save_triplet(
  pC,
  "Fig3C_Geminiviridae_phylogram_endpoint_traits",
  OUT_PANEL,
  width_mm = 172,
  height_mm = 132
)

# ============================================================
# Panel D
# Full-width short summary strip.
# Point = observed Spearman rho.
# Horizontal line = pair-sampling 95% interval.
# P values remain in the figure legend.
# ============================================================
d26 <- read_tsv(f26)
d27 <- read_tsv(f27)

need26 <- c(
  "virus_family",
  "rho_median",
  "rho_q025",
  "rho_q975"
)

need27 <- c(
  "virus_family",
  "observed_rho",
  "empirical_p_two_sided_abs"
)

if (
  !all(
    need26 %in% names(d26)
  )
) {
  stop(
    "Unexpected 26_03G schema."
  )
}

if (
  !all(
    need27 %in% names(d27)
  )
) {
  stop(
    "Unexpected 27_03G schema."
  )
}

families_order <- c(
  "Geminiviridae",
  "Caulimoviridae",
  "Potyviridae",
  "Tombusviridae"
)

d <- merge(
  d27[
    ,
    .(
      virus_family,
      observed_rho,
      empirical_p_two_sided_abs
    )
  ],
  d26[
    ,
    .(
      virus_family,
      rho_median,
      rho_q025,
      rho_q975
    )
  ],
  by = "virus_family",
  all.x = TRUE
)

d <- d[
  virus_family %in%
    families_order
]

d[, y_pos := c(
  Geminiviridae = 4,
  Caulimoviridae = 3,
  Potyviridae = 2,
  Tombusviridae = 1
)[virus_family]]

pD <- ggplot(d) +
  geom_segment(
    aes(
      x = rho_q025,
      xend = rho_q975,
      y = y_pos,
      yend = y_pos,
      colour = virus_family
    ),
    linewidth = 0.95,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      x = observed_rho,
      y = y_pos,
      colour = virus_family
    ),
    size = 2.60
  ) +
  scale_colour_manual(
    values = FAMILY_COLORS,
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = 1:4,
    labels = c(
      "Tombusviridae",
      "Potyviridae",
      "Caulimoviridae",
      "Geminiviridae"
    ),
    limits = c(
      0.55,
      4.45
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_x_continuous(
    limits = c(
      0,
      0.26
    ),
    breaks = seq(
      0,
      0.25,
      by = 0.05
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  labs(
    x = "Phylogeny–host spectrum correlation (ρ)",
    y = NULL
  ) +
  PUB_THEME() +
  theme(
    axis.text.y = element_text(
      face = "italic",
      size = TEXT_SMALL
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    axis.line.x = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    plot.margin = margin(
      2,
      6,
      2,
      8
    )
  )

save_triplet(
  pD,
  "Fig3D_phylogeny_host_spectrum_correlation",
  OUT_PANEL,
  width_mm = 172,
  height_mm = 48
)

# ============================================================
# MAIN FIGURE ASSEMBLY
# A/B/C/D are overlaid in reserved panel margins.
# ============================================================
A_tag <-
  tag_plot_inset_v211(
    pA,
    "A",
    left = 0.006,
    bottom = 0.925,
    right = 0.060,
    top = 0.995
  )

B_tag <-
  tag_plot_inset_v211(
    pB,
    "B",
    left = 0.006,
    bottom = 0.925,
    right = 0.060,
    top = 0.995
  )

C_tag <-
  tag_plot_inset_v211(
    pC,
    "C",
    left = 0.002,
    bottom = 0.932,
    right = 0.070,
    top = 0.997
  )

D_tag <-
  tag_plot_inset_v211(
    pD,
    "D",
    left = 0.004,
    bottom = 0.855,
    right = 0.075,
    top = 0.990
  )

top_row <-
  patchwork::wrap_plots(
    list(
      patchwork::wrap_elements(
        full = A_tag
      ),
      patchwork::wrap_elements(
        full = B_tag
      )
    ),
    ncol = 2,
    widths = c(
      1.10,
      0.90
    )
  )

main_fig <-
  patchwork::wrap_plots(
    list(
      patchwork::wrap_elements(
        full = top_row
      ),
      patchwork::wrap_elements(
        full = C_tag
      ),
      patchwork::wrap_elements(
        full = D_tag
      )
    ),
    ncol = 1,
    heights = c(
      0.34,
      0.50,
      0.16
    )
  )

save_triplet(
  main_fig,
  "Fig3_main_publication_v2_11_1",
  OUT_MAIN,
  width_mm = 180,
  height_mm = 222
)

# ============================================================
# SUPPLEMENTARY A
# RdRP65 host-sharing decay
# ============================================================
if (file.exists(f04)) {
  s <- read_tsv(f04)
  cat("Supplementary RdRP65 columns: ", paste(names(s), collapse = ", "), "\n", sep = "")

  pS_A <- make_distance_decay_plot(
    s,
    line_colour = COL_ORANGE,
    y_title = "Probability of sharing RdRP65 units (%)",
    show_pair_size = TRUE
  )

  save_triplet(
    pS_A,
    "FigS3A_RdRP65_host_taxonomic_decay",
    OUT_SUPP,
    width_mm = 98,
    height_mm = 75
  )
} else {
  warning("RdRP65 frozen table not found; supplementary RdRP65 panel skipped: ", f04)
}

# ============================================================
# SUPPLEMENTARY B-E
# Complete uncompressed branch-length phylograms.
# Each family has its own Top 8 host orders + Other,
# while recurring orders retain a shared global colour.
# ============================================================
supp_families <- c(
  "Geminiviridae",
  "Caulimoviridae",
  "Potyviridae",
  "Tombusviridae"
)

supp_letters <- c(
  "B",
  "C",
  "D",
  "E"
)

SUPP_BREADTH_MAX <-
  max(
    unlist(
      lapply(
        all_family_traits,
        function(x) {
          x$host_family_breadth[
            is.finite(
              x$host_family_breadth
            )
          ]
        }
      )
    ),
    na.rm = TRUE
  )

cat(
  "Supplementary shared host-family breadth maximum: ",
  SUPP_BREADTH_MAX,
  "\n",
  sep = ""
)

supp_tree_plots <-
  list()

long_branch_qc <-
  list()

for (
  i in seq_along(
    supp_families
  )
) {
  fam <-
    supp_families[i]

  tr_raw <-
    ape::read.tree(
      get_tree_file(
        fam
      )
    )

  tr_full <-
    midpoint_for_display(
      tr_raw
    )

  tr_full <-
    ape::ladderize(
      tr_full,
      right = FALSE
    )

  trait_full <-
    copy(
      all_family_traits[[fam]]
    )

  trait_full[, collapsed := FALSE]
  trait_full[, collapsed_n := 1L]

  family_orders <-
    FAMILY_TOP8[[fam]]

  family_cols <-
    GLOBAL_ORDER_COLS[
      family_orders
    ]

  family_max <-
    max(
      trait_full$host_family_breadth,
      na.rm = TRUE
    )

  # Regular, human-readable legend landmarks.
  # Examples with current data:
  # Geminiviridae (max 17): 1, 5, 10, 15
  # Caulimoviridae (max 5): 1, 2, 3, 4, 5
  # Potyviridae (max 13):   1, 5, 10
  # Tombusviridae (max 14): 1, 5, 10
  family_breaks <-
    balanced_breadth_breaks_v211(
      family_max
    )

  scale_offset <-
    scale_y_offset_v211(
      length(
        tr_full$tip.label
      ),
      minimum = 9
    )

  cat(
    "Supplementary ",
    fam,
    " Top dominant host orders: ",
    paste(
      family_orders,
      collapse = ", "
    ),
    " | breadth max = ",
    family_max,
    " | breadth legend = ",
    paste(
      family_breaks,
      collapse = ", "
    ),
    " | scale-bar y offset = ",
    scale_offset,
    "\n",
    sep = ""
  )

  pp <-
    make_endpoint_phylogram_v210(
      family =
        fam,
      tr =
        tr_full,
      trait =
        trait_full,
      top_host_orders =
        family_orders,
      order_palette =
        family_cols,
      show_legend =
        TRUE,
      folded =
        FALSE,
      breadth_breaks =
        family_breaks,
      # All four Supplementary trees use the SAME size mapping up to 17.
      breadth_limit_max =
        SUPP_BREADTH_MAX,
      inside_legend =
        FALSE,
      add_scale_bar =
        TRUE,
      scale_y_offset =
        scale_offset,
      legend_position =
        "bottom"
    ) +
    theme(
      plot.margin = margin(
        6,
        18,
        24,
        4
      )
    )

  supp_tree_plots[[fam]] <- pp

  save_triplet(
    pp,
    paste0(
      "FigS3",
      supp_letters[i],
      "_",
      fam,
      "_full_phylogram_endpoint_traits"
    ),
    OUT_SUPP,
    width_mm = 180,
    height_mm = 152
  )

  # Diagnostic only. No automatic pruning/filtering.
  long_branch_qc[[fam]] <-
    terminal_branch_qc_v211(
      tr =
        tr_raw,
      family =
        fam,
      trait =
        trait_full,
      top_fraction =
        0.01
    )
}

long_branch_qc_tbl <-
  rbindlist(
    long_branch_qc,
    fill = TRUE
  )

fwrite(
  long_branch_qc_tbl,
  file.path(
    OUT_LOG,
    "FigS3_terminal_branch_top1pct_QC.tsv"
  ),
  sep = "\t"
)

# Long-format composite retained because Supplementary is allowed to be tall
# and each family intentionally keeps its own explanatory legend.
supp_tree_composite <-
  patchwork::wrap_plots(
    supp_tree_plots[
      supp_families
    ],
    ncol = 1
  )

save_triplet(
  supp_tree_composite,
  "FigS3B_E_four_family_phylogenies_v2_11",
  OUT_SUPP,
  width_mm = 190,
  height_mm = 565
)

# ============================================================
# Reproducibility manifests
# ============================================================
fwrite(
  input_manifest,
  file.path(OUT_LOG, "Fig3_v2_11_1_input_manifest.tsv"),
  sep = "\t"
)

fwrite(
  tree_manifest,
  file.path(OUT_LOG, "Fig3_v2_11_1_tree_manifest.tsv"),
  sep = "\t"
)

writeLines(
  c(
    paste0("Resolved_font\t", FONT_FAMILY),
    paste0("R_version\t", R.version.string),
    paste0("FIG3_ROOT\t", FIG3_ROOT)
  ),
  file.path(OUT_LOG, "Fig3_v2_11_1_environment.txt")
)

cat("Fig.3 publication plotting v2.11.1 completed: ", format(Sys.time()), "\n", sep = "")
cat("Main outputs: ", OUT_MAIN, "\n", sep = "")
cat("Supplementary outputs: ", OUT_SUPP, "\n", sep = "")
cat("Logs/manifests: ", OUT_LOG, "\n", sep = "")
