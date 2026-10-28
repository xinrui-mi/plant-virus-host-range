suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ape)
  library(ggtree)
  library(patchwork)
})

source("00_paths.R")
source("00_theme.R")
source("00_helpers.R")

cat("===== Fig7 publication FINAL v1.3.6 =====\n")
cat("Plotting only; all analytical inputs and formal trees are frozen.\n")
cat("Main: biology -> segment discordance -> coding-region recombination/selection landscape.\n")
cat("Supplement: five full coding-region cladograms + full L/M/S branch-length trees.\n\n")

required <- c(
  GENE_SUMMARY,
  SEGMENT_SUMMARY,
  GARD_STATUS,
  SELECTION_SUMMARY,
  SELECTION_SITES
)

for (f in required) {
  if (!file.exists(f) || file.info(f)$size <= 0) {
    stop("Missing/non-empty required input: ", f)
  }
}


gene_order <- c(
  "RdRp",
  "GnGc",
  "NSm",
  "N",
  "NSs"
)

segment_order <- c(
  "L",
  "M",
  "S"
)

gene_summary <- read_tsv(GENE_SUMMARY)
seg_sum <- read_tsv(SEGMENT_SUMMARY)
gard <- read_tsv(GARD_STATUS)
sel <- read_tsv(SELECTION_SUMMARY)
sites <- read_tsv(SELECTION_SITES)

# ============================================================
# Input gates.
# ============================================================
for (g in gene_order) {
  if (nrow(gene_summary[gene == g]) != 1L) {
    stop("Gene summary missing/duplicated: ", g)
  }

  if (nrow(gard[gene == g]) != 1L) {
    stop("GARD status missing/duplicated: ", g)
  }

  if (nrow(sel[gene == g]) != 1L) {
    stop("Selection summary missing/duplicated: ", g)
  }

  tf <- file.path(
    GENE_TREE_ROOT,
    g,
    paste0(g, ".codon_partitioned.treefile")
  )

  if (!file.exists(tf) || file.info(tf)$size <= 0) {
    stop("Missing gene tree: ", tf)
  }
}

for (s in segment_order) {
  tf <- file.path(
    SEGMENT_TREE_ROOT,
    s,
    paste0(s, ".treefile")
  )

  if (!file.exists(tf) || file.info(tf)$size <= 0) {
    stop("Missing segment tree: ", tf)
  }
}

if (
  !all(
    gard$status %in% c(
      "FINAL_GARD",
      "PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1",
      "PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1"
    )
  )
) {
  stop("Unexpected GARD final status in table 16.")
}

if (
  !all(sel$BUSTED_status == "present_valid") ||
  !all(sel$MEME_status == "present_valid") ||
  !all(sel$FUBAR_status == "present_valid")
) {
  stop("Selection summary contains non-valid method status.")
}

# ============================================================
# Shared five-gene coordinates and formal selection values.
# ============================================================
gard_plot <- gard[match(gene_order, gene)]
gard_plot[, y := rev(seq_len(.N))]
gard_plot[, segment := gene_segment(gene)]
gard_plot[, alignment_nt := as.numeric(alignment_nt)]
gard_plot[, final_breakpoint_n := as.integer(final_breakpoint_n)]
gard_plot[, display_gene := gene_display(gene)]

sel_plot <- sel[match(gene_order, gene)]
sel_plot[, gene_key := gene]
sel_plot[, display_gene := gene_display(gene_key)]
sel_plot[, segment := gene_segment(gene_key)]
sel_plot[, y := gard_plot$y[match(gene_key, gard_plot$gene)]]

sel_plot[, MEME_pct := 100 *
  as.numeric(MEME_sites_p_le_0.10) /
  as.numeric(MEME_sites_parsed_n)
]

sel_plot[, FUBAR_pos_pct := 100 *
  as.numeric(FUBAR_positive_sites_posterior_ge_0.90) /
  as.numeric(FUBAR_sites_parsed_n)
]

sel_plot[, FUBAR_neg_pct := 100 *
  as.numeric(FUBAR_negative_sites_posterior_ge_0.90) /
  as.numeric(FUBAR_sites_parsed_n)
]

# ============================================================
# Publication plotting integrity checks.
# These checks do not recalculate any Fig.7 analysis. They only ensure that
# frozen analytical tables can be represented without missing pairs or silent
# display-range clipping.
# ============================================================
required_segment_pairs <- c(
  "L_vs_M",
  "L_vs_S",
  "M_vs_S"
)

pair_counts <- seg_sum[
  segment_pair %in% required_segment_pairs,
  .N,
  by = segment_pair
]

if (
  nrow(pair_counts) != length(required_segment_pairs) ||
  !setequal(pair_counts$segment_pair, required_segment_pairs) ||
  any(pair_counts$N != 1L)
) {
  stop(
    "Segment-summary plotting gate failed: L_vs_M, L_vs_S and M_vs_S ",
    "must each occur exactly once."
  )
}

matched_isolate_values <- unique(
  as.integer(
    seg_sum[
      segment_pair %in% required_segment_pairs,
      matched_isolate_n
    ]
  )
)
matched_isolate_values <- matched_isolate_values[
  is.finite(matched_isolate_values)
]

if (length(matched_isolate_values) != 1L) {
  stop(
    "Segment-summary plotting gate failed: the three L/M/S pairwise ",
    "comparisons must use one identical matched-isolate count."
  )
}
matched_n_validated <- matched_isolate_values[[1]]

check_display_range <- function(x, upper, label) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(label, " contains non-finite values.")
  }
  if (any(x < 0 | x > upper)) {
    stop(
      label,
      " exceeds the frozen publication display range 0-",
      upper,
      "%. Refuse silent clipping; revise the plotting scale explicitly."
    )
  }
}

check_display_range(sel_plot$MEME_pct, 10, "MEME episodic selected-site fraction")
check_display_range(sel_plot$FUBAR_pos_pct, 2, "FUBAR diversifying selected-site fraction")
check_display_range(sel_plot$FUBAR_neg_pct, 80, "FUBAR purifying selected-site fraction")

# ============================================================
# MAIN A
# Minimal TSWV genome framework + ring-shaped transmission flow.
# ============================================================
genome_blocks <- data.table(
  segment = c(
    "L",
    "M", "M",
    "S", "S"
  ),
  y = c(
    3,
    2, 2,
    1, 1
  ),
  xmin = c(
    18,
    10, 66,
    10, 65
  ),
  xmax = c(
    90,
    54, 93,
    41, 93
  ),
  gene = c(
    "RdRp",
    "Gn/Gc", "NSm",
    "N", "NSs"
  ),
  function_label = c(
    "Replication",
    "Glycoproteins",
    "Movement protein",
    "Nucleocapsid",
    "Silencing suppressor"
  )
)

seg_lines <- data.table(
  segment = c("L", "M", "S"),
  y = c(3, 2, 1),
  xmin = 7,
  xmax = 98,
  genome_type = c(
    "negative-sense",
    "ambisense",
    "ambisense"
  )
)

pA_genome <- ggplot() +
  geom_segment(
    data = seg_lines,
    aes(
      x = xmin,
      xend = xmax,
      y = y,
      yend = y,
      colour = segment
    ),
    linewidth = 2.0,
    lineend = "round"
  ) +
  geom_rect(
    data = genome_blocks,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = y - 0.18,
      ymax = y + 0.18,
      fill = segment
    ),
    colour = "white",
    linewidth = 0.30
  ) +
  geom_text(
    data = genome_blocks,
    aes(
      x = (xmin + xmax) / 2,
      y = y,
      label = gene
    ),
    family = FONT_FAMILY,
    fontface = "bold",
    size = TEXT_NORMAL / ggplot2::.pt,
    colour = "white"
  ) +
  geom_text(
    data = genome_blocks,
    aes(
      x = (xmin + xmax) / 2,
      y = y - 0.31,
      label = function_label
    ),
    family = FONT_FAMILY,
    size = TEXT_TINY / ggplot2::.pt,
    colour = COL_MID
  ) +
  geom_text(
    data = seg_lines,
    aes(
      x = 2.2,
      y = y,
      label = paste0(
        segment,
        "  ",
        genome_type
      ),
      colour = segment
    ),
    family = FONT_FAMILY,
    fontface = "bold",
    hjust = 1,
    size = TEXT_SMALL / ggplot2::.pt
  ) +
  scale_colour_manual(
    values = SEGMENT_COLORS,
    guide = "none"
  ) +
  scale_fill_manual(
    values = SEGMENT_COLORS,
    guide = "none"
  ) +
  coord_cartesian(
    # Wider left limit prevents the segment/genome-type labels being clipped.
    xlim = c(-18, 101),
    ylim = c(0.35, 3.55),
    clip = "off"
  ) +
  theme_void(
    base_family = FONT_FAMILY
  ) +
  theme(
    plot.margin = margin(3, 4, 3, 8)
  )

# Five-node plant-to-plant cycle.
flow_angles <- c(
  150,
  78,
  0,
  -78,
  -150
) * pi / 180

flow_nodes <- data.table(
  node = seq_along(flow_angles),
  angle = flow_angles,
  x = 1.08 * cos(flow_angles),
  y = 1.08 * sin(flow_angles),
  label = c(
    "Plant",
    "Larval\nacquisition",
    "Replication\nin thrips",
    "Adult\ntransmission",
    "Plant"
  )
)

# Arrow endpoints terminate before the node circles.
node_radius <- 0.115
edge_rows <- list()

for (i in 1:4) {
  x1 <- flow_nodes$x[i]
  y1 <- flow_nodes$y[i]
  x2 <- flow_nodes$x[i + 1]
  y2 <- flow_nodes$y[i + 1]

  dx <- x2 - x1
  dy <- y2 - y1
  d <- sqrt(dx^2 + dy^2)

  ux <- dx / d
  uy <- dy / d

  edge_rows[[i]] <- data.table(
    x = x1 + node_radius * ux,
    y = y1 + node_radius * uy,
    xend = x2 - node_radius * ux,
    yend = y2 - node_radius * uy
  )
}

flow_edges <- rbindlist(edge_rows)

# Manual label positions keep text outside the arrow ring.
flow_nodes[, `:=`(
  label_x = c(-1.48, 0.18, 1.48, 0.18, -1.48),
  label_y = c( 0.57, 1.38, 0.00,-1.38, -0.57),
  hjust = c(1.00, 0.50, 0.00, 0.50, 1.00),
  vjust = c(0.50, 0.00, 0.50, 1.00, 0.50)
)]

pA_vector <- ggplot() +
  geom_curve(
    data = flow_edges,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    curvature = -0.18,
    linewidth = 0.52,
    colour = "grey48",
    arrow = grid::arrow(
      length = grid::unit(0.92, "mm"),
      type = "closed"
    )
  ) +
  geom_point(
    data = flow_nodes,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 3.10,
    stroke = 0.52,
    fill = "white",
    colour = COL_DARK
  ) +
  geom_text(
    data = flow_nodes,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      hjust = hjust,
      vjust = vjust
    ),
    family = FONT_FAMILY,
    size = TEXT_SMALL / ggplot2::.pt,
    lineheight = 0.95,
    colour = COL_DARK
  ) +
  annotate(
    "text",
    x = 0,
    y = 0.03,
    label = "Persistent\npropagative\ntransmission",
    family = FONT_FAMILY,
    fontface = "bold",
    size = TEXT_NORMAL / ggplot2::.pt,
    lineheight = 0.92,
    colour = COL_DARK
  ) +
  coord_fixed(
    xlim = c(-1.88, 1.88),
    ylim = c(-1.72, 1.72),
    clip = "off"
  ) +
  theme_void(
    base_family = FONT_FAMILY
  ) +
  theme(
    plot.margin = margin(3, 6, 3, 6)
  )

pA <- wrap_plots(
  pA_genome,
  pA_vector,
  ncol = 2,
  widths = c(
    1.58,
    1.02
  )
)

# ============================================================

# MAIN B
# Segment-level phylogenetic discordance.
# No arbitrary L1-L6 display groups.
# ============================================================
seg_trees <- setNames(
  lapply(
    segment_order,
    function(s) {
      read.tree(
        file.path(
          SEGMENT_TREE_ROOT,
          s,
          paste0(s, ".treefile")
        )
      )
    }
  ),
  segment_order
)

tip_sets <- lapply(
  seg_trees,
  function(x) sort(x$tip.label)
)

if (
  length(
    unique(
      vapply(
        tip_sets,
        paste,
        collapse = "|",
        FUN.VALUE = ""
      )
    )
  ) != 1L
) {
  stop(
    "L/M/S segment trees do not contain identical isolate labels."
  )
}

mini_L <- make_mini_segment_tree(
  seg_trees$L,
  "L"
)

mini_M <- make_mini_segment_tree(
  seg_trees$M,
  "M"
)

mini_S <- make_mini_segment_tree(
  seg_trees$S,
  "S"
)

pB_trees <- wrap_plots(
  mini_L,
  mini_M,
  mini_S,
  ncol = 3,
  widths = c(1, 1, 1)
)

pB_table <- make_segment_metric_table(
  seg_sum,
  include_incompatible = FALSE
)

pB <- wrap_plots(
  pB_trees,
  pB_table,
  ncol = 1,
  heights = c(1.74, 0.62)
)

# ============================================================
# MAIN C
# Unified five-gene evolutionary landscape.
# A single coordinate system aligns WHERE signals occur with HOW MUCH
# signal is detected for each gene.
# ============================================================

# ----- Site-level signal positions -----
site_sel <- sites[
  truthy(selected) &
  (
    (
      method == "MEME" &
      statistic == "p_value"
    ) |
    (
      method == "FUBAR" &
      statistic == "posterior_positive"
    )
  )
]

site_sel[, codon := as.numeric(global_alignment_codon)]
site_sel <- site_sel[
  is.finite(codon)
]

site_tracks <- sel_plot[
  ,
  .(
    gene = gene_key,
    display_gene,
    segment,
    y,
    codon_n = as.numeric(
      MEME_sites_parsed_n
    )
  )
]

site_sel <- merge(
  site_sel,
  site_tracks[
    ,
    .(
      gene,
      display_gene,
      segment,
      y,
      codon_n
    )
  ],
  by = "gene",
  all.x = TRUE
)

if (
  any(!is.finite(site_sel$y)) ||
  any(!is.finite(site_sel$codon_n)) ||
  any(!site_sel$gene %in% gene_order)
) {
  stop(
    "Unified main C site-to-gene merge failed. ",
    "Formal gene keys must remain RdRp/GnGc/NSm/N/NSs."
  )
}

site_sel[, rel := 100 * codon / codon_n]

site_sel[
  ,
  method_display := fifelse(
    method == "MEME",
    "MEME episodic",
    "FUBAR diversifying"
  )
]

expected_meme <- sum(
  as.integer(
    sel$MEME_sites_p_le_0.10
  ),
  na.rm = TRUE
)

expected_fubar_pos <- sum(
  as.integer(
    sel$FUBAR_positive_sites_posterior_ge_0.90
  ),
  na.rm = TRUE
)

observed_meme <- nrow(
  site_sel[
    method_display == "MEME episodic"
  ]
)

observed_fubar_pos <- nrow(
  site_sel[
    method_display == "FUBAR diversifying"
  ]
)

if (
  observed_meme != expected_meme ||
  observed_fubar_pos != expected_fubar_pos
) {
  stop(
    "Unified main C selected-site count gate failed: ",
    "MEME observed/expected = ",
    observed_meme, "/", expected_meme,
    "; FUBAR+ observed/expected = ",
    observed_fubar_pos, "/", expected_fubar_pos
  )
}

# ----- GARD breakpoint positions -----
bp_site <- list()
bp_k <- 0L

for (i in seq_len(nrow(gard_plot))) {
  bps <- parse_breakpoints(
    gard_plot$raw_breakpoint_positions_nt[i]
  )

  if (!length(bps)) next

  for (bp in bps) {
    bp_k <- bp_k + 1L

    bp_site[[bp_k]] <- data.table(
      gene = gard_plot$gene[i],
      y = gard_plot$y[i],
      rel =
        100 *
        bp /
        gard_plot$alignment_nt[i],
      raw_nt = bp
    )
  }
}

bp_site <- if (length(bp_site)) {
  rbindlist(bp_site)
} else {
  data.table()
}

# ----- Map three selection summaries into fixed right-hand display bands -----
# The three bands use independent scientific scales, but share the same
# gene rows.  This avoids separate mini-panels and duplicate y-axes.
BAND_MEME_START <- 116
BAND_MEME_END   <- 140

BAND_FDIV_START <- 148
BAND_FDIV_END   <- 172

BAND_FPUR_START <- 180
BAND_FPUR_END   <- 204

map_band <- function(value, start, end, max_value) {
  start + (pmin(pmax(value, 0), max_value) / max_value) * (end - start)
}

summary_points <- rbindlist(
  list(
    sel_plot[
      ,
      .(
        gene = gene_key,
        display_gene,
        segment,
        y,
        metric = "MEME",
        value = MEME_pct,
        x0 = BAND_MEME_START,
        x = map_band(
          MEME_pct,
          BAND_MEME_START,
          BAND_MEME_END,
          10
        )
      )
    ],
    sel_plot[
      ,
      .(
        gene = gene_key,
        display_gene,
        segment,
        y,
        metric = "FUBAR_div",
        value = FUBAR_pos_pct,
        x0 = BAND_FDIV_START,
        x = map_band(
          FUBAR_pos_pct,
          BAND_FDIV_START,
          BAND_FDIV_END,
          2
        )
      )
    ],
    sel_plot[
      ,
      .(
        gene = gene_key,
        display_gene,
        segment,
        y,
        metric = "FUBAR_pur",
        value = FUBAR_neg_pct,
        x0 = BAND_FPUR_START,
        x = map_band(
          FUBAR_neg_pct,
          BAND_FPUR_START,
          BAND_FPUR_END,
          80
        )
      )
    ]
  )
)

# Manual scientific ticks for the three independent right-hand scales.
summary_ticks <- rbindlist(
  list(
    data.table(
      x = map_band(
        c(0, 5, 10),
        BAND_MEME_START,
        BAND_MEME_END,
        10
      ),
      label = c("0", "5", "10")
    ),
    data.table(
      x = map_band(
        c(0, 1, 2),
        BAND_FDIV_START,
        BAND_FDIV_END,
        2
      ),
      label = c("0", "1", "2")
    ),
    data.table(
      x = map_band(
        c(0, 40, 80),
        BAND_FPUR_START,
        BAND_FPUR_END,
        80
      ),
      label = c("0", "40", "80")
    )
  )
)

summary_baselines <- data.table(
  x = c(
    BAND_MEME_START,
    BAND_FDIV_START,
    BAND_FPUR_START
  ),
  xend = c(
    BAND_MEME_END,
    BAND_FDIV_END,
    BAND_FPUR_END
  )
)

summary_headers <- data.table(
  x = c(
    (BAND_MEME_START + BAND_MEME_END) / 2,
    (BAND_FDIV_START + BAND_FDIV_END) / 2,
    (BAND_FPUR_START + BAND_FPUR_END) / 2
  ),
  y = 5.72,
  label = c(
    "MEME episodic\n(%)",
    "FUBAR diversifying\n(%)",
    "FUBAR purifying\n(%)"
  )
)

left_axis_ticks <- data.table(
  x = c(0, 50, 100),
  label = c("0", "50", "100")
)

# ----- Unified Panel C -----
pC <- ggplot() +

  # Left: gene tracks.
  geom_segment(
    data = site_tracks,
    aes(
      x = 0,
      xend = 100,
      y = y,
      yend = y,
      colour = segment
    ),
    linewidth = 2.05,
    lineend = "round"
  ) +

  # MEME episodic sites.
  geom_segment(
    data = site_sel[
      method_display == "MEME episodic"
    ],
    aes(
      x = rel,
      xend = rel,
      y = y + 0.06,
      yend = y + 0.29
    ),
    linewidth = 0.43,
    colour = COL_EVENT
  ) +

  # FUBAR diversifying sites.
  geom_point(
    data = site_sel[
      method_display == "FUBAR diversifying"
    ],
    aes(
      x = rel,
      y = y - 0.18
    ),
    shape = 2,
    stroke = 0.70,
    size = 1.85,
    colour = COL_EVENT
  ) +

  # GARD breakpoints.
  geom_segment(
    data = bp_site,
    aes(
      x = rel,
      xend = rel,
      y = y - 0.36,
      yend = y + 0.36
    ),
    linewidth = 0.80,
    linetype = "dashed",
    colour = COL_BREAKPOINT
  ) +

  # Gene names: one shared y-axis for the whole panel.
  geom_text(
    data = site_tracks,
    aes(
      x = -5.0,
      y = y,
      label = display_gene
    ),
    hjust = 1,
    family = FONT_FAMILY,
    fontface = "bold",
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  # Right: selection summary lollipops.
  geom_segment(
    data = summary_points,
    aes(
      x = x0,
      xend = x,
      y = y,
      yend = y,
      colour = segment
    ),
    linewidth = 0.68,
    lineend = "round"
  ) +

  geom_point(
    data = summary_points,
    aes(
      x = x,
      y = y,
      colour = segment
    ),
    size = 2.15
  ) +

  # Visual separators between site landscape and summaries.
  annotate(
    "segment",
    x = 108,
    xend = 108,
    y = 0.45,
    yend = 5.45,
    linewidth = 0.32,
    colour = "grey82"
  ) +

  # Compact site-level key at the top of the left region.
  annotate(
    "segment",
    x = 5,
    xend = 5,
    y = 5.57,
    yend = 5.80,
    linewidth = 0.43,
    colour = COL_EVENT
  ) +

  annotate(
    "text",
    x = 8,
    y = 5.69,
    label = "MEME episodic",
    hjust = 0,
    family = FONT_FAMILY,
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  annotate(
    "point",
    x = 38,
    y = 5.69,
    shape = 2,
    size = 1.80,
    colour = COL_EVENT
  ) +

  annotate(
    "text",
    x = 41,
    y = 5.69,
    label = "FUBAR diversifying",
    hjust = 0,
    family = FONT_FAMILY,
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  annotate(
    "segment",
    x = 79,
    xend = 79,
    y = 5.54,
    yend = 5.82,
    linewidth = 0.80,
    linetype = "dashed",
    colour = COL_BREAKPOINT
  ) +

  annotate(
    "text",
    x = 82,
    y = 5.69,
    label = "GARD breakpoint",
    hjust = 0,
    family = FONT_FAMILY,
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  # Right-side summary headers.
  geom_text(
    data = summary_headers,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    family = FONT_FAMILY,
    fontface = "bold",
    lineheight = 0.95,
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  # Manual bottom axes for the three independent summary scales.
  geom_segment(
    data = summary_baselines,
    aes(
      x = x,
      xend = xend,
      y = 0.50,
      yend = 0.50
    ),
    inherit.aes = FALSE,
    linewidth = 0.38,
    colour = "black"
  ) +

  geom_segment(
    data = summary_ticks,
    aes(
      x = x,
      xend = x,
      y = 0.50,
      yend = 0.42
    ),
    inherit.aes = FALSE,
    linewidth = 0.34,
    colour = "black"
  ) +

  geom_text(
    data = summary_ticks,
    aes(
      x = x,
      y = 0.27,
      label = label
    ),
    inherit.aes = FALSE,
    family = FONT_FAMILY,
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  scale_colour_manual(
    values = SEGMENT_COLORS,
    guide = "none"
  ) +

  scale_x_continuous(
    limits = c(-14, 207),
    breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +

  scale_y_continuous(
    limits = c(-0.28, 5.95),
    breaks = NULL
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  PUB_THEME() +

  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 4, 5, 10)
  ) +

  # Local left axis for relative gene position only.
  annotate(
    "segment",
    x = 0,
    xend = 100,
    y = 0.50,
    yend = 0.50,
    linewidth = 0.38,
    colour = "black"
  ) +

  annotate(
    "segment",
    x = c(0, 50, 100),
    xend = c(0, 50, 100),
    y = 0.50,
    yend = 0.42,
    linewidth = 0.34,
    colour = "black"
  ) +

  geom_text(
    data = left_axis_ticks,
    aes(
      x = x,
      y = 0.27,
      label = label
    ),
    inherit.aes = FALSE,
    family = FONT_FAMILY,
    size = TEXT_SMALL / ggplot2::.pt,
    colour = COL_DARK
  ) +

  annotate(
    "text",
    x = 50,
    y = -0.10,
    label = "Relative gene position (%)",
    family = FONT_FAMILY,
    fontface = "bold",
    size = TEXT_NORMAL / ggplot2::.pt,
    colour = COL_DARK
  )

# ============================================================

# Save main panels.
# ============================================================
save_triplet(
  pA,
  "Fig7A_TSWV_genome_transmission_FINAL_v1_3_6",
  OUT_PANEL,
  180,
  55
)

save_triplet(
  pB,
  "Fig7B_segment_phylogenetic_discordance_FINAL_v1_3_6",
  OUT_PANEL,
  180,
  89
)

save_triplet(
  pC,
  "Fig7C_unified_gene_evolution_FINAL_v1_3_6",
  OUT_PANEL,
  180,
  74
)

# ============================================================
# Main Fig.7: biological context -> segment scale -> gene/site scale.
# ============================================================
main_final <- wrap_plots(
  pA,
  pB,
  pC,
  ncol = 1,
  heights = c(
    0.60,
    0.92,
    0.90
  )
)

save_triplet(
  main_final,
  "Fig7_TSWV_multiscale_evolution_FINAL_v1_3_6",
  OUT_MAIN,
  180,
  210
)

# ============================================================
# Supplementary Fig.S7
# Five equal-sized topology-oriented cladograms.
# No midpoint rooting and no host/geography colouring.
# ============================================================
gene_trees <- setNames(
  lapply(
    gene_order,
    function(g) {
      read.tree(
        file.path(
          GENE_TREE_ROOT,
          g,
          paste0(
            g,
            ".codon_partitioned.treefile"
          )
        )
      )
    }
  ),
  gene_order
)

s7_RdRp <- make_full_gene_tree(
  gene_trees[["RdRp"]],
  "RdRp",
  TRUE,
  0.28,
  0.17
)

s7_GnGc <- make_full_gene_tree(
  gene_trees[["GnGc"]],
  "GnGc",
  TRUE,
  0.28,
  0.17
)

s7_NSm <- make_full_gene_tree(
  gene_trees[["NSm"]],
  "NSm",
  TRUE,
  0.28,
  0.17
)

s7_N <- make_full_gene_tree(
  gene_trees[["N"]],
  "N",
  TRUE,
  0.28,
  0.17
)

s7_NSs <- make_full_gene_tree(
  gene_trees[["NSs"]],
  "NSs",
  TRUE,
  0.28,
  0.17
)

design_s7 <- "
AABB
CCDD
#EE#
"

figS7 <- wrap_plots(
  A = s7_RdRp,
  B = s7_GnGc,
  C = s7_NSm,
  D = s7_N,
  E = s7_NSs,
  design = design_s7,
  heights = c(
    1,
    1
  )
)

save_triplet(
  figS7,
  "FigS7_five_gene_phylogenies_FINAL_v1_3_6",
  OUT_SUPP,
  190,
  272
)

# Individual full-size cladograms remain available for detailed review.
for (g in gene_order) {
  pp <- make_full_gene_tree(
    gene_trees[[g]],
    g,
    show_support = TRUE,
    support_size = 0.28,
    tree_linewidth = 0.20
  )

  save_triplet(
    pp,
    paste0(
      "FigS7_",
      g,
      "_full_gene_cladogram_FINAL_v1_3_6"
    ),
    OUT_SUPP,
    178,
    178
  )
}

# ============================================================

# Supplementary Fig.S8
# Full branch-length L/M/S phylogenies + detailed discordance metrics.
# No arbitrary L1-L6 display groups.
# ============================================================
s8_plots <- setNames(
  lapply(
    segment_order,
    function(s) {
      make_full_segment_tree(
        seg_trees[[s]],
        s,
        show_labels = FALSE
      )
    }
  ),
  segment_order
)

s8_metric_table <- make_segment_metric_table(
  seg_sum,
  include_incompatible = TRUE
)

design_s8 <- "
L
M
S
T
"

figS8 <- wrap_plots(
  L = s8_plots[["L"]],
  M = s8_plots[["M"]],
  S = s8_plots[["S"]],
  T = s8_metric_table,
  design = design_s8,
  heights = c(
    1.20,
    1.20,
    1.20,
    0.58
  )
)

save_triplet(
  figS8,
  "FigS8_LMS_segment_phylogenies_FINAL_v1_3_6",
  OUT_SUPP,
  180,
  286
)

# ============================================================
# Plotting audit.
# ============================================================
matched_n <- matched_n_validated

audit <- data.table(
  metric = c(
    "status",
    "matched_LMS_isolate_n",
    "RdRp_GARD_breakpoint_n",
    "BUSTED_significant_gene_n",
    "MEME_selected_site_expected_n",
    "MEME_selected_site_plotted_mainC_n",
    "FUBAR_positive_site_expected_n",
    "FUBAR_positive_site_plotted_mainC_n",
    "FUBAR_negative_site_n",
    "mainC_integrates_site_and_selection_summary",
    "FigS7_host_continent_annotation",
    "FigS8_L_reference_display_groups"
  ),
  value = c(
    "PASS",
    matched_n,
    gard[
      gene == "RdRp",
      final_breakpoint_n
    ],
    sum(
      truthy(
        sel$BUSTED_significant_p_lt_0.05
      )
    ),
    expected_meme,
    observed_meme,
    expected_fubar_pos,
    observed_fubar_pos,
    sum(
      as.integer(
        sel$FUBAR_negative_sites_posterior_ge_0.90
      ),
      na.rm = TRUE
    ),
    "TRUE",
    "FALSE",
    "FALSE"
  )
)

fwrite(
  audit,
  file.path(
    OUT_LOG,
    "Fig7_FINAL_v1_3_6_plotting_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

cat("\nFig7 publication FINAL v1.3.6 plotting PASS\n")
cat("Main C integrates site-level signal positions and gene-level selection fractions.\n")
cat("Fig.S7 shows five coding-region topology-oriented cladograms with subtle UFBoot >=95 support dots.\n")
cat("Fig.S8 contains full branch-length L/M/S trees with increased vertical spacing and detailed topology metrics.\n")
cat("No midpoint rooting, host/geography colouring, or L1-L6 display groups are introduced.\n")
cat("All outputs saved as PDF, PNG and TIFF.\n")
