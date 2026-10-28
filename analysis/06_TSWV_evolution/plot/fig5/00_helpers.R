suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ape)
  library(ggtree)
  library(patchwork)
})

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/non-empty required input: ", path)
  }
  data.table::fread(
    path,
    sep = "\t",
    quote = "",
    na.strings = c("", "NA", "N/A", "na", "n/a"),
    showProgress = FALSE
  )
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  z <- tolower(clean_text(x))
  z %in% c("true", "t", "1", "yes", "y")
}

gene_display <- function(x) {
  ifelse(as.character(x) == "GnGc", "Gn/Gc", as.character(x))
}

gene_segment <- function(gene) {
  out <- c(
    "RdRp" = "L",
    "GnGc" = "M",
    "NSm" = "M",
    "N" = "S",
    "NSs" = "S"
  )
  unname(out[as.character(gene)])
}

save_triplet <- function(plot, stem, outdir, width_mm, height_mm, dpi = 600) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggsave(
    file.path(outdir, paste0(stem, ".pdf")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    file.path(outdir, paste0(stem, ".png")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    file.path(outdir, paste0(stem, ".tiff")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
}

parse_breakpoints <- function(x) {
  z <- clean_text(x)
  if (!nzchar(z)) return(numeric())

  vals <- suppressWarnings(
    as.numeric(unlist(strsplit(z, ",", fixed = TRUE)))
  )

  vals[is.finite(vals)]
}

parse_ufboot <- function(x) {
  x <- clean_text(x)

  vapply(
    x,
    function(z) {
      if (!nzchar(z)) return(NA_real_)

      nums <- regmatches(
        z,
        gregexpr("[0-9]+(?:\\.[0-9]+)?", z, perl = TRUE)
      )[[1]]

      nums <- suppressWarnings(as.numeric(nums))
      if (!length(nums)) return(NA_real_)

      if (length(nums) >= 2L) nums[2] else nums[1]
    },
    numeric(1)
  )
}

nice_scale_width <- function(tr) {
  dd <- ape::node.depth.edgelength(tr)
  nt <- length(tr$tip.label)
  mx <- max(dd[seq_len(nt)], na.rm = TRUE)

  if (!is.finite(mx) || mx <= 0) return(0.01)

  target <- mx / 5
  pow <- 10 ^ floor(log10(target))
  m <- target / pow

  nice_m <- if (m < 1.5) {
    1
  } else if (m < 3.5) {
    2
  } else if (m < 7.5) {
    5
  } else {
    10
  }

  nice_m * pow
}

make_mini_segment_tree <- function(tr, segment) {
  # IMPORTANT:
  # In ggtree circular layout, changing the tree-data y coordinate mainly
  # changes angular position, so it does NOT behave like a simple visual
  # up/down nudge.  Therefore the segment letter is now overlaid in
  # normalized panel coordinates, independent of the circular tree geometry.

  tree_plot <- ggtree(
    tr,
    layout = "circular",
    branch.length = "none"
  ) +
    geom_tree(
      linewidth = 0.24,
      colour = COL_TREE
    ) +
    TREE_THEME()

  label_plot <- ggplot() +
    annotate(
      "text",
      x = 0.50,
      y = 0.425,
      label = segment,
      family = FONT_FAMILY,
      fontface = "bold",
      size = TEXT_NORMAL / ggplot2::.pt,
      colour = unname(SEGMENT_COLORS[segment])
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(
      base_family = FONT_FAMILY
    ) +
    theme(
      plot.background = element_rect(
        fill = NA,
        colour = NA
      ),
      panel.background = element_rect(
        fill = NA,
        colour = NA
      )
    )

  tree_plot +
    patchwork::inset_element(
      label_plot,
      left = 0,
      bottom = 0,
      right = 1,
      top = 1,
      align_to = "panel",
      on_top = TRUE,
      clip = FALSE
    )
}

make_full_gene_tree <- function(
  tr,
  gene,
  show_support = TRUE,
  support_size = 0.26,
  tree_linewidth = 0.17
) {
  p <- ggtree(
    tr,
    layout = "circular",
    branch.length = "none"
  ) +
    geom_tree(
      linewidth = tree_linewidth,
      colour = COL_TREE
    )

  p$data$UFBoot <- NA_real_

  node_idx <- which(
    !p$data$isTip &
    !is.na(p$data$label)
  )

  if (length(node_idx)) {
    p$data$UFBoot[node_idx] <- parse_ufboot(
      p$data$label[node_idx]
    )
  }

  if (isTRUE(show_support)) {
    p <- p +
      ggtree::geom_point2(
        aes(
          subset =
            !isTip &
            is.finite(UFBoot) &
            UFBoot >= 95
        ),
        size = support_size,
        colour = COL_SUPPORT
      )
  }

  seg <- gene_segment(gene)

  p +
    ggtitle(gene_display(gene)) +
    TREE_THEME(
      title_colour = unname(
        SEGMENT_COLORS[seg]
      )
    ) +
    theme(
      plot.title = element_text(
        margin = margin(0, 0, -13.0, 0)
      ),
      plot.margin = margin(-1.0, 0.2, 0.2, 0.2)
    )
}


make_full_segment_tree <- function(
  tr,
  segment,
  show_labels = FALSE
) {
  p <- ggtree(
    tr,
    layout = "rectangular",
    branch.length = "branch.length"
  ) +
    geom_tree(
      linewidth = 0.25,
      colour = COL_TREE
    )

  p$data$UFBoot <- NA_real_

  node_idx <- which(
    !p$data$isTip &
    !is.na(p$data$label)
  )

  if (length(node_idx)) {
    p$data$UFBoot[node_idx] <- parse_ufboot(
      p$data$label[node_idx]
    )
  }

  p <- p +
    ggtree::geom_point2(
      aes(
        subset =
          !isTip &
          is.finite(UFBoot) &
          UFBoot >= 95
      ),
      size = 0.34,
      colour = COL_SUPPORT
    ) +
    ggtitle(segment) +
    TREE_THEME(
      title_colour = unname(
        SEGMENT_COLORS[segment]
      )
    )

  sw <- nice_scale_width(tr)

  # Put the scale bar in the lower-left corner rather than in the visual
  # centre between stacked trees.
  xmax <- max(p$data$x, na.rm = TRUE)
  ntip <- length(tr$tip.label)

  # Manual scale bar: explicit Open Sans text avoids font fallback from
  # ggtree's default scale-bar label rendering.
  scale_x0 <- 0.02 * xmax
  scale_y0 <- -2.15
  scale_tick_h <- 0.42

  scale_label <- format(
    sw,
    scientific = FALSE,
    trim = TRUE
  )

  p <- p +
    annotate(
      "segment",
      x = scale_x0,
      xend = scale_x0 + sw,
      y = scale_y0,
      yend = scale_y0,
      linewidth = 0.34,
      colour = "black"
    ) +
    annotate(
      "segment",
      x = scale_x0,
      xend = scale_x0,
      y = scale_y0 - scale_tick_h / 2,
      yend = scale_y0 + scale_tick_h / 2,
      linewidth = 0.34,
      colour = "black"
    ) +
    annotate(
      "segment",
      x = scale_x0 + sw,
      xend = scale_x0 + sw,
      y = scale_y0 - scale_tick_h / 2,
      yend = scale_y0 + scale_tick_h / 2,
      linewidth = 0.34,
      colour = "black"
    ) +
    annotate(
      "text",
      x = scale_x0 + sw / 2,
      y = scale_y0 - 0.62,
      label = scale_label,
      family = FONT_FAMILY,
      size = TEXT_TINY / ggplot2::.pt,
      colour = COL_DARK,
      vjust = 1
    ) +
    scale_y_continuous(
      limits = c(-4.6, ntip + 2.5),
      expand = expansion(mult = c(0, 0))
    )

  if (isTRUE(show_labels)) {
    p <- p +
      geom_tiplab(
        size = TEXT_TINY / ggplot2::.pt,
        family = FONT_FAMILY,
        align = FALSE,
        offset = sw * 0.12
      )
  }

  p
}


make_segment_metric_table <- function(
  seg_sum,
  include_incompatible = FALSE
) {
  required <- c(
    "segment_pair",
    "matched_isolate_n",
    "patristic_spearman_rho",
    "normalized_RF_distance",
    "high_support_split_jaccard",
    "incompatible_high_support_split_pair_n"
  )

  miss <- setdiff(required, names(seg_sum))

  if (length(miss)) {
    stop(
      "Segment summary lacks: ",
      paste(miss, collapse = ", ")
    )
  }

  pair_order <- c(
    "L_vs_M",
    "L_vs_S",
    "M_vs_S"
  )

  pair_label <- c(
    "L_vs_M" = "L–M",
    "L_vs_S" = "L–S",
    "M_vs_S" = "M–S"
  )

  ss <- copy(seg_sum)
  ss <- ss[match(pair_order, segment_pair)]

  metrics <- list(
    list(
      key = "patristic_spearman_rho",
      label = "Patristic ρ",
      digits = 3L
    ),
    list(
      key = "normalized_RF_distance",
      label = "Normalized RF",
      digits = 3L
    ),
    list(
      key = "high_support_split_jaccard",
      label = "High-support split Jaccard",
      digits = 3L
    )
  )

  if (isTRUE(include_incompatible)) {
    metrics <- c(
      metrics,
      list(
        list(
          key = "incompatible_high_support_split_pair_n",
          label = "Incompatible high-support split pairs",
          digits = 0L
        )
      )
    )
  }

  nr <- length(metrics)
  rows <- vector("list", nr * length(pair_order))
  k <- 0L

  for (ri in seq_along(metrics)) {
    mm <- metrics[[ri]]

    for (ci in seq_along(pair_order)) {
      k <- k + 1L
      val <- as.numeric(ss[[mm$key]][ci])

      lab <- if (mm$digits == 0L) {
        formatC(val, format = "f", digits = 0)
      } else {
        formatC(val, format = "f", digits = mm$digits)
      }

      rows[[k]] <- data.table(
        row = ri,
        col = ci,
        metric = mm$label,
        pair = pair_label[pair_order[ci]],
        label = lab
      )
    }
  }

  dt <- rbindlist(rows)
  dt[, y := rev(seq_len(nr))[row]]
  dt[, x := col]

  row_labels <- unique(
    dt[, .(row, metric, y)]
  )

  col_labels <- unique(
    dt[, .(col, pair)]
  )

  ggplot() +
    geom_text(
      data = col_labels,
      aes(
        x = col,
        y = nr + 0.72,
        label = pair
      ),
      family = BASE_FAMILY,
      fontface = "bold",
      size = STRIP_SIZE / ggplot2::.pt,
      colour = COL_DARK
    ) +
    geom_text(
      data = row_labels,
      aes(
        x = 0.58,
        y = y,
        label = metric
      ),
      family = BASE_FAMILY,
      hjust = 1,
      size = TEXT_SMALL / ggplot2::.pt,
      colour = COL_DARK
    ) +
    geom_text(
      data = dt,
      aes(
        x = x,
        y = y,
        label = label
      ),
      family = BASE_FAMILY,
      size = AXIS_TEXT_SIZE / ggplot2::.pt,
      colour = COL_DARK
    ) +
    annotate(
      "segment",
      x = 0.72,
      xend = 3.28,
      y = nr + 0.38,
      yend = nr + 0.38,
      linewidth = 0.34,
      colour = "grey70"
    ) +
    coord_cartesian(
      xlim = c(-1.18, 4.18),
      ylim = c(0.45, nr + 0.98),
      clip = "off"
    ) +
    theme_void(base_family = BASE_FAMILY) +
    theme(
      plot.margin = margin(0, 2, 0, 2)
    )
}


