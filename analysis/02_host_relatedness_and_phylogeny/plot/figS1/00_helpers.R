suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ape)
  library(ggtree)
  library(patchwork)
})

find_col <- function(nms, exact = character(), regex = character(), required = TRUE) {
  low <- tolower(nms)

  for (x in exact) {
    hit <- which(low == tolower(x))
    if (length(hit)) return(nms[hit[1]])
  }

  for (r in regex) {
    hit <- grep(r, low, perl = TRUE)
    if (length(hit)) return(nms[hit[1]])
  }

  if (required) {
    stop(
      "Could not identify required column.\nAvailable columns:\n",
      paste(nms, collapse = ", ")
    )
  }

  NULL
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
}

save_triplet <- function(plot, stem, outdir, width_mm, height_mm, dpi = 600) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  pdf_path  <- file.path(outdir, paste0(stem, ".pdf"))
  tiff_path <- file.path(outdir, paste0(stem, ".tiff"))
  png_path  <- file.path(outdir, paste0(stem, ".png"))

  ggsave(
    pdf_path, plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white", limitsize = FALSE
  )

  ggsave(
    tiff_path, plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, compression = "lzw",
    bg = "white", limitsize = FALSE
  )

  ggsave(
    png_path, plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi,
    bg = "white", limitsize = FALSE
  )

  invisible(c(pdf_path, tiff_path, png_path))
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}

distance_level_map <- c(
  "0" = "Same species",
  "1" = "Same genus",
  "2" = "Same family",
  "3" = "Same order",
  "4" = "Same class",
  "5" = "Same phylum",
  "6" = "Same kingdom"
)

standardize_distance_table <- function(x) {
  nms <- names(x)

  dcol <- find_col(
    nms,
    exact = c(
      "taxonomy_distance",
      "host_distance",
      "distance",
      "host_taxonomic_distance",
      "host_phylogenetic_distance"
    ),
    regex = c(
      "^taxonomy.*distance$",
      "host.*distance",
      "^distance$"
    )
  )

  pcol <- find_col(
    nms,
    exact = c(
      "sharing_probability",
      "shared_probability",
      "probability",
      "shared_virus_probability"
    ),
    regex = c("sharing.*prob", "shared.*prob", "^prob"),
    required = FALSE
  )

  if (is.null(pcol)) {
    shared_col <- find_col(
      nms,
      exact = c(
        "positive_pairs",
        "shared_pairs",
        "shared_pair_n",
        "shared_pairs_n"
      ),
      regex = c(
        "^positive.*pairs?$",
        "shared.*pair.*n",
        "shared.*n"
      )
    )
    pair_col <- find_col(
      nms,
      exact = c(
        "total_pairs",
        "pair_n",
        "host_pair_n",
        "host_pair_count"
      ),
      regex = c(
        "^total.*pairs?$",
        "pair.*(n|count)"
      )
    )
    x[, .plot_probability := as.numeric(get(shared_col)) / as.numeric(get(pair_col))]
    pcol <- ".plot_probability"
  }

  pair_col <- find_col(
    nms,
    exact = c(
      "total_pairs",
      "pair_n",
      "host_pair_n",
      "host_pair_count",
      "host_pairs_n"
    ),
    regex = c(
      "^total.*pairs?$",
      "^pair.*(n|count)$",
      "host.*pair.*(n|count)"
    ),
    required = FALSE
  )

  x[, .distance := as.numeric(get(dcol))]
  x[, .prob := as.numeric(get(pcol))]

  if (!is.null(pair_col)) {
    x[, .pair_n := as.numeric(get(pair_col))]
  } else {
    x[, .pair_n := NA_real_]
  }

  x <- x[is.finite(.distance) & is.finite(.prob)]
  setorder(x, .distance)

  x[, .distance_label := factor(
    as.character(.distance),
    levels = names(distance_level_map),
    labels = unname(distance_level_map)
  )]

  if (all(is.na(x$.distance_label))) {
    x[, .distance_label := factor(.distance)]
  }

  x
}

make_distance_decay_plot <- function(
    x,
    line_colour,
    y_title,
    show_pair_size = TRUE,
    legend_inside = TRUE
) {
  x <- standardize_distance_table(copy(x))

  ymax <- max(x$.prob, na.rm = TRUE)
  upper_break <- max(0.2, ceiling(ymax / 0.2) * 0.2)
  y_breaks <- seq(0, upper_break, by = 0.2)

  has_pairs <- show_pair_size && any(is.finite(x$.pair_n) & x$.pair_n > 0)

  p <- ggplot(
    x,
    aes(x = .distance_label, y = .prob, group = 1)
  ) +
    geom_line(
      linewidth = 0.72,
      colour = line_colour
    )

  if (has_pairs) {
    pair_vals <- x$.pair_n[is.finite(x$.pair_n) & x$.pair_n > 0]
    pmin <- min(pair_vals)
    pmax <- max(pair_vals)

    power_seq <- seq(
      floor(log10(pmin)),
      ceiling(log10(pmax)),
      by = 1
    )
    pair_breaks <- 10^power_seq
    pair_breaks <- pair_breaks[
      pair_breaks >= pmin * 0.8 &
      pair_breaks <= pmax * 1.25
    ]

    if (length(pair_breaks) < 2) {
      pair_breaks <- pretty(pair_vals, n = 3)
      pair_breaks <- pair_breaks[pair_breaks > 0]
    }

    p <- p +
      geom_point(
        data = x,
        aes(size = .pair_n),
        shape = 21,
        stroke = 0.58,
        fill = "white",
        colour = line_colour
      ) +
      scale_size_continuous(
        trans = "log10",
        range = c(1.6, 4.2),
        name = "Host-pair count",
        breaks = pair_breaks,
        labels = scales::label_number(
          accuracy = 1,
          big.mark = ","
        )
      ) +
      guides(
        size = guide_legend(
          title.position = "top",
          title.hjust = 0,
          override.aes = list(
            fill = "white",
            colour = line_colour,
            stroke = 0.55
          )
        )
      )
  } else {
    p <- p +
      geom_point(
        size = 2.3,
        shape = 21,
        stroke = 0.58,
        fill = "white",
        colour = line_colour
      )
  }

  base_plot <- p +
    scale_y_continuous(
      breaks = y_breaks,
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.05, 0.07))
    ) +
    labs(
      x = "Host taxonomic distance",
      y = y_title
    ) +
    PUB_THEME() +
    theme(
      axis.text.x = element_text(
        angle = 23,
        hjust = 1,
        vjust = 1,
        face = "plain"
      )
    )

  if (has_pairs && legend_inside) {
    # Numeric legend.position is retained for broad ggplot2 compatibility.
    # It places the size legend in the otherwise empty upper-right region.
    base_plot <- base_plot +
      theme(
        legend.position = c(0.76, 0.67),
        legend.justification = c(0.5, 0.5),
        legend.background = element_rect(
          fill = scales::alpha("white", 0.88),
          colour = NA
        ),
        legend.margin = margin(1.5, 1.5, 1.5, 1.5),
        legend.key.height = grid::unit(4.4, "mm"),
        legend.key.width = grid::unit(4.4, "mm")
      )
  } else {
    base_plot <- base_plot +
      theme(
        legend.position = if (has_pairs) "right" else "none"
      )
  }

  base_plot
}

find_tree_for_family <- function(fig3_root, family) {
  roots <- unique(c(
    file.path(fig3_root, "work"),
    file.path(fig3_root, "pipeline"),
    file.path(fig3_root, "tables"),
    fig3_root
  ))

  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(NA_character_)

  allf <- unique(unlist(lapply(
    roots,
    function(z) list.files(
      z,
      pattern = "\\.(treefile|nwk|newick|tree|tre)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  )))

  if (!length(allf)) return(NA_character_)

  fam_full <- tolower(family)
  fam_stem <- sub("idae$", "", fam_full)

  hit <- allf[
    grepl(fam_full, tolower(allf), fixed = TRUE) |
      grepl(fam_stem, tolower(allf), fixed = TRUE)
  ]

  if (!length(hit)) return(NA_character_)

  low <- tolower(hit)
  score <- integer(length(hit))
  score <- score + 5L * grepl("\\.treefile$", low)
  score <- score + 3L * grepl("final|formal", low)
  score <- score + 2L * grepl("iqtree|ml", low)
  score <- score + 1L * grepl("trim", low)

  hit[order(-score, nchar(hit))][1]
}

find_tip_host_family_file <- function(fig3_root) {
  roots <- unique(c(
    file.path(fig3_root, "tables"),
    file.path(fig3_root, "work"),
    file.path(fig3_root, "pipeline"),
    file.path(dirname(fig3_root), "inputs")
  ))
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(NA_character_)

  cand <- unique(unlist(lapply(
    roots,
    function(z) list.files(
      z,
      pattern = "\\.(tsv|txt|csv)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  )))

  if (!length(cand)) return(NA_character_)

  # Avoid obvious pairwise/permutation tables for this task.
  cand <- cand[
    !grepl(
      "pair_sampling|permutation|null_beta|iterations|jaccard.*pair",
      basename(cand),
      ignore.case = TRUE
    )
  ]

  out <- list()

  for (f in cand) {
    hdr <- tryCatch(
      fread(f, nrows = 0, showProgress = FALSE),
      error = function(e) NULL
    )
    if (is.null(hdr)) next

    nms <- names(hdr)
    low <- tolower(nms)

    id_ok <- any(low %in% c(
      "representative_sequence_id",
      "tip",
      "tip_id",
      "label",
      "sequence_id"
    ))
    hf_ok <- any(low == "host_family")

    if (!id_ok || !hf_ok) next

    fam_ok <- any(low %in% c("virus_family", "family"))
    score <- 0L
    score <- score + 10L * fam_ok
    score <- score + 5L * grepl("03e|03f|tip|annotation|host", basename(f), ignore.case = TRUE)
    score <- score - as.integer(file.info(f)$size / 1e9)

    out[[length(out) + 1L]] <- data.table(file = f, score = score)
  }

  if (!length(out)) return(NA_character_)

  tab <- rbindlist(out)
  tab[order(-score, nchar(file))]$file[1]
}

read_tip_host_family_long <- function(path) {
  if (is.na(path) || !file.exists(path)) return(data.table())

  x <- fread(path, showProgress = FALSE)
  nms <- names(x)

  idcol <- find_col(
    nms,
    exact = c(
      "representative_sequence_id",
      "tip",
      "tip_id",
      "label",
      "sequence_id"
    ),
    regex = c("representative.*sequence", "^tip", "^label$", "sequence.*id")
  )

  hfcol <- find_col(
    nms,
    exact = c("host_family"),
    regex = c("^host.*family$")
  )

  famcol <- find_col(
    nms,
    exact = c("virus_family"),
    regex = c("^virus.*family$"),
    required = FALSE
  )

  y <- data.table(
    label = as.character(x[[idcol]]),
    host_family = as.character(x[[hfcol]])
  )

  if (!is.null(famcol)) {
    y[, virus_family := as.character(x[[famcol]])]
  }

  y <- y[
    !is.na(label) & nzchar(label) &
      !is.na(host_family) & nzchar(host_family)
  ]

  unique(y)
}

extract_tip_breadth_from_pair_table <- function(table_dir) {
  cand <- file.path(table_dir, "25_03F_phylogeny_host_distance_pairs.tsv")
  if (!file.exists(cand)) {
    allf <- list.files(
      table_dir,
      pattern = "03F.*\\.tsv$|25_03F.*\\.tsv$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (!length(allf)) return(data.table())
    cand <- allf[1]
  }

  x <- tryCatch(fread(cand, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(data.table())

  req <- c("tip1", "tip2", "host_breadth_1", "host_breadth_2")
  if (!all(req %in% names(x))) return(data.table())

  a <- x[, .(
    label = as.character(tip1),
    host_family_breadth = as.numeric(host_breadth_1)
  )]
  b <- x[, .(
    label = as.character(tip2),
    host_family_breadth = as.numeric(host_breadth_2)
  )]

  if ("virus_family" %in% names(x)) {
    a[, virus_family := as.character(x$virus_family)]
    b[, virus_family := as.character(x$virus_family)]
  }

  y <- rbind(a, b, fill = TRUE)
  y <- y[is.finite(host_family_breadth)]

  if ("virus_family" %in% names(y)) {
    y[, .(host_family_breadth = max(host_family_breadth)), by = .(virus_family, label)]
  } else {
    y[, .(host_family_breadth = max(host_family_breadth)), by = label]
  }
}

make_tree_host_matrix <- function(
    family,
    tree_file,
    host_long,
    top_n = 5,
    add_panel_tag = NULL,
    family_colour = "#3182BD"
) {
  if (is.na(tree_file) || !file.exists(tree_file)) {
    stop("Tree file not found for ", family)
  }

  tr <- ape::read.tree(tree_file)
  tr <- ape::ladderize(tr, right = FALSE)

  hf <- copy(host_long)

  if ("virus_family" %in% names(hf)) {
    tmp <- hf[virus_family == family]
    if (nrow(tmp)) hf <- tmp
  }

  hf <- hf[label %in% tr$tip.label]

  if (!nrow(hf)) {
    gp <- suppressWarnings(
      ggtree::ggtree(tr, layout = "rectangular")
    ) +
      theme_void(base_family = FONT_FAMILY) +
      labs(
        title = paste0(family, " (n = ", length(tr$tip.label), ")")
      ) +
      theme(
        plot.title = element_text(
          family = FONT_FAMILY,
          size = TEXT_SMALL,
          face = "italic",
          hjust = 0
        ),
        plot.margin = margin(2, 2, 2, 2)
      )

    if (!is.null(add_panel_tag)) {
      gp <- gp + labs(tag = add_panel_tag) + PANEL_TAG_THEME
    }

    return(gp)
  }

  hf <- unique(hf[, .(label, host_family)])

  hf_count <- hf[, .(tip_n = uniqueN(label)), by = host_family]
  setorder(hf_count, -tip_n, host_family)
  top_hosts <- head(hf_count$host_family, top_n)

  tree_plot <- suppressWarnings(
    ggtree::ggtree(tr, layout = "rectangular")
  ) +
    labs(
      title = paste0(family, " (n = ", length(tr$tip.label), ")")
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 1.5)
      ),
      plot.margin = margin(2, 0, 2, 2)
    )

  if (!is.null(add_panel_tag)) {
    tree_plot <- tree_plot +
      labs(tag = add_panel_tag) +
      PANEL_TAG_THEME
  }

  tip_y <- as.data.table(tree_plot$data)[
    isTip == TRUE,
    .(
      label = as.character(label),
      y = as.numeric(y)
    )
  ]

  presence <- unique(
    hf[host_family %in% top_hosts, .(label, host_family)]
  )
  presence[, present := 1L]

  grid <- CJ(
    label = tr$tip.label,
    host_family = top_hosts,
    unique = TRUE
  )
  grid <- merge(
    grid,
    presence,
    by = c("label", "host_family"),
    all.x = TRUE
  )
  grid[is.na(present), present := 0L]
  grid <- merge(
    grid,
    tip_y,
    by = "label",
    all.x = TRUE
  )

  grid <- grid[is.finite(y)]
  grid[, host_family := factor(
    host_family,
    levels = top_hosts
  )]

  y_rng <- range(tip_y$y, finite = TRUE)

  heat <- ggplot(
    grid,
    aes(
      x = host_family,
      y = y,
      fill = factor(present)
    )
  ) +
    geom_tile(
      width = 0.92,
      height = 0.92
    ) +
    scale_fill_manual(
      values = c(
        "0" = "white",
        "1" = family_colour
      ),
      guide = "none"
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    coord_cartesian(
      ylim = c(y_rng[1] - 0.55, y_rng[2] + 0.55),
      expand = FALSE,
      clip = "off"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      axis.text.x = element_text(
        size = TEXT_TINY,
        colour = "black",
        angle = 52,
        hjust = 0,
        vjust = 0.5,
        face = "plain"
      ),
      axis.text.y = element_blank(),
      plot.margin = margin(12, 2, 2, 1)
    )

  # Explicit layout prevents patchwork from flattening the nested tree+matrix
  # object when it is later combined with another representative family.
  patchwork::wrap_plots(
    list(tree_plot, heat),
    ncol = 2,
    widths = c(3.8, 1.45)
  )
}

extract_observed_beta <- function(global_model_table) {
  x <- read_tsv(global_model_table)
  beta_col <- find_col(
    names(x),
    exact = c("beta", "estimate", "coefficient"),
    regex = c("beta", "estimate", "coef")
  )
  as.numeric(x[[beta_col]][1])
}

extract_empirical_p <- function(global_model_table) {
  x <- read_tsv(global_model_table)
  pcol <- find_col(
    names(x),
    exact = c(
      "empirical_p",
      "empirical_p_value",
      "permutation_p",
      "p_empirical"
    ),
    regex = c("empirical.*p", "permut.*p"),
    required = FALSE
  )
  if (is.null(pcol)) return(NA_real_)
  as.numeric(x[[pcol]][1])
}


# ============================================================
# v2.3 helpers
# ============================================================

tag_strip <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.5,
      label = label,
      hjust = 0,
      vjust = 0.5,
      family = FONT_FAMILY,
      fontface = "bold",
      size = 4.0,
      colour = "black"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.margin = margin(0, 0, 0, 1)
    )
}

add_external_panel_tag <- function(plot, label, tag_height = 0.08) {
  patchwork::wrap_plots(
    list(tag_strip(label), patchwork::wrap_elements(full = plot)),
    ncol = 1,
    heights = c(tag_height, 1 - tag_height)
  )
}

find_tip_annotation_master <- function(fig3_root) {
  candidates <- c(
    file.path(fig3_root, "tables", "21_03E_tree_tip_annotation_master.tsv"),
    file.path(fig3_root, "tables", "22_03E_cluster_member_host_long.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[1])

  allf <- list.files(
    file.path(fig3_root, "tables"),
    pattern = "03E.*annotation.*\\.tsv$|tree_tip.*\\.tsv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(allf)) return(allf[1])
  NA_character_
}

read_tip_annotation_master <- function(path) {
  if (is.na(path) || !file.exists(path)) return(data.table())
  fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
}

normalize_tip_id_col <- function(x) {
  idcol <- find_col(
    names(x),
    exact = c(
      "representative_sequence_id",
      "tip",
      "tip_id",
      "label",
      "sequence_id",
      "representative_id"
    ),
    regex = c(
      "representative.*sequence",
      "^tip",
      "^label$",
      "sequence.*id",
      "representative.*id"
    )
  )
  x[, .tip_id := as.character(get(idcol))]
  x
}

extract_tip_genus <- function(ann) {
  if (!nrow(ann)) return(data.table())
  x <- normalize_tip_id_col(copy(ann))
  gcol <- find_col(
    names(x),
    exact = c(
      "virus_genus",
      "genus",
      "current_genus",
      "ictv_genus"
    ),
    regex = c(
      "^virus.*genus$",
      "^genus$",
      "ictv.*genus",
      "current.*genus"
    ),
    required = FALSE
  )
  if (is.null(gcol)) return(data.table())
  y <- x[
    !is.na(.tip_id) & nzchar(.tip_id),
    .(
      virus_genus = {
        z <- as.character(get(gcol))
        z <- z[!is.na(z) & nzchar(z)]
        if (length(z)) names(sort(table(z), decreasing = TRUE))[1] else NA_character_
      }
    ),
    by = .tip_id
  ]
  y
}

extract_tip_host_order <- function(ann, host_long = data.table()) {
  candidates <- list()

  if (nrow(ann)) {
    x <- normalize_tip_id_col(copy(ann))
    ocol <- find_col(
      names(x),
      exact = c(
        "dominant_host_order",
        "host_order",
        "host_order_name"
      ),
      regex = c(
        "dominant.*host.*order",
        "^host.*order$"
      ),
      required = FALSE
    )
    if (!is.null(ocol)) {
      candidates[[length(candidates) + 1L]] <- x[
        !is.na(.tip_id) & nzchar(.tip_id),
        .(
          host_order = {
            z <- as.character(get(ocol))
            z <- z[!is.na(z) & nzchar(z)]
            if (length(z)) names(sort(table(z), decreasing = TRUE))[1] else NA_character_
          }
        ),
        by = .tip_id
      ]
    }
  }

  if (nrow(host_long)) {
    hl <- copy(host_long)
    idcol <- find_col(
      names(hl),
      exact = c("label", "representative_sequence_id", "tip", "tip_id", "sequence_id"),
      regex = c("^label$", "representative.*sequence", "^tip", "sequence.*id"),
      required = FALSE
    )
    ocol <- find_col(
      names(hl),
      exact = c("host_order", "dominant_host_order"),
      regex = c("^host.*order$", "dominant.*host.*order"),
      required = FALSE
    )
    if (!is.null(idcol) && !is.null(ocol)) {
      candidates[[length(candidates) + 1L]] <- hl[
        !is.na(get(idcol)) & nzchar(as.character(get(idcol))),
        .(
          host_order = {
            z <- as.character(get(ocol))
            z <- z[!is.na(z) & nzchar(z)]
            if (length(z)) names(sort(table(z), decreasing = TRUE))[1] else NA_character_
          }
        ),
        by = .( .tip_id = as.character(get(idcol)) )
      ]
    }
  }

  if (!length(candidates)) return(data.table())
  y <- rbindlist(candidates, fill = TRUE)
  y <- y[
    !is.na(.tip_id) & nzchar(.tip_id),
    .(
      host_order = {
        z <- host_order[!is.na(host_order) & nzchar(host_order)]
        if (length(z)) z[1] else NA_character_
      }
    ),
    by = .tip_id
  ]
  y
}

extract_tip_breadth <- function(ann, host_long) {
  # Prefer a direct numeric breadth field in the formal tip annotation table.
  if (nrow(ann)) {
    x <- normalize_tip_id_col(copy(ann))
    bcol <- find_col(
      names(x),
      exact = c(
        "host_family_breadth",
        "host_family_n",
        "host_family_count",
        "n_host_family",
        "n_host_families"
      ),
      regex = c(
        "host.*family.*breadth",
        "host.*family.*(count|_n$)",
        "n.*host.*famil"
      ),
      required = FALSE
    )
    if (!is.null(bcol)) {
      y <- x[
        !is.na(.tip_id) & nzchar(.tip_id),
        .(
          host_family_breadth = suppressWarnings(
            max(as.numeric(get(bcol)), na.rm = TRUE)
          )
        ),
        by = .tip_id
      ]
      y[!is.finite(host_family_breadth), host_family_breadth := NA_real_]
      if (any(is.finite(y$host_family_breadth))) return(y)
    }
  }

  # Fallback: calculate the exact number of unique host families per tree tip.
  if (nrow(host_long)) {
    hl <- copy(host_long)
    idcol <- find_col(
      names(hl),
      exact = c("label", "representative_sequence_id", "tip", "tip_id", "sequence_id"),
      regex = c("^label$", "representative.*sequence", "^tip", "sequence.*id")
    )
    hfcol <- find_col(
      names(hl),
      exact = c("host_family"),
      regex = c("^host.*family$")
    )
    return(
      hl[
        !is.na(get(idcol)) &
          nzchar(as.character(get(idcol))) &
          !is.na(get(hfcol)) &
          nzchar(as.character(get(hfcol))),
        .(
          host_family_breadth = uniqueN(as.character(get(hfcol)))
        ),
        by = .( .tip_id = as.character(get(idcol)) )
      ]
    )
  }

  data.table()
}

midpoint_for_display <- function(tr) {
  # Do not add a hard package dependency only for visualization.
  # If phangorn is available, midpoint-root. Otherwise retain the
  # unrooted ML topology and log the fallback.
  if (requireNamespace("phangorn", quietly = TRUE)) {
    out <- tryCatch(
      phangorn::midpoint(tr),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      attr(out, "display_rooting") <- "midpoint"
      return(out)
    }
  }
  attr(tr, "display_rooting") <- "unrooted_fallback"
  tr
}

host_order_palette <- function(levels_vec) {
  levels_vec <- unique(levels_vec[!is.na(levels_vec) & nzchar(levels_vec)])
  if (!length(levels_vec)) return(character())

  # Fixed, publication-friendly qualitative palette; surplus categories become grey.
  base <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
    "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7",
    "#9C755F", "#BAB0AC", "#6B8E23", "#8C6D31"
  )

  vals <- rep("#BDBDBD", length(levels_vec))
  vals[seq_len(min(length(base), length(levels_vec)))] <- base[seq_len(min(length(base), length(levels_vec)))]
  stats::setNames(vals, levels_vec)
}

genus_palette <- function(levels_vec) {
  levels_vec <- unique(levels_vec[!is.na(levels_vec) & nzchar(levels_vec)])
  if (!length(levels_vec)) return(character())

  base <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
    "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7",
    "#9C755F", "#BAB0AC", "#2F4B7C", "#A05195",
    "#665191", "#D45087", "#F95D6A", "#FFA600"
  )
  vals <- rep("#9E9E9E", length(levels_vec))
  vals[seq_len(min(length(base), length(levels_vec)))] <- base[seq_len(min(length(base), length(levels_vec)))]
  stats::setNames(vals, levels_vec)
}

make_phylogeny_trait_panel <- function(
    family,
    tree_file,
    ann_master,
    host_long,
    show_genus_legend = FALSE,
    show_host_order_legend = TRUE,
    breadth_legend = TRUE
) {
  if (is.na(tree_file) || !file.exists(tree_file)) {
    stop("Tree file not found for ", family)
  }

  tr_raw <- ape::read.tree(tree_file)
  tr <- midpoint_for_display(tr_raw)
  tr <- ape::ladderize(tr, right = FALSE)
  rooting_mode <- attr(tr, "display_rooting")
  if (is.null(rooting_mode)) rooting_mode <- "unknown"

  genus_dt <- extract_tip_genus(ann_master)
  order_dt <- extract_tip_host_order(ann_master, host_long)
  breadth_dt <- extract_tip_breadth(ann_master, host_long)

  tip_anno <- data.table(.tip_id = tr$tip.label)
  if (nrow(genus_dt)) tip_anno <- merge(tip_anno, genus_dt, by = ".tip_id", all.x = TRUE)
  if (nrow(order_dt)) tip_anno <- merge(tip_anno, order_dt, by = ".tip_id", all.x = TRUE)
  if (nrow(breadth_dt)) tip_anno <- merge(tip_anno, breadth_dt, by = ".tip_id", all.x = TRUE)

  if (!"virus_genus" %in% names(tip_anno)) tip_anno[, virus_genus := NA_character_]
  if (!"host_order" %in% names(tip_anno)) tip_anno[, host_order := NA_character_]
  if (!"host_family_breadth" %in% names(tip_anno)) tip_anno[, host_family_breadth := NA_real_]

  # Limit genus legend to major genera; terminal branches still retain colours
  # for all recovered genera. "Other" is used only for palette/legend compression.
  genus_counts <- tip_anno[
    !is.na(virus_genus) & nzchar(virus_genus),
    .N,
    by = virus_genus
  ][order(-N)]

  major_genera <- head(genus_counts$virus_genus, 10)
  tip_anno[, genus_display := fifelse(
    !is.na(virus_genus) & virus_genus %in% major_genera,
    virus_genus,
    "Other"
  )]

  glevels <- c(major_genera, "Other")
  glevels <- glevels[glevels %in% unique(tip_anno$genus_display)]
  gpal <- genus_palette(glevels)

  # ggtree associates branch aesthetics with child nodes.
  # Supplying genus only for tips colours terminal branches; internal branches
  # remain grey, avoiding any implication of ancestral-state reconstruction.
  p_tree <- suppressWarnings(
    ggtree::ggtree(tr, layout = "rectangular")
  )

  td <- as.data.table(p_tree$data)
  td[, genus_display := NA_character_]
  m <- match(td$label, tip_anno$.tip_id)
  tip_rows <- which(td$isTip %in% TRUE & !is.na(m))
  td$genus_display[tip_rows] <- tip_anno$genus_display[m[tip_rows]]

  p_tree$data <- as.data.frame(td)

  p_tree <- p_tree +
    geom_tree(
      aes(colour = genus_display),
      linewidth = 0.27
    ) +
    scale_colour_manual(
      values = c(gpal, "NA" = "#6F6F6F"),
      na.value = "#6F6F6F",
      name = "Virus genus",
      guide = if (show_genus_legend) "legend" else "none"
    ) +
    labs(
      title = paste0(family, " (n = ", length(tr$tip.label), ")")
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 1)
      ),
      legend.title = element_text(size = TEXT_TINY, face = "bold"),
      legend.text = element_text(size = TEXT_TINY),
      plot.margin = margin(2, 1, 2, 1)
    )

  # y positions are inherited from the actual ggtree layout.
  tip_y <- as.data.table(p_tree$data)[
    isTip %in% TRUE,
    .( .tip_id = as.character(label), y = as.numeric(y) )
  ]
  trait <- merge(tip_y, tip_anno, by = ".tip_id", all.x = TRUE)

  # Host-order track.
  order_levels <- names(sort(table(trait$host_order), decreasing = TRUE))
  order_levels <- order_levels[!is.na(order_levels) & nzchar(order_levels)]
  opal <- host_order_palette(order_levels)

  p_order <- ggplot(trait, aes(x = 1, y = y)) +
    geom_tile(
      aes(fill = host_order),
      width = 0.95,
      height = 0.92,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = opal,
      na.value = "#E5E5E5",
      name = "Host order",
      guide = if (show_host_order_legend) "legend" else "none"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(
      ylim = range(tip_y$y, finite = TRUE) + c(-0.55, 0.55),
      expand = FALSE,
      clip = "off"
    ) +
    labs(x = NULL, y = NULL) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.title = element_text(size = TEXT_TINY, face = "bold"),
      legend.text = element_text(size = TEXT_TINY),
      plot.margin = margin(12, 0.5, 2, 0.5)
    )

  # Numeric host-family breadth track.
  breadth_vals <- trait$host_family_breadth[
    is.finite(trait$host_family_breadth) &
      trait$host_family_breadth > 0
  ]

  if (length(breadth_vals)) {
    bbreaks <- unique(round(pretty(breadth_vals, n = 4)))
    bbreaks <- bbreaks[bbreaks > 0]
    if (!length(bbreaks)) bbreaks <- sort(unique(breadth_vals))[1:min(4, length(unique(breadth_vals)))]

    p_breadth <- ggplot(trait, aes(x = 1, y = y)) +
      geom_point(
        aes(size = host_family_breadth),
        shape = 21,
        stroke = 0.25,
        fill = "black",
        colour = "black",
        na.rm = TRUE
      ) +
      scale_size_continuous(
        range = c(0.45, 3.0),
        breaks = bbreaks,
        name = "Host-family breadth",
        guide = if (breadth_legend) "legend" else "none"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(
        ylim = range(tip_y$y, finite = TRUE) + c(-0.55, 0.55),
        expand = FALSE,
        clip = "off"
      ) +
      labs(x = NULL, y = NULL) +
      theme_void(base_family = FONT_FAMILY) +
      theme(
        legend.title = element_text(size = TEXT_TINY, face = "bold"),
        legend.text = element_text(size = TEXT_TINY),
        plot.margin = margin(12, 1, 2, 0.5)
      )
  } else {
    p_breadth <- ggplot() + theme_void()
  }

  combined <- patchwork::wrap_plots(
    list(
      patchwork::wrap_elements(full = p_tree),
      patchwork::wrap_elements(full = p_order),
      patchwork::wrap_elements(full = p_breadth)
    ),
    ncol = 3,
    widths = c(4.5, 0.42, 0.65)
  )

  attr(combined, "display_rooting") <- rooting_mode
  attr(combined, "tip_n") <- length(tr$tip.label)
  attr(combined, "annotated_host_order_n") <- sum(!is.na(trait$host_order) & nzchar(trait$host_order))
  attr(combined, "annotated_breadth_n") <- sum(is.finite(trait$host_family_breadth))

  combined
}


# ============================================================
# v2.4 tree-display helpers
# ============================================================

major_host_orders <- function(trait_list, top_n = 6) {
  z <- rbindlist(
    lapply(
      trait_list,
      function(x) {
        if (!nrow(x) || !"host_order" %in% names(x)) return(data.table())
        x[
          !is.na(host_order) & nzchar(host_order),
          .(host_order)
        ]
      }
    ),
    fill = TRUE
  )

  if (!nrow(z)) return(character())

  cnt <- z[, .N, by = host_order][order(-N, host_order)]
  head(cnt$host_order, top_n)
}

prepare_family_tip_traits <- function(
    family,
    tree_file,
    ann_master,
    host_long
) {
  tr <- ape::read.tree(tree_file)

  genus_dt <- extract_tip_genus(ann_master)
  order_dt <- extract_tip_host_order(ann_master, host_long)
  breadth_dt <- extract_tip_breadth(ann_master, host_long)

  x <- data.table(.tip_id = tr$tip.label)

  if (nrow(genus_dt)) x <- merge(x, genus_dt, by = ".tip_id", all.x = TRUE)
  if (nrow(order_dt)) x <- merge(x, order_dt, by = ".tip_id", all.x = TRUE)
  if (nrow(breadth_dt)) x <- merge(x, breadth_dt, by = ".tip_id", all.x = TRUE)

  if (!"virus_genus" %in% names(x)) x[, virus_genus := NA_character_]
  if (!"host_order" %in% names(x)) x[, host_order := NA_character_]
  if (!"host_family_breadth" %in% names(x)) x[, host_family_breadth := NA_real_]

  x[, virus_family := family]
  x
}

is_monophyletic_safe <- function(tr, tips) {
  tips <- intersect(as.character(tips), tr$tip.label)
  if (length(tips) < 2L) return(FALSE)

  out <- tryCatch(
    ape::is.monophyletic(tr, tips),
    error = function(e) FALSE
  )

  isTRUE(out)
}

compress_monophyletic_genera <- function(
    tr,
    trait,
    min_clade_n = 20
) {
  if (!nrow(trait) || !"virus_genus" %in% names(trait)) {
    return(list(tree = tr, trait = trait, collapsed = data.table()))
  }

  genus_tab <- trait[
    !is.na(virus_genus) & nzchar(virus_genus),
    .(
      tips = list(.tip_id),
      n = uniqueN(.tip_id)
    ),
    by = virus_genus
  ][n >= min_clade_n]

  if (!nrow(genus_tab)) {
    return(list(tree = tr, trait = trait, collapsed = data.table()))
  }

  genus_tab[, monophyletic := vapply(
    tips,
    function(z) is_monophyletic_safe(tr, z),
    logical(1)
  )]

  genus_tab <- genus_tab[monophyletic == TRUE]
  if (!nrow(genus_tab)) {
    return(list(tree = tr, trait = trait, collapsed = data.table()))
  }

  setorder(genus_tab, -n, virus_genus)

  drop_tips <- character()
  collapsed_rows <- list()

  for (i in seq_len(nrow(genus_tab))) {
    genus <- genus_tab$virus_genus[i]
    tips <- intersect(genus_tab$tips[[i]], tr$tip.label)

    if (length(tips) < min_clade_n) next
    if (!is_monophyletic_safe(tr, tips)) next

    rep_tip <- tips[1]
    drop_tips <- c(drop_tips, setdiff(tips, rep_tip))

    tx <- trait[.tip_id %in% tips]

    dom_order <- NA_character_
    if ("host_order" %in% names(tx)) {
      zo <- tx$host_order[!is.na(tx$host_order) & nzchar(tx$host_order)]
      if (length(zo)) dom_order <- names(sort(table(zo), decreasing = TRUE))[1]
    }

    med_breadth <- NA_real_
    if ("host_family_breadth" %in% names(tx)) {
      zb <- tx$host_family_breadth[is.finite(tx$host_family_breadth)]
      if (length(zb)) med_breadth <- stats::median(zb)
    }

    collapsed_rows[[length(collapsed_rows) + 1L]] <- data.table(
      .tip_id = rep_tip,
      collapsed_genus = genus,
      collapsed_n = length(tips),
      collapsed_host_order = dom_order,
      collapsed_median_breadth = med_breadth
    )
  }

  collapsed <- if (length(collapsed_rows)) {
    rbindlist(collapsed_rows, fill = TRUE)
  } else {
    data.table()
  }

  if (!length(drop_tips)) {
    return(list(tree = tr, trait = trait, collapsed = collapsed))
  }

  tr2 <- ape::drop.tip(
    tr,
    unique(drop_tips),
    trim.internal = TRUE,
    collapse.singles = TRUE
  )

  trait2 <- trait[.tip_id %in% tr2$tip.label]

  if (nrow(collapsed)) {
    trait2 <- merge(trait2, collapsed, by = ".tip_id", all.x = TRUE)

    trait2[!is.na(collapsed_genus), virus_genus := collapsed_genus]
    trait2[!is.na(collapsed_host_order), host_order := collapsed_host_order]
    trait2[
      is.finite(collapsed_median_breadth),
      host_family_breadth := collapsed_median_breadth
    ]
  } else {
    trait2[, collapsed_genus := NA_character_]
    trait2[, collapsed_n := NA_integer_]
  }

  list(tree = tr2, trait = trait2, collapsed = collapsed)
}

count_scaled_density_df <- function(x, bins = 30, n = 512) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(data.table())

  h <- hist(x, breaks = bins, plot = FALSE)
  bin_width <- mean(diff(h$breaks))
  d <- density(x, n = n)

  data.table(
    x = d$x,
    y_count_scaled = d$y * length(x) * bin_width
  )
}

make_cladogram_trait_panel_v24 <- function(
    family,
    tree_file,
    ann_master,
    host_long,
    top_host_orders,
    host_order_colors,
    min_collapse_n = 20,
    show_order_legend = TRUE,
    show_breadth_legend = TRUE
) {
  tr_raw <- ape::read.tree(tree_file)
  tr_disp <- midpoint_for_display(tr_raw)

  rooting_mode <- attr(tr_disp, "display_rooting")
  if (is.null(rooting_mode)) rooting_mode <- "unknown"

  trait_raw <- prepare_family_tip_traits(
    family,
    tree_file,
    ann_master,
    host_long
  )

  comp <- compress_monophyletic_genera(
    tr_disp,
    trait_raw,
    min_clade_n = min_collapse_n
  )

  tr <- ape::ladderize(comp$tree, right = FALSE)
  trait <- comp$trait
  collapsed <- comp$collapsed

  trait[, host_order_display := fifelse(
    !is.na(host_order) & host_order %in% top_host_orders,
    host_order,
    fifelse(
      is.na(host_order) | !nzchar(host_order),
      "Unknown",
      "Other"
    )
  )]

  order_levels <- unique(c(top_host_orders, "Other", "Unknown"))

  pal <- host_order_colors
  if (!"Other" %in% names(pal)) pal["Other"] <- "#BDBDBD"
  if (!"Unknown" %in% names(pal)) pal["Unknown"] <- "#E5E5E5"

  p_tree <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular",
      branch.length = "none"
    )
  ) +
    geom_tree(
      linewidth = 0.32,
      colour = "#444444"
    ) +
    labs(
      title = paste0(
        family,
        " (",
        length(tr_raw$tip.label),
        " tips)"
      )
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 1)
      ),
      plot.margin = margin(2, 2, 2, 1)
    )

  tree_dt <- as.data.table(p_tree$data)
  tip_y <- tree_dt[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      y = as.numeric(y),
      x = as.numeric(x)
    )
  ]

  trait_plot <- merge(
    tip_y,
    trait,
    by = ".tip_id",
    all.x = TRUE
  )

  if (!"collapsed_genus" %in% names(trait_plot)) {
    trait_plot[, collapsed_genus := NA_character_]
  }
  if (!"collapsed_n" %in% names(trait_plot)) {
    trait_plot[, collapsed_n := NA_integer_]
  }

  collapse_points <- trait_plot[
    !is.na(collapsed_genus) & !is.na(collapsed_n)
  ]

  if (nrow(collapse_points)) {
    p_tree <- p_tree +
      geom_point(
        data = collapse_points,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        shape = 24,
        size = 2.0,
        stroke = 0.35,
        fill = "white",
        colour = "black"
      ) +
      geom_text(
        data = collapse_points,
        aes(
          x = x,
          y = y,
          label = paste0(collapsed_genus, " (n=", collapsed_n, ")")
        ),
        inherit.aes = FALSE,
        hjust = -0.05,
        vjust = 0.5,
        size = 1.8,
        family = FONT_FAMILY,
        colour = "black"
      ) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(2, 34, 2, 1))
  }

  p_order <- ggplot(
    trait_plot,
    aes(x = 1, y = y)
  ) +
    geom_tile(
      aes(fill = host_order_display),
      width = 0.94,
      height = 0.90
    ) +
    scale_fill_manual(
      values = pal,
      limits = order_levels,
      drop = TRUE,
      name = "Host order",
      guide = if (show_order_legend) "legend" else "none"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(
      ylim = range(tip_y$y, finite = TRUE) + c(-0.55, 0.55),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.title = element_text(size = TEXT_TINY, face = "bold"),
      legend.text = element_text(size = TEXT_TINY),
      legend.key.height = grid::unit(2.6, "mm"),
      plot.margin = margin(12, 0.5, 2, 0.5)
    )

  bvals <- trait_plot$host_family_breadth[
    is.finite(trait_plot$host_family_breadth) &
      trait_plot$host_family_breadth > 0
  ]

  if (length(bvals)) {
    bbreaks <- unique(round(pretty(bvals, n = 4)))
    bbreaks <- bbreaks[bbreaks > 0]

    p_breadth <- ggplot(
      trait_plot,
      aes(x = 1, y = y)
    ) +
      geom_point(
        aes(size = host_family_breadth),
        shape = 21,
        stroke = 0.22,
        fill = "black",
        colour = "black",
        na.rm = TRUE
      ) +
      scale_size_continuous(
        range = c(0.45, 2.8),
        breaks = bbreaks,
        name = "Host-family breadth",
        guide = if (show_breadth_legend) "legend" else "none"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(
        ylim = range(tip_y$y, finite = TRUE) + c(-0.55, 0.55),
        expand = FALSE,
        clip = "off"
      ) +
      theme_void(base_family = FONT_FAMILY) +
      theme(
        legend.title = element_text(size = TEXT_TINY, face = "bold"),
        legend.text = element_text(size = TEXT_TINY),
        legend.key.height = grid::unit(3.0, "mm"),
        plot.margin = margin(12, 1, 2, 0.5)
      )
  } else {
    p_breadth <- ggplot() + theme_void()
  }

  combined <- patchwork::wrap_plots(
    list(
      patchwork::wrap_elements(full = p_tree),
      patchwork::wrap_elements(full = p_order),
      patchwork::wrap_elements(full = p_breadth)
    ),
    ncol = 3,
    widths = c(5.7, 0.45, 0.80)
  )

  attr(combined, "display_rooting") <- rooting_mode
  attr(combined, "original_tip_n") <- length(tr_raw$tip.label)
  attr(combined, "display_tip_n") <- length(tr$tip.label)
  attr(combined, "collapsed_genus_n") <- nrow(collapsed)

  combined
}


# ============================================================
# v2.5 helpers
# ============================================================

parse_iqtree_support <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    return(c(sh_alrt = NA_real_, ufboot = NA_real_))
  }

  m <- regmatches(
    x,
    gregexpr("[0-9]+(?:\\.[0-9]+)?", x, perl = TRUE)
  )[[1]]

  if (!length(m)) {
    return(c(sh_alrt = NA_real_, ufboot = NA_real_))
  }

  z <- suppressWarnings(as.numeric(m))
  z <- z[is.finite(z)]

  if (!length(z)) {
    return(c(sh_alrt = NA_real_, ufboot = NA_real_))
  }

  if (length(z) >= 2L) {
    return(c(sh_alrt = z[1], ufboot = z[2]))
  }

  c(sh_alrt = NA_real_, ufboot = z[1])
}

node_descendant_tips <- function(tr, node) {
  out <- tryCatch(
    ape::extract.clade(tr, node = node)$tip.label,
    error = function(e) character()
  )
  as.character(out)
}

build_node_support_table <- function(tr) {
  ntip <- length(tr$tip.label)
  nodes <- seq.int(ntip + 1L, ntip + tr$Nnode)

  labels <- tr$node.label
  if (is.null(labels)) labels <- rep(NA_character_, tr$Nnode)
  if (length(labels) < tr$Nnode) {
    labels <- c(labels, rep(NA_character_, tr$Nnode - length(labels)))
  }

  sup <- t(vapply(
    labels,
    parse_iqtree_support,
    numeric(2)
  ))

  dt <- data.table(
    node = nodes,
    node_label = labels,
    sh_alrt = sup[, "sh_alrt"],
    ufboot = sup[, "ufboot"]
  )

  dt[, tips := lapply(node, function(z) node_descendant_tips(tr, z))]
  dt[, clade_n := lengths(tips)]

  dt
}

supported_fold_candidates <- function(
    tr,
    min_clade_n = 18,
    max_clade_fraction = 0.18,
    sh_alrt_min = 80,
    ufboot_min = 95
) {
  dt <- build_node_support_table(tr)
  ntip <- length(tr$tip.label)

  dt <- dt[
    clade_n >= min_clade_n &
      clade_n <= floor(ntip * max_clade_fraction)
  ]

  if (!nrow(dt)) return(dt)

  dt[, support_available := is.finite(sh_alrt) | is.finite(ufboot)]

  dt[, support_ok := fifelse(
    is.finite(sh_alrt) & is.finite(ufboot),
    sh_alrt >= sh_alrt_min & ufboot >= ufboot_min,
    fifelse(
      is.finite(ufboot),
      ufboot >= ufboot_min,
      FALSE
    )
  )]

  dt[support_ok == TRUE]
}

select_nonoverlapping_folds <- function(
    tr,
    candidates,
    target_display_n
) {
  ntip <- length(tr$tip.label)

  if (!nrow(candidates) || target_display_n >= ntip) {
    return(data.table())
  }

  # Prefer larger well-supported clades while preserving major topology.
  cand <- copy(candidates)
  cand[, support_score := fifelse(
    is.finite(sh_alrt) & is.finite(ufboot),
    sh_alrt + ufboot,
    fifelse(is.finite(ufboot), ufboot, 0)
  )]

  setorder(cand, -clade_n, -support_score, node)

  selected <- list()
  occupied <- character()
  display_n <- ntip

  for (i in seq_len(nrow(cand))) {
    tips <- cand$tips[[i]]

    if (!length(tips)) next
    if (any(tips %in% occupied)) next

    reduction <- length(tips) - 1L
    if (reduction <= 0L) next

    # Avoid collapsing so aggressively that the display falls far below target.
    if (display_n - reduction < target_display_n * 0.80) next

    selected[[length(selected) + 1L]] <- cand[i]
    occupied <- c(occupied, tips)
    display_n <- display_n - reduction

    if (display_n <= target_display_n) break
  }

  if (!length(selected)) return(data.table())

  rbindlist(selected, fill = TRUE)
}

choose_farthest_representative <- function(tr, tips) {
  tips <- intersect(as.character(tips), tr$tip.label)
  if (!length(tips)) return(NA_character_)
  if (length(tips) == 1L) return(tips)

  depths <- tryCatch(
    ape::node.depth.edgelength(tr),
    error = function(e) NULL
  )

  if (is.null(depths)) return(sort(tips)[1])

  idx <- match(tips, tr$tip.label)
  idx <- idx[is.finite(idx)]
  if (!length(idx)) return(sort(tips)[1])

  tips[which.max(depths[idx])]
}

aggregate_fold_traits <- function(trait, tips) {
  x <- trait[.tip_id %in% tips]

  dom_order <- NA_character_
  z <- x$host_order[
    !is.na(x$host_order) &
      nzchar(x$host_order)
  ]
  if (length(z)) {
    dom_order <- names(sort(table(z), decreasing = TRUE))[1]
  }

  med_breadth <- NA_real_
  b <- x$host_family_breadth[
    is.finite(x$host_family_breadth)
  ]
  if (length(b)) {
    med_breadth <- stats::median(b)
  }

  list(
    dominant_host_order = dom_order,
    median_host_family_breadth = med_breadth
  )
}

make_supported_folded_display <- function(
    tr,
    trait,
    target_display_n,
    min_clade_n = 18,
    max_clade_fraction = 0.18,
    sh_alrt_min = 80,
    ufboot_min = 95
) {
  candidates <- supported_fold_candidates(
    tr,
    min_clade_n = min_clade_n,
    max_clade_fraction = max_clade_fraction,
    sh_alrt_min = sh_alrt_min,
    ufboot_min = ufboot_min
  )

  support_tab <- build_node_support_table(tr)
  support_present_n <- support_tab[
    is.finite(sh_alrt) | is.finite(ufboot),
    .N
  ]

  selected <- select_nonoverlapping_folds(
    tr,
    candidates,
    target_display_n = target_display_n
  )

  if (!nrow(selected)) {
    trait2 <- copy(trait)
    trait2[, collapsed := FALSE]
    trait2[, collapsed_n := 1L]

    return(list(
      tree = tr,
      trait = trait2,
      selected = data.table(),
      support_present_n = support_present_n,
      display_tip_n = length(tr$tip.label)
    ))
  }

  depth <- tryCatch(
    ape::node.depth.edgelength(tr),
    error = function(e) NULL
  )

  collapsed_rows <- list()
  drop_tips <- character()

  for (i in seq_len(nrow(selected))) {
    tips <- selected$tips[[i]]
    rep_tip <- choose_farthest_representative(tr, tips)
    if (is.na(rep_tip)) next

    ag <- aggregate_fold_traits(trait, tips)

    drop_tips <- c(
      drop_tips,
      setdiff(tips, rep_tip)
    )

    collapsed_rows[[length(collapsed_rows) + 1L]] <- data.table(
      .tip_id = rep_tip,
      collapsed = TRUE,
      collapsed_n = length(tips),
      host_order_fold = ag$dominant_host_order,
      host_family_breadth_fold = ag$median_host_family_breadth,
      fold_node = selected$node[i],
      fold_sh_alrt = selected$sh_alrt[i],
      fold_ufboot = selected$ufboot[i]
    )
  }

  collapsed <- if (length(collapsed_rows)) {
    rbindlist(collapsed_rows, fill = TRUE)
  } else {
    data.table()
  }

  if (!length(drop_tips)) {
    trait2 <- copy(trait)
    trait2[, collapsed := FALSE]
    trait2[, collapsed_n := 1L]

    return(list(
      tree = tr,
      trait = trait2,
      selected = selected,
      support_present_n = support_present_n,
      display_tip_n = length(tr$tip.label)
    ))
  }

  tr2 <- ape::drop.tip(
    tr,
    unique(drop_tips),
    trim.internal = TRUE,
    collapse.singles = TRUE
  )

  trait2 <- trait[.tip_id %in% tr2$tip.label]
  trait2[, collapsed := FALSE]
  trait2[, collapsed_n := 1L]

  if (nrow(collapsed)) {
    trait2 <- merge(
      trait2,
      collapsed,
      by = ".tip_id",
      all.x = TRUE,
      suffixes = c("", ".fold")
    )

    trait2[
      is.na(collapsed.fold),
      collapsed.fold := FALSE
    ]

    trait2[
      collapsed.fold == TRUE,
      `:=`(
        collapsed = TRUE,
        collapsed_n = collapsed_n.fold,
        host_order = host_order_fold,
        host_family_breadth = host_family_breadth_fold
      )
    ]

    drop_cols <- intersect(
      c(
        "collapsed.fold",
        "collapsed_n.fold",
        "host_order_fold",
        "host_family_breadth_fold"
      ),
      names(trait2)
    )

    if (length(drop_cols)) {
      trait2[, (drop_cols) := NULL]
    }
  }

  list(
    tree = tr2,
    trait = trait2,
    selected = selected,
    support_present_n = support_present_n,
    display_tip_n = length(tr2$tip.label)
  )
}

shared_breadth_breaks <- function(trait_list, n = 4) {
  z <- unlist(lapply(
    trait_list,
    function(x) {
      if (!nrow(x) || !"host_family_breadth" %in% names(x)) return(numeric())
      x$host_family_breadth[
        is.finite(x$host_family_breadth) &
          x$host_family_breadth > 0
      ]
    }
  ))

  if (!length(z)) return(NULL)

  br <- unique(round(pretty(z, n = n)))
  br <- br[br > 0 & br <= max(z, na.rm = TRUE)]

  if (!length(br)) {
    br <- sort(unique(z))[seq_len(min(4, length(unique(z))))]
  }

  br
}

main_host_order_palette <- function(orders) {
  orders <- unique(orders[!is.na(orders) & nzchar(orders)])

  base <- c(
    "#4E79A7",
    "#F28E2B",
    "#59A14F",
    "#E15759",
    "#B07AA1",
    "#76B7B2",
    "#EDC948",
    "#9C755F"
  )

  vals <- base[seq_len(min(length(base), length(orders)))]
  stats::setNames(vals, orders[seq_along(vals)])
}

prepare_order_display <- function(
    trait,
    top_orders
) {
  x <- copy(trait)

  x[, host_order_display := fifelse(
    !is.na(host_order) &
      nzchar(host_order) &
      host_order %in% top_orders,
    host_order,
    fifelse(
      is.na(host_order) | !nzchar(host_order),
      NA_character_,
      "Other"
    )
  )]

  x
}

make_endpoint_phylogram <- function(
    family,
    tr,
    trait,
    top_host_orders,
    order_palette,
    breadth_breaks,
    show_legend = TRUE,
    folded = FALSE,
    label_large_folds_n = 45
) {
  trait <- prepare_order_display(
    trait,
    top_host_orders
  )

  pal <- c(
    order_palette,
    Other = "#BDBDBD"
  )

  p <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular"
    )
  ) +
    geom_tree(
      linewidth = 0.30,
      colour = "#444444"
    ) +
    labs(
      title = paste0(
        family,
        " (n = ",
        if (folded && "collapsed_n" %in% names(trait)) {
          sum(trait$collapsed_n, na.rm = TRUE)
        } else {
          length(tr$tip.label)
        },
        ")"
      )
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 1)
      ),
      plot.margin = margin(2, 16, 2, 2)
    )

  td <- as.data.table(p$data)[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      x = as.numeric(x),
      y = as.numeric(y)
    )
  ]

  td <- merge(
    td,
    trait,
    by = ".tip_id",
    all.x = TRUE
  )

  if (!"collapsed" %in% names(td)) td[, collapsed := FALSE]
  if (!"collapsed_n" %in% names(td)) td[, collapsed_n := 1L]

  # Endpoint points: colour = dominant host order, size = host-family breadth.
  p <- p +
    geom_point(
      data = td,
      aes(
        x = x,
        y = y,
        fill = host_order_display,
        size = host_family_breadth,
        shape = collapsed
      ),
      inherit.aes = FALSE,
      colour = "black",
      stroke = 0.22,
      alpha = 0.95,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = pal,
      na.value = "#E6E6E6",
      breaks = c(top_host_orders, "Other"),
      name = "Dominant host order",
      guide = if (show_legend) "legend" else "none"
    ) +
    scale_size_continuous(
      range = c(0.65, 3.0),
      breaks = breadth_breaks,
      name = "Host-family breadth",
      guide = if (show_legend) "legend" else "none"
    ) +
    scale_shape_manual(
      values = c(
        `FALSE` = 21,
        `TRUE` = 24
      ),
      guide = "none"
    ) +
    coord_cartesian(clip = "off")

  if (folded) {
    fold_lab <- td[
      collapsed %in% TRUE &
        is.finite(collapsed_n) &
        collapsed_n >= label_large_folds_n
    ]

    if (nrow(fold_lab)) {
      p <- p +
        geom_text(
          data = fold_lab,
          aes(
            x = x,
            y = y,
            label = paste0("n=", collapsed_n)
          ),
          inherit.aes = FALSE,
          hjust = -0.30,
          vjust = 0.5,
          family = FONT_FAMILY,
          size = 1.75,
          colour = "#555555"
        )
    }
  }

  p +
    guides(
      fill = guide_legend(
        order = 1,
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = 21,
          size = 2.4,
          colour = "black"
        )
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        nrow = 1
      )
    ) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(
        size = TEXT_TINY,
        face = "bold"
      ),
      legend.text = element_text(
        size = TEXT_TINY
      ),
      legend.key.height = grid::unit(3.1, "mm"),
      legend.key.width = grid::unit(3.6, "mm")
    )
}

count_scaled_density_explicit <- function(
    x,
    bins = 30,
    n = 512
) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(list(
    df = data.table(),
    binwidth = NA_real_,
    boundary = NA_real_
  ))

  xr <- range(x)
  binwidth <- diff(xr) / bins
  if (!is.finite(binwidth) || binwidth <= 0) {
    return(list(
      df = data.table(),
      binwidth = NA_real_,
      boundary = NA_real_
    ))
  }

  boundary <- xr[1]

  d <- density(
    x,
    from = xr[1],
    to = xr[2],
    n = n
  )

  df <- data.table(
    x = d$x,
    y_count_scaled = d$y * length(x) * binwidth
  )

  list(
    df = df,
    binwidth = binwidth,
    boundary = boundary
  )
}


# ============================================================
# v2.6 helpers
# ============================================================

prepare_endpoint_trait_display_v26 <- function(
    trait,
    top_orders
) {
  x <- copy(trait)

  x[, host_order_display := fifelse(
    !is.na(host_order) &
      nzchar(host_order) &
      host_order %in% top_orders,
    host_order,
    fifelse(
      is.na(host_order) | !nzchar(host_order),
      NA_character_,
      "Other"
    )
  )]

  # Log1p visual scaling makes breadth 1/2/3/4/5 clearly distinguishable
  # while retaining continuous values above 5 without hard categorisation.
  x[, breadth_size_value := fifelse(
    is.finite(host_family_breadth) &
      host_family_breadth > 0,
    log1p(host_family_breadth),
    NA_real_
  )]

  x
}

make_fold_triangle_df <- function(
    tree_data,
    trait_data
) {
  if (!nrow(trait_data)) return(data.table())

  td <- as.data.table(tree_data)

  tips <- td[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      node,
      parent,
      x,
      y
    )
  ]

  parent_x <- td[, .(
    parent = node,
    x_parent = x
  )]

  tips <- merge(
    tips,
    parent_x,
    by = "parent",
    all.x = TRUE
  )

  z <- merge(
    tips,
    trait_data[
      collapsed %in% TRUE,
      .(
        .tip_id,
        collapsed_n,
        host_order_display,
        host_family_breadth,
        breadth_size_value
      )
    ],
    by = ".tip_id"
  )

  if (!nrow(z)) return(data.table())

  z[
    !is.finite(x_parent),
    x_parent := x
  ]

  z[, edge_len := pmax(x - x_parent, 0)]
  z[, x_base := x_parent + 0.62 * edge_len]

  # Keep triangles readable but not so tall that they obscure neighbouring clades.
  z[, half_height := pmin(
    4.0,
    pmax(
      0.80,
      0.95 * log10(pmax(collapsed_n, 2))
    )
  )]

  # Three vertices per collapsed display unit:
  # vertical base + right-pointing tip.
  out <- rbindlist(
    lapply(
      seq_len(nrow(z)),
      function(i) {
        data.table(
          fold_id = i,
          .tip_id = z$.tip_id[i],
          collapsed_n = z$collapsed_n[i],
          host_order_display = z$host_order_display[i],
          host_family_breadth = z$host_family_breadth[i],
          breadth_size_value = z$breadth_size_value[i],
          x_poly = c(
            z$x_base[i],
            z$x_base[i],
            z$x[i]
          ),
          y_poly = c(
            z$y[i] - z$half_height[i],
            z$y[i] + z$half_height[i],
            z$y[i]
          ),
          x_tip = z$x[i],
          y_tip = z$y[i]
        )
      }
    )
  )

  out
}

make_endpoint_phylogram_v26 <- function(
    family,
    tr,
    trait,
    top_host_orders,
    order_palette,
    show_legend = TRUE,
    folded = FALSE,
    breadth_legend_breaks = 1:5,
    label_fold_n_min = 35
) {
  trait <- prepare_endpoint_trait_display_v26(
    trait,
    top_host_orders
  )

  pal <- c(
    order_palette,
    Other = "#BDBDBD"
  )

  # First build the tree to recover branch/node coordinates.
  p0 <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular"
    )
  )

  td <- as.data.table(p0$data)

  endpoint <- trait[, .(
    label = .tip_id,
    branch_order = host_order_display,
    host_order_display,
    host_family_breadth,
    breadth_size_value,
    collapsed,
    collapsed_n
  )]

  td <- merge(
    td,
    endpoint,
    by = "label",
    all.x = TRUE,
    sort = FALSE
  )

  # Internal branches have NA branch_order and remain dark grey.
  # Only terminal edges inherit the dominant-host-order colour.
  p0$data <- as.data.frame(td)
  p0$layers <- list()

  p <- p0 +
    geom_tree(
      aes(colour = branch_order),
      linewidth = 0.33
    ) +
    scale_colour_manual(
      values = pal,
      breaks = c(top_host_orders, "Other"),
      na.value = "#444444",
      name = "Dominant host order",
      guide = if (show_legend) "legend" else "none"
    ) +
    labs(
      title = paste0(
        family,
        " (n = ",
        if (folded && "collapsed_n" %in% names(trait)) {
          sum(trait$collapsed_n, na.rm = TRUE)
        } else {
          length(tr$tip.label)
        },
        ")"
      )
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 1)
      ),
      plot.margin = margin(2, 26, 2, 2)
    )

  tree_tip_xy <- as.data.table(p$data)[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      x = as.numeric(x),
      y = as.numeric(y)
    )
  ]

  tip_display <- merge(
    tree_tip_xy,
    trait,
    by = ".tip_id",
    all.x = TRUE
  )

  # Actual host-family breadth remains continuous.
  # log1p is used only as a visual size transform so 1-5 are visibly distinct.
  size_breaks_plot <- log1p(breadth_legend_breaks)

  p <- p +
    geom_point(
      data = tip_display,
      aes(
        x = x,
        y = y,
        size = breadth_size_value,
        colour = host_order_display
      ),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      stroke = 0.45,
      na.rm = TRUE
    ) +
    scale_size_continuous(
      range = c(1.0, 4.8),
      breaks = size_breaks_plot,
      labels = as.character(breadth_legend_breaks),
      name = "Host-family breadth",
      guide = if (show_legend) "legend" else "none"
    )

  # True horizontal triangular wedges for folded supported clades.
  if (folded) {
    tri <- make_fold_triangle_df(
      p$data,
      tip_display
    )

    if (nrow(tri)) {
      p <- p +
        geom_polygon(
          data = tri,
          aes(
            x = x_poly,
            y = y_poly,
            group = fold_id
          ),
          inherit.aes = FALSE,
          fill = "white",
          colour = "#555555",
          linewidth = 0.38
        )

      fold_label <- unique(
        tri[
          collapsed_n >= label_fold_n_min,
          .(
            fold_id,
            collapsed_n,
            x_tip,
            y_tip
          )
        ]
      )

      if (nrow(fold_label)) {
        p <- p +
          geom_text(
            data = fold_label,
            aes(
              x = x_tip,
              y = y_tip,
              label = paste0("n=", collapsed_n)
            ),
            inherit.aes = FALSE,
            hjust = -0.22,
            vjust = 0.5,
            family = FONT_FAMILY,
            size = 1.8,
            colour = "#555555"
          )
      }
    }
  }

  p +
    coord_cartesian(clip = "off") +
    guides(
      colour = guide_legend(
        order = 1,
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.2,
          shape = NA
        )
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        nrow = 1,
        override.aes = list(
          colour = "black",
          fill = "white",
          shape = 21,
          stroke = 0.45
        )
      )
    ) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(
        size = TEXT_TINY,
        face = "bold"
      ),
      legend.text = element_text(
        size = TEXT_TINY
      ),
      legend.key.height = grid::unit(3.3, "mm"),
      legend.key.width = grid::unit(4.1, "mm")
    )
}


# ============================================================
# v2.7 helpers
# ============================================================

topology_fold_candidates_v27 <- function(
    tr,
    min_clade_n = 4,
    max_clade_fraction = 0.30
) {
  ntip <- length(tr$tip.label)
  ntot <- ntip + tr$Nnode

  nodes <- seq.int(ntip + 1L, ntot)

  out <- rbindlist(
    lapply(
      nodes,
      function(node) {
        tips <- node_descendant_tips(tr, node)
        n <- length(tips)

        if (
          n < min_clade_n ||
          n >= ntip ||
          n > floor(ntip * max_clade_fraction)
        ) {
          return(NULL)
        }

        data.table(
          node = node,
          clade_n = n,
          tips = list(tips)
        )
      }
    ),
    fill = TRUE
  )

  if (!nrow(out)) return(data.table())

  setorder(out, -clade_n, node)
  out
}

select_topology_folds_v27 <- function(
    tr,
    target_display_n = 330,
    min_clade_n = 4,
    max_clade_fraction = 0.30,
    lower_target_fraction = 0.88
) {
  ntip <- length(tr$tip.label)

  if (target_display_n >= ntip) {
    return(data.table())
  }

  cand <- topology_fold_candidates_v27(
    tr,
    min_clade_n = min_clade_n,
    max_clade_fraction = max_clade_fraction
  )

  if (!nrow(cand)) return(data.table())

  selected <- list()
  occupied <- character()
  display_n <- ntip
  lower_target <- target_display_n * lower_target_fraction

  # Pass 1: largest non-overlapping subtrees first.
  for (i in seq_len(nrow(cand))) {
    tips <- cand$tips[[i]]

    if (!length(tips)) next
    if (any(tips %in% occupied)) next

    reduction <- length(tips) - 1L
    if (reduction <= 0L) next

    proposed <- display_n - reduction

    # Avoid a large overshoot below the desired display size.
    if (proposed < lower_target) next

    selected[[length(selected) + 1L]] <- cand[i]
    occupied <- c(occupied, tips)
    display_n <- proposed

    if (display_n <= target_display_n) break
  }

  # Pass 2: if still above target, add smaller compatible clades
  # ordered by how close they move the display toward target.
  if (display_n > target_display_n) {
    remain <- cand[
      !vapply(
        tips,
        function(z) any(z %in% occupied),
        logical(1)
      )
    ]

    if (nrow(remain)) {
      remain[, reduction := clade_n - 1L]
      remain[, distance_after := abs(
        (display_n - reduction) - target_display_n
      )]

      setorder(
        remain,
        distance_after,
        -clade_n,
        node
      )

      for (i in seq_len(nrow(remain))) {
        tips <- remain$tips[[i]]

        if (any(tips %in% occupied)) next

        reduction <- length(tips) - 1L
        proposed <- display_n - reduction

        if (proposed < lower_target) next

        selected[[length(selected) + 1L]] <- remain[i]
        occupied <- c(occupied, tips)
        display_n <- proposed

        if (display_n <= target_display_n) break
      }
    }
  }

  if (!length(selected)) return(data.table())

  ans <- rbindlist(selected, fill = TRUE)
  ans[, display_n_after_selection := display_n]
  ans
}

make_topology_folded_display_v27 <- function(
    tr,
    trait,
    target_display_n = 330,
    min_clade_n = 4,
    max_clade_fraction = 0.30
) {
  selected <- select_topology_folds_v27(
    tr,
    target_display_n = target_display_n,
    min_clade_n = min_clade_n,
    max_clade_fraction = max_clade_fraction
  )

  if (!nrow(selected)) {
    trait2 <- copy(trait)
    trait2[, collapsed := FALSE]
    trait2[, collapsed_n := 1L]

    return(list(
      tree = tr,
      trait = trait2,
      selected = data.table(),
      display_tip_n = length(tr$tip.label)
    ))
  }

  collapsed_rows <- list()
  drop_tips <- character()

  for (i in seq_len(nrow(selected))) {
    tips <- selected$tips[[i]]
    rep_tip <- choose_farthest_representative(
      tr,
      tips
    )

    if (is.na(rep_tip)) next

    ag <- aggregate_fold_traits(
      trait,
      tips
    )

    drop_tips <- c(
      drop_tips,
      setdiff(tips, rep_tip)
    )

    collapsed_rows[[length(collapsed_rows) + 1L]] <- data.table(
      .tip_id = rep_tip,
      collapsed = TRUE,
      collapsed_n = length(tips),
      host_order_fold = ag$dominant_host_order,
      host_family_breadth_fold = ag$median_host_family_breadth,
      fold_node = selected$node[i]
    )
  }

  collapsed <- if (length(collapsed_rows)) {
    rbindlist(collapsed_rows, fill = TRUE)
  } else {
    data.table()
  }

  if (!length(drop_tips)) {
    trait2 <- copy(trait)
    trait2[, collapsed := FALSE]
    trait2[, collapsed_n := 1L]

    return(list(
      tree = tr,
      trait = trait2,
      selected = selected,
      display_tip_n = length(tr$tip.label)
    ))
  }

  tr2 <- ape::drop.tip(
    tr,
    unique(drop_tips),
    trim.internal = TRUE,
    collapse.singles = TRUE
  )

  trait2 <- trait[
    .tip_id %in% tr2$tip.label
  ]
  trait2[, collapsed := FALSE]
  trait2[, collapsed_n := 1L]

  if (nrow(collapsed)) {
    trait2 <- merge(
      trait2,
      collapsed,
      by = ".tip_id",
      all.x = TRUE,
      suffixes = c("", ".fold")
    )

    trait2[
      is.na(collapsed.fold),
      collapsed.fold := FALSE
    ]

    trait2[
      collapsed.fold == TRUE,
      `:=`(
        collapsed = TRUE,
        collapsed_n = collapsed_n.fold,
        host_order = host_order_fold,
        host_family_breadth = host_family_breadth_fold
      )
    ]

    drop_cols <- intersect(
      c(
        "collapsed.fold",
        "collapsed_n.fold",
        "host_order_fold",
        "host_family_breadth_fold"
      ),
      names(trait2)
    )

    if (length(drop_cols)) {
      trait2[, (drop_cols) := NULL]
    }
  }

  list(
    tree = tr2,
    trait = trait2,
    selected = selected,
    display_tip_n = length(tr2$tip.label)
  )
}

make_fold_triangle_df_v27 <- function(
    tree_data,
    trait_data
) {
  if (!nrow(trait_data)) return(data.table())

  td <- as.data.table(tree_data)

  tips <- td[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      node,
      parent,
      x,
      y
    )
  ]

  parent_x <- td[, .(
    parent = node,
    x_parent = x
  )]

  tips <- merge(
    tips,
    parent_x,
    by = "parent",
    all.x = TRUE
  )

  z <- merge(
    tips,
    trait_data[
      collapsed %in% TRUE,
      .(
        .tip_id,
        collapsed_n,
        host_order_display,
        host_family_breadth,
        breadth_size_value
      )
    ],
    by = ".tip_id"
  )

  if (!nrow(z)) return(data.table())

  z[
    !is.finite(x_parent),
    x_parent := x
  ]

  z[, edge_len := pmax(x - x_parent, 0)]

  # Make the triangle visually read as a collapsed clade rather than a special tip.
  z[, x_base := x_parent + 0.40 * edge_len]

  z[, half_height := pmin(
    7.0,
    pmax(
      1.15,
      sqrt(pmax(collapsed_n, 2)) / 2.25
    )
  )]

  out <- rbindlist(
    lapply(
      seq_len(nrow(z)),
      function(i) {
        data.table(
          fold_id = i,
          .tip_id = z$.tip_id[i],
          collapsed_n = z$collapsed_n[i],
          host_order_display = z$host_order_display[i],
          host_family_breadth = z$host_family_breadth[i],
          breadth_size_value = z$breadth_size_value[i],
          x_poly = c(
            z$x_base[i],
            z$x_base[i],
            z$x[i]
          ),
          y_poly = c(
            z$y[i] - z$half_height[i],
            z$y[i] + z$half_height[i],
            z$y[i]
          ),
          x_tip = z$x[i],
          y_tip = z$y[i]
        )
      }
    )
  )

  out
}

make_endpoint_phylogram_v27 <- function(
    family,
    tr,
    trait,
    top_host_orders,
    order_palette,
    show_legend = TRUE,
    folded = FALSE,
    breadth_legend_breaks = 1:5,
    label_top_fold_n = 5
) {
  trait <- prepare_endpoint_trait_display_v26(
    trait,
    top_host_orders
  )

  pal <- c(
    order_palette,
    Other = "#BDBDBD"
  )

  p0 <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular"
    )
  )

  td <- as.data.table(p0$data)

  endpoint <- trait[, .(
    label = .tip_id,
    branch_order = host_order_display,
    host_order_display,
    host_family_breadth,
    breadth_size_value,
    collapsed,
    collapsed_n
  )]

  # Preserve the exact ggtree row order.
  m <- match(
    as.character(td$label),
    endpoint$label
  )

  td[, branch_order := NA_character_]
  td[, host_order_display := NA_character_]
  td[, host_family_breadth := NA_real_]
  td[, breadth_size_value := NA_real_]
  td[, collapsed := FALSE]
  td[, collapsed_n := 1L]

  hit <- which(
    !is.na(m) &
      td$isTip %in% TRUE
  )

  if (length(hit)) {
    td$branch_order[hit] <- endpoint$branch_order[m[hit]]
    td$host_order_display[hit] <- endpoint$host_order_display[m[hit]]
    td$host_family_breadth[hit] <- endpoint$host_family_breadth[m[hit]]
    td$breadth_size_value[hit] <- endpoint$breadth_size_value[m[hit]]
    td$collapsed[hit] <- endpoint$collapsed[m[hit]]
    td$collapsed_n[hit] <- endpoint$collapsed_n[m[hit]]
  }

  p0$data <- as.data.frame(td)
  p0$layers <- list()

  p <- p0 +
    geom_tree(
      aes(colour = branch_order),
      linewidth = 0.34
    ) +
    scale_colour_manual(
      values = pal,
      breaks = c(top_host_orders, "Other"),
      na.value = "#444444",
      name = "Dominant host order",
      guide = if (show_legend) "legend" else "none"
    ) +
    labs(
      title = paste0(
        family,
        " (n = ",
        if (folded && "collapsed_n" %in% names(trait)) {
          sum(trait$collapsed_n, na.rm = TRUE)
        } else {
          length(tr$tip.label)
        },
        ")"
      )
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 2)
      ),
      plot.margin = margin(2, 28, 2, 2)
    )

  tip_display <- as.data.table(p$data)[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      x = as.numeric(x),
      y = as.numeric(y),
      host_order_display,
      host_family_breadth,
      breadth_size_value,
      collapsed,
      collapsed_n
    )
  ]

  # Stronger point-size separation for breadth 1-5.
  size_breaks_plot <- log1p(breadth_legend_breaks)

  regular_tip <- tip_display[
    collapsed %in% FALSE |
      is.na(collapsed)
  ]

  if (nrow(regular_tip)) {
    p <- p +
      geom_point(
        data = regular_tip,
        aes(
          x = x,
          y = y,
          size = breadth_size_value,
          colour = host_order_display
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        stroke = 0.46,
        na.rm = TRUE
      )
  }

  p <- p +
    scale_size_continuous(
      range = c(1.15, 5.6),
      breaks = size_breaks_plot,
      labels = as.character(
        breadth_legend_breaks
      ),
      name = "Host-family breadth",
      guide = if (show_legend) "legend" else "none"
    )

  if (folded) {
    tri <- make_fold_triangle_df_v27(
      p$data,
      tip_display
    )

    if (nrow(tri)) {
      p <- p +
        geom_polygon(
          data = tri,
          aes(
            x = x_poly,
            y = y_poly,
            group = fold_id
          ),
          inherit.aes = FALSE,
          fill = "white",
          colour = "#555555",
          linewidth = 0.42
        )

      # Put the median breadth point at the tip of each collapsed triangle.
      fold_point <- unique(
        tri[, .(
          fold_id,
          x_tip,
          y_tip,
          host_order_display,
          breadth_size_value,
          collapsed_n
        )]
      )

      p <- p +
        geom_point(
          data = fold_point,
          aes(
            x = x_tip,
            y = y_tip,
            size = breadth_size_value,
            colour = host_order_display
          ),
          inherit.aes = FALSE,
          shape = 21,
          fill = "white",
          stroke = 0.50,
          na.rm = TRUE
        )

      # Label only the largest few folded clades.
      fold_label <- copy(fold_point)
      setorder(
        fold_label,
        -collapsed_n
      )

      fold_label <- head(
        fold_label,
        label_top_fold_n
      )

      if (nrow(fold_label)) {
        p <- p +
          geom_text(
            data = fold_label,
            aes(
              x = x_tip,
              y = y_tip,
              label = paste0(
                "n=",
                collapsed_n
              )
            ),
            inherit.aes = FALSE,
            hjust = -0.24,
            vjust = 0.5,
            family = FONT_FAMILY,
            size = 1.85,
            colour = "#555555"
          )
      }
    }
  }

  p +
    coord_cartesian(
      clip = "off"
    ) +
    guides(
      colour = guide_legend(
        order = 1,
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.5,
          shape = NA
        )
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        nrow = 1,
        override.aes = list(
          colour = "black",
          fill = "white",
          shape = 21,
          stroke = 0.45
        )
      )
    ) +
    theme(
      # Top-left legend block, outside the tree body.
      legend.position = if (show_legend) "top" else "none",
      legend.justification = "left",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(
        size = TEXT_TINY,
        face = "bold",
        hjust = 0
      ),
      legend.text = element_text(
        size = TEXT_TINY
      ),
      legend.key.height = grid::unit(
        3.4,
        "mm"
      ),
      legend.key.width = grid::unit(
        4.3,
        "mm"
      )
    )
}


# ============================================================
# v2.8 rendering helpers
# ============================================================

make_fold_triangle_df_v28 <- function(
    tree_data,
    trait_data
) {
  if (!nrow(trait_data)) return(data.table())

  td <- as.data.table(tree_data)

  tips <- td[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      node,
      parent,
      x,
      y
    )
  ]

  parent_x <- td[, .(
    parent = node,
    x_parent = x
  )]

  tips <- merge(
    tips,
    parent_x,
    by = "parent",
    all.x = TRUE
  )

  z <- merge(
    tips,
    trait_data[
      collapsed %in% TRUE,
      .(
        .tip_id,
        collapsed_n,
        host_order_display,
        host_family_breadth,
        breadth_size_value
      )
    ],
    by = ".tip_id"
  )

  if (!nrow(z)) return(data.table())

  z[
    !is.finite(x_parent),
    x_parent := x
  ]

  z[, edge_len := pmax(x - x_parent, 0)]

  # v2.8: shorter triangles and much smaller vertical span.
  z[, x_base := x_parent + 0.60 * edge_len]

  z[, half_height := pmin(
    2.3,
    pmax(
      0.70,
      0.35 + 0.55 * log10(
        pmax(collapsed_n, 2)
      )
    )
  )]

  rbindlist(
    lapply(
      seq_len(nrow(z)),
      function(i) {
        data.table(
          fold_id = i,
          .tip_id = z$.tip_id[i],
          collapsed_n = z$collapsed_n[i],
          host_order_display = z$host_order_display[i],
          host_family_breadth = z$host_family_breadth[i],
          breadth_size_value = z$breadth_size_value[i],
          x_poly = c(
            z$x_base[i],
            z$x_base[i],
            z$x[i]
          ),
          y_poly = c(
            z$y[i] - z$half_height[i],
            z$y[i] + z$half_height[i],
            z$y[i]
          ),
          x_tip = z$x[i],
          y_tip = z$y[i]
        )
      }
    )
  )
}

make_endpoint_phylogram_v28 <- function(
    family,
    tr,
    trait,
    top_host_orders,
    order_palette,
    show_legend = TRUE,
    folded = FALSE,
    breadth_legend_breaks = 1:5,
    label_top_fold_n = 4
) {
  trait <- prepare_endpoint_trait_display_v26(
    trait,
    top_host_orders
  )

  pal <- c(
    order_palette,
    Other = "#BDBDBD"
  )

  p0 <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular"
    )
  )

  td <- as.data.table(p0$data)

  endpoint <- trait[, .(
    label = .tip_id,
    branch_order = host_order_display,
    host_order_display,
    host_family_breadth,
    breadth_size_value,
    collapsed,
    collapsed_n
  )]

  m <- match(
    as.character(td$label),
    endpoint$label
  )

  td[, branch_order := NA_character_]
  td[, host_order_display := NA_character_]
  td[, host_family_breadth := NA_real_]
  td[, breadth_size_value := NA_real_]
  td[, collapsed := FALSE]
  td[, collapsed_n := 1L]

  hit <- which(
    !is.na(m) &
      td$isTip %in% TRUE
  )

  if (length(hit)) {
    td$branch_order[hit] <- endpoint$branch_order[m[hit]]
    td$host_order_display[hit] <- endpoint$host_order_display[m[hit]]
    td$host_family_breadth[hit] <- endpoint$host_family_breadth[m[hit]]
    td$breadth_size_value[hit] <- endpoint$breadth_size_value[m[hit]]
    td$collapsed[hit] <- endpoint$collapsed[m[hit]]
    td$collapsed_n[hit] <- endpoint$collapsed_n[m[hit]]
  }

  p0$data <- as.data.frame(td)
  p0$layers <- list()

  p <- p0 +
    geom_tree(
      aes(colour = branch_order),
      linewidth = 0.32
    ) +
    scale_colour_manual(
      values = pal,
      breaks = c(
        top_host_orders,
        "Other"
      ),
      na.value = "#444444",
      name = "Dominant host order",
      guide = if (show_legend) "legend" else "none"
    ) +
    labs(
      title = paste0(
        family,
        " (n = ",
        if (
          folded &&
          "collapsed_n" %in% names(trait)
        ) {
          sum(
            trait$collapsed_n,
            na.rm = TRUE
          )
        } else {
          length(tr$tip.label)
        },
        ")"
      )
    ) +
    theme_void(
      base_family = FONT_FAMILY
    ) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 2)
      ),
      plot.margin = margin(
        2,
        24,
        2,
        2
      )
    )

  tip_display <- as.data.table(p$data)[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      x = as.numeric(x),
      y = as.numeric(y),
      host_order_display,
      host_family_breadth,
      breadth_size_value,
      collapsed,
      collapsed_n
    )
  ]

  # Keep 1–5 visibly separated, but shrink the entire point scale.
  size_breaks_plot <- log1p(
    breadth_legend_breaks
  )

  regular_tip <- tip_display[
    collapsed %in% FALSE |
      is.na(collapsed)
  ]

  if (nrow(regular_tip)) {
    p <- p +
      geom_point(
        data = regular_tip,
        aes(
          x = x,
          y = y,
          size = breadth_size_value,
          colour = host_order_display
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        stroke = 0.38,
        na.rm = TRUE
      )
  }

  p <- p +
    scale_size_continuous(
      range = c(
        0.65,
        3.75
      ),
      breaks = size_breaks_plot,
      labels = as.character(
        breadth_legend_breaks
      ),
      name = "Host-family breadth",
      guide = if (show_legend) "legend" else "none"
    )

  if (folded) {
    tri <- make_fold_triangle_df_v28(
      p$data,
      tip_display
    )

    if (nrow(tri)) {
      p <- p +
        geom_polygon(
          data = tri,
          aes(
            x = x_poly,
            y = y_poly,
            group = fold_id
          ),
          inherit.aes = FALSE,
          fill = "white",
          colour = "#555555",
          linewidth = 0.38
        )

      fold_point <- unique(
        tri[, .(
          fold_id,
          x_tip,
          y_tip,
          host_order_display,
          breadth_size_value,
          collapsed_n
        )]
      )

      p <- p +
        geom_point(
          data = fold_point,
          aes(
            x = x_tip,
            y = y_tip,
            size = breadth_size_value,
            colour = host_order_display
          ),
          inherit.aes = FALSE,
          shape = 21,
          fill = "white",
          stroke = 0.42,
          na.rm = TRUE
        )

      fold_label <- copy(
        fold_point
      )

      setorder(
        fold_label,
        -collapsed_n
      )

      fold_label <- head(
        fold_label,
        label_top_fold_n
      )

      if (nrow(fold_label)) {
        p <- p +
          geom_text(
            data = fold_label,
            aes(
              x = x_tip,
              y = y_tip,
              label = paste0(
                "n=",
                collapsed_n
              )
            ),
            inherit.aes = FALSE,
            hjust = -0.22,
            vjust = 0.5,
            family = FONT_FAMILY,
            size = 1.70,
            colour = "#555555"
          )
      }
    }
  }

  p +
    coord_cartesian(
      clip = "off"
    ) +
    guides(
      colour = guide_legend(
        order = 1,
        title.position = "top",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.4,
          shape = NA
        )
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        nrow = 1,
        override.aes = list(
          colour = "black",
          fill = "white",
          shape = 21,
          stroke = 0.40
        )
      )
    ) +
    theme(
      # Top-left, wrapped rather than clipped.
      legend.position = if (
        show_legend
      ) "top" else "none",
      legend.justification = "left",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(
        size = TEXT_TINY,
        face = "bold",
        hjust = 0
      ),
      legend.text = element_text(
        size = TEXT_TINY
      ),
      legend.spacing.x = grid::unit(
        1.2,
        "mm"
      ),
      legend.spacing.y = grid::unit(
        0.6,
        "mm"
      ),
      legend.key.height = grid::unit(
        3.0,
        "mm"
      ),
      legend.key.width = grid::unit(
        3.6,
        "mm"
      )
    )
}


# ============================================================
# v2.9 helpers
# ============================================================

make_fold_triangle_df_v29 <- function(
    tree_data,
    trait_data
) {
  if (!nrow(trait_data)) return(data.table())

  td <- as.data.table(tree_data)

  tips <- td[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      node,
      parent,
      x,
      y
    )
  ]

  parent_x <- td[, .(
    parent = node,
    x_parent = x
  )]

  tips <- merge(
    tips,
    parent_x,
    by = "parent",
    all.x = TRUE
  )

  z <- merge(
    tips,
    trait_data[
      collapsed %in% TRUE,
      .(
        .tip_id,
        collapsed_n,
        host_order_display,
        host_family_breadth,
        breadth_size_value
      )
    ],
    by = ".tip_id"
  )

  if (!nrow(z)) return(data.table())

  z[
    !is.finite(x_parent),
    x_parent := x
  ]

  z[, edge_len := pmax(x - x_parent, 0)]

  # Smaller, shorter triangles than v2.8
  z[, x_base := x_parent + 0.72 * edge_len]

  z[, half_height := pmin(
    1.55,
    pmax(
      0.55,
      0.22 + 0.38 * log10(
        pmax(collapsed_n, 2)
      )
    )
  )]

  rbindlist(
    lapply(
      seq_len(nrow(z)),
      function(i) {
        data.table(
          fold_id = i,
          .tip_id = z$.tip_id[i],
          collapsed_n = z$collapsed_n[i],
          host_order_display = z$host_order_display[i],
          host_family_breadth = z$host_family_breadth[i],
          breadth_size_value = z$breadth_size_value[i],
          x_poly = c(
            z$x_base[i],
            z$x_base[i],
            z$x[i]
          ),
          y_poly = c(
            z$y[i] - z$half_height[i],
            z$y[i] + z$half_height[i],
            z$y[i]
          ),
          x_tip = z$x[i],
          y_tip = z$y[i]
        )
      }
    )
  )
}

make_tree_scale_bar_v29 <- function(
    tree_plot_data,
    y_offset = 6
) {
  td <- as.data.table(tree_plot_data)
  rng_x <- range(td$x, na.rm = TRUE)
  rng_y <- range(td$y, na.rm = TRUE)
  tree_span <- diff(rng_x)

  width <- signif(tree_span * 0.08, 2)
  if (!is.finite(width) || width <= 0) width <- 0.05

  list(
    x = rng_x[1] + tree_span * 0.02,
    y = rng_y[1] - y_offset,
    width = width
  )
}

make_endpoint_phylogram_v29 <- function(
    family,
    tr,
    trait,
    top_host_orders,
    order_palette,
    show_legend = TRUE,
    folded = FALSE,
    breadth_legend_breaks = 1:5,
    label_top_fold_n = 0,
    inside_legend = FALSE,
    add_scale_bar = TRUE,
    scale_y_offset = 6
) {
  trait <- prepare_endpoint_trait_display_v26(
    trait,
    top_host_orders
  )

  pal <- c(
    order_palette,
    Other = "#BDBDBD"
  )

  p0 <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular"
    )
  )

  td <- as.data.table(p0$data)

  endpoint <- trait[, .(
    label = .tip_id,
    branch_order = host_order_display,
    host_order_display,
    host_family_breadth,
    breadth_size_value,
    collapsed,
    collapsed_n
  )]

  m <- match(
    as.character(td$label),
    endpoint$label
  )

  td[, branch_order := NA_character_]
  td[, host_order_display := NA_character_]
  td[, host_family_breadth := NA_real_]
  td[, breadth_size_value := NA_real_]
  td[, collapsed := FALSE]
  td[, collapsed_n := 1L]

  hit <- which(
    !is.na(m) &
      td$isTip %in% TRUE
  )

  if (length(hit)) {
    td$branch_order[hit] <- endpoint$branch_order[m[hit]]
    td$host_order_display[hit] <- endpoint$host_order_display[m[hit]]
    td$host_family_breadth[hit] <- endpoint$host_family_breadth[m[hit]]
    td$breadth_size_value[hit] <- endpoint$breadth_size_value[m[hit]]
    td$collapsed[hit] <- endpoint$collapsed[m[hit]]
    td$collapsed_n[hit] <- endpoint$collapsed_n[m[hit]]
  }

  p0$data <- as.data.frame(td)
  p0$layers <- list()

  p <- p0 +
    geom_tree(
      aes(colour = branch_order),
      linewidth = 0.30
    ) +
    scale_colour_manual(
      values = pal,
      breaks = c(
        top_host_orders,
        "Other"
      ),
      na.value = "#444444",
      name = "Dominant host order",
      guide = if (show_legend) "legend" else "none"
    ) +
    labs(
      title = paste0(
        family,
        " (n = ",
        if (
          folded &&
          "collapsed_n" %in% names(trait)
        ) {
          sum(
            trait$collapsed_n,
            na.rm = TRUE
          )
        } else {
          length(tr$tip.label)
        },
        ")"
      )
    ) +
    theme_void(
      base_family = FONT_FAMILY
    ) +
    theme(
      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "italic",
        hjust = 0,
        margin = margin(b = 2)
      ),
      plot.margin = margin(
        5,
        18,
        16,
        2
      )
    )

  tip_display <- as.data.table(p$data)[
    isTip %in% TRUE,
    .(
      .tip_id = as.character(label),
      x = as.numeric(x),
      y = as.numeric(y),
      host_order_display,
      host_family_breadth,
      breadth_size_value,
      collapsed,
      collapsed_n
    )
  ]

  regular_tip <- tip_display[
    collapsed %in% FALSE |
      is.na(collapsed)
  ]

  size_breaks_plot <- log1p(
    breadth_legend_breaks
  )

  if (nrow(regular_tip)) {
    p <- p +
      geom_point(
        data = regular_tip,
        aes(
          x = x,
          y = y,
          size = breadth_size_value,
          colour = host_order_display
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        stroke = 0.34,
        na.rm = TRUE
      )
  }

  p <- p +
    scale_size_continuous(
      range = c(
        0.55,
        3.10
      ),
      breaks = size_breaks_plot,
      labels = as.character(
        breadth_legend_breaks
      ),
      name = "Host-family breadth",
      guide = if (show_legend) "legend" else "none"
    )

  if (folded) {
    tri <- make_fold_triangle_df_v29(
      p$data,
      tip_display
    )

    if (nrow(tri)) {
      p <- p +
        geom_polygon(
          data = tri,
          aes(
            x = x_poly,
            y = y_poly,
            group = fold_id
          ),
          inherit.aes = FALSE,
          fill = "white",
          colour = "#555555",
          linewidth = 0.34
        )

      fold_point <- unique(
        tri[, .(
          fold_id,
          x_tip,
          y_tip,
          host_order_display,
          breadth_size_value,
          collapsed_n
        )]
      )

      p <- p +
        geom_point(
          data = fold_point,
          aes(
            x = x_tip,
            y = y_tip,
            size = breadth_size_value,
            colour = host_order_display
          ),
          inherit.aes = FALSE,
          shape = 21,
          fill = "white",
          stroke = 0.36,
          na.rm = TRUE
        )

      if (
        isTRUE(label_top_fold_n > 0)
      ) {
        fold_label <- copy(
          fold_point
        )

        setorder(
          fold_label,
          -collapsed_n
        )

        fold_label <- head(
          fold_label,
          label_top_fold_n
        )

        if (nrow(fold_label)) {
          p <- p +
            geom_text(
              data = fold_label,
              aes(
                x = x_tip,
                y = y_tip,
                label = paste0(
                  "n=",
                  collapsed_n
                )
              ),
              inherit.aes = FALSE,
              hjust = -0.22,
              vjust = 0.5,
              family = FONT_FAMILY,
              size = 1.60,
              colour = "#555555"
            )
        }
      }
    }
  }

  if (add_scale_bar) {
    sb <- make_tree_scale_bar_v29(
      p$data,
      y_offset = scale_y_offset
    )

    p <- p +
      ggtree::geom_treescale(
        x = sb$x,
        y = sb$y,
        width = sb$width,
        fontsize = 2.0,
        linesize = 0.35,
        offset = 0
      )
  }

  p <- p +
    coord_cartesian(
      clip = "off"
    ) +
    guides(
      colour = guide_legend(
        order = 1,
        title.position = "top",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.3,
          shape = NA
        )
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        nrow = 1,
        override.aes = list(
          colour = "black",
          fill = "white",
          shape = 21,
          stroke = 0.34
        )
      )
    )

  if (inside_legend) {
    p <- p +
      theme(
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.background = element_rect(
          fill = scales::alpha("white", 0.0),
          colour = NA
        ),
        legend.title = element_text(
          size = TEXT_TINY,
          face = "bold",
          hjust = 0
        ),
        legend.text = element_text(
          size = TEXT_TINY
        ),
        legend.spacing.x = grid::unit(
          0.8,
          "mm"
        ),
        legend.spacing.y = grid::unit(
          0.5,
          "mm"
        ),
        legend.key.height = grid::unit(
          2.7,
          "mm"
        ),
        legend.key.width = grid::unit(
          3.0,
          "mm"
        )
      )
  } else {
    p <- p +
      theme(
        legend.position = if (
          show_legend
        ) "top" else "none",
        legend.justification = "left",
        legend.box.just = "left",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.title = element_text(
          size = TEXT_TINY,
          face = "bold",
          hjust = 0
        ),
        legend.text = element_text(
          size = TEXT_TINY
        )
      )
  }

  p
}


# ============================================================
# v2.10 final-publication helpers
# ============================================================

family_top_orders_v210 <- function(
    trait,
    top_n = 8
) {
  if (!nrow(trait)) return(character())

  x <- trait[
    !is.na(host_order) &
      nzchar(host_order),
    .N,
    by = host_order
  ]

  if (!nrow(x)) return(character())

  setorder(
    x,
    -N,
    host_order
  )

  head(
    x$host_order,
    top_n
  )
}

shared_order_palette_v210 <- function(
    ordered_levels
) {
  ordered_levels <- unique(
    ordered_levels[
      !is.na(ordered_levels) &
        nzchar(ordered_levels)
    ]
  )

  if (!length(ordered_levels)) {
    return(character())
  }

  # Keep the manuscript's first six colours stable, then extend with
  # qualitative HCL colours if Supplementary requires more host orders.
  fixed <- c(
    "#4E79A7",
    "#F28E2B",
    "#59A14F",
    "#E15759",
    "#B07AA1",
    "#76B7B2"
  )

  n <- length(ordered_levels)
  cols <- character(n)

  n_fixed <- min(
    length(fixed),
    n
  )

  cols[
    seq_len(n_fixed)
  ] <- fixed[
    seq_len(n_fixed)
  ]

  if (n > n_fixed) {
    extra_n <- n - n_fixed

    extra <- grDevices::hcl.colors(
      max(
        extra_n,
        3
      ),
      palette = "Dark 3"
    )[
      seq_len(extra_n)
    ]

    # Avoid accidental near-duplication of the first fixed colours by
    # rotating the generated extension deterministically.
    if (length(extra) > 1) {
      extra <- c(
        extra[-1],
        extra[1]
      )
    }

    cols[
      (n_fixed + 1):n
    ] <- extra
  }

  stats::setNames(
    cols,
    ordered_levels
  )
}

breadth_breaks_v210 <- function(
    trait
) {
  if (
    !nrow(trait) ||
      !"host_family_breadth" %in% names(trait)
  ) {
    return(c(1, 2, 4))
  }

  vals <- trait$host_family_breadth[
    is.finite(
      trait$host_family_breadth
    ) &
      trait$host_family_breadth > 0
  ]

  if (!length(vals)) {
    return(c(1, 2, 4))
  }

  vmax <- max(
    vals,
    na.rm = TRUE
  )

  if (vmax <= 5) {
    return(
      seq_len(
        max(
          1L,
          as.integer(
            ceiling(vmax)
          )
        )
      )
    )
  }

  core <- c(
    1,
    2,
    4,
    8
  )

  core <- core[
    core < vmax
  ]

  unique(
    c(
      core,
      vmax
    )
  )
}

nice_125_value_v210 <- function(
    x
) {
  if (
    !is.finite(x) ||
      x <= 0
  ) {
    return(0.1)
  }

  p <- floor(
    log10(x)
  )

  candidates <- as.vector(
    outer(
      c(
        10^(p - 1),
        10^p,
        10^(p + 1)
      ),
      c(
        1,
        2,
        5
      )
    )
  )

  candidates <- sort(
    unique(
      candidates[
        is.finite(candidates) &
          candidates > 0
      ]
    )
  )

  candidates[
    which.min(
      abs(
        candidates - x
      )
    )
  ]
}

tree_scale_bar_v210 <- function(
    tree_plot_data,
    target_fraction = 0.08,
    y_offset = 6
) {
  td <- as.data.table(
    tree_plot_data
  )

  xr <- range(
    td$x,
    na.rm = TRUE
  )

  yr <- range(
    td$y,
    na.rm = TRUE
  )

  span <- diff(xr)

  raw_width <- span *
    target_fraction

  width <- nice_125_value_v210(
    raw_width
  )

  list(
    x = xr[1] +
      span * 0.02,
    y = yr[1] -
      y_offset,
    width = width
  )
}

tag_plot_inset_v210 <- function(
    plot,
    label,
    left = 0.002,
    bottom = 0.935,
    right = 0.052,
    top = 0.997
) {
  patchwork::wrap_plots(
    plot
  ) +
    patchwork::inset_element(
      tag_strip(
        label
      ),
      left = left,
      bottom = bottom,
      right = right,
      top = top,
      align_to = "full",
      clip = FALSE
    )
}

make_endpoint_phylogram_v210 <- function(
    family,
    tr,
    trait,
    top_host_orders,
    order_palette,
    show_legend = TRUE,
    folded = FALSE,
    breadth_breaks = NULL,
    breadth_limit_max = NULL,
    inside_legend = FALSE,
    add_scale_bar = TRUE,
    scale_y_offset = 6,
    legend_position = "top"
) {
  trait <- prepare_endpoint_trait_display_v26(
    trait,
    top_host_orders
  )

  pal <- c(
    order_palette,
    Other = "#BDBDBD"
  )

  if (is.null(breadth_breaks)) {
    breadth_breaks <- breadth_breaks_v210(
      trait
    )
  }

  if (
    is.null(breadth_limit_max) ||
      !is.finite(breadth_limit_max)
  ) {
    raw_vals <- trait$host_family_breadth[
      is.finite(
        trait$host_family_breadth
      ) &
        trait$host_family_breadth > 0
    ]

    breadth_limit_max <- if (
      length(raw_vals)
    ) {
      max(
        raw_vals,
        na.rm = TRUE
      )
    } else {
      max(
        breadth_breaks
      )
    }
  }

  p0 <- suppressWarnings(
    ggtree::ggtree(
      tr,
      layout = "rectangular"
    )
  )

  td <- as.data.table(
    p0$data
  )

  endpoint <- trait[, .(
    label = .tip_id,
    branch_order = host_order_display,
    host_order_display,
    host_family_breadth,
    breadth_size_value,
    collapsed,
    collapsed_n
  )]

  m <- match(
    as.character(
      td$label
    ),
    endpoint$label
  )

  td[, branch_order := NA_character_]
  td[, host_order_display := NA_character_]
  td[, host_family_breadth := NA_real_]
  td[, breadth_size_value := NA_real_]
  td[, collapsed := FALSE]
  td[, collapsed_n := 1L]

  hit <- which(
    !is.na(m) &
      td$isTip %in% TRUE
  )

  if (length(hit)) {
    td$branch_order[hit] <- endpoint$branch_order[
      m[hit]
    ]
    td$host_order_display[hit] <- endpoint$host_order_display[
      m[hit]
    ]
    td$host_family_breadth[hit] <- endpoint$host_family_breadth[
      m[hit]
    ]
    td$breadth_size_value[hit] <- endpoint$breadth_size_value[
      m[hit]
    ]
    td$collapsed[hit] <- endpoint$collapsed[
      m[hit]
    ]
    td$collapsed_n[hit] <- endpoint$collapsed_n[
      m[hit]
    ]
  }

  p0$data <- as.data.frame(
    td
  )
  p0$layers <- list()

  p <- p0 +
    geom_tree(
      aes(
        colour = branch_order
      ),
      linewidth = 0.30
    ) +
    scale_colour_manual(
      values = pal,
      breaks = c(
        top_host_orders,
        "Other"
      ),
      na.value = "#444444",
      name = "Dominant host order",
      guide = if (
        show_legend
      ) "legend" else "none"
    ) +
    labs(
      title = paste0(
        family,
        " (n = ",
        if (
          folded &&
          "collapsed_n" %in%
            names(trait)
        ) {
          sum(
            trait$collapsed_n,
            na.rm = TRUE
          )
        } else {
          length(
            tr$tip.label
          )
        },
        ")"
      )
    ) +
    theme_void(
      base_family =
        FONT_FAMILY
    ) +
    theme(
      plot.title =
        element_text(
          family =
            FONT_FAMILY,
          size =
            TEXT_SMALL,
          face =
            "italic",
          hjust =
            0,
          margin =
            margin(
              b = 2
            )
        ),
      plot.margin =
        margin(
          5,
          18,
          16,
          2
        )
    )

  tip_display <- as.data.table(
    p$data
  )[
    isTip %in% TRUE,
    .(
      .tip_id =
        as.character(
          label
        ),
      x =
        as.numeric(x),
      y =
        as.numeric(y),
      host_order_display,
      host_family_breadth,
      breadth_size_value,
      collapsed,
      collapsed_n
    )
  ]

  regular_tip <- tip_display[
    collapsed %in% FALSE |
      is.na(collapsed)
  ]

  size_breaks_plot <- log1p(
    breadth_breaks
  )

  if (nrow(regular_tip)) {
    p <- p +
      geom_point(
        data =
          regular_tip,
        aes(
          x = x,
          y = y,
          size =
            breadth_size_value,
          colour =
            host_order_display
        ),
        inherit.aes =
          FALSE,
        shape = 21,
        fill = "white",
        stroke = 0.34,
        na.rm = TRUE
      )
  }

  p <- p +
    scale_size_continuous(
      range = c(
        0.55,
        3.10
      ),
      limits = log1p(
        c(
          1,
          breadth_limit_max
        )
      ),
      breaks =
        size_breaks_plot,
      labels =
        as.character(
          breadth_breaks
        ),
      name =
        "Host-family breadth",
      guide = if (
        show_legend
      ) "legend" else "none"
    )

  if (folded) {
    tri <- make_fold_triangle_df_v29(
      p$data,
      tip_display
    )

    if (nrow(tri)) {
      p <- p +
        geom_polygon(
          data = tri,
          aes(
            x = x_poly,
            y = y_poly,
            group = fold_id
          ),
          inherit.aes =
            FALSE,
          fill = "white",
          colour = "#555555",
          linewidth = 0.34
        )

      fold_point <- unique(
        tri[, .(
          fold_id,
          x_tip,
          y_tip,
          host_order_display,
          breadth_size_value,
          collapsed_n
        )]
      )

      p <- p +
        geom_point(
          data =
            fold_point,
          aes(
            x = x_tip,
            y = y_tip,
            size =
              breadth_size_value,
            colour =
              host_order_display
          ),
          inherit.aes =
            FALSE,
          shape = 21,
          fill = "white",
          stroke = 0.36,
          na.rm = TRUE
        )
    }
  }

  if (add_scale_bar) {
    sb <- tree_scale_bar_v210(
      p$data,
      target_fraction =
        0.08,
      y_offset =
        scale_y_offset
    )

    p <- p +
      ggtree::geom_treescale(
        x = sb$x,
        y = sb$y,
        width = sb$width,
        fontsize = 2.0,
        linesize = 0.35,
        offset = 0
      )
  }

  p <- p +
    coord_cartesian(
      clip = "off"
    ) +
    guides(
      colour =
        guide_legend(
          order = 1,
          title.position =
            "top",
          nrow = 2,
          byrow = TRUE,
          override.aes =
            list(
              linewidth =
                1.30,
              shape =
                NA
            )
        ),
      size =
        guide_legend(
          order = 2,
          title.position =
            "top",
          nrow = 1,
          override.aes =
            list(
              colour =
                "black",
              fill =
                "white",
              shape =
                21,
              stroke =
                0.34
            )
        )
    )

  common_legend_theme <- theme(
    # IMPORTANT: match A's PUB_THEME legend typography.
    legend.title =
      element_text(
        size =
          TEXT_SMALL,
        face =
          "bold",
        colour =
          "black",
        hjust = 0
      ),
    legend.text =
      element_text(
        size =
          TEXT_SMALL,
        face =
          "plain",
        colour =
          "black"
      ),
    legend.spacing.x =
      grid::unit(
        0.9,
        "mm"
      ),
    legend.spacing.y =
      grid::unit(
        0.5,
        "mm"
      ),
    legend.key.height =
      grid::unit(
        3.0,
        "mm"
      ),
    legend.key.width =
      grid::unit(
        3.5,
        "mm"
      )
  )

  if (inside_legend) {
    # Start farther to the right than the C panel tag.
    p <- p +
      theme(
        legend.position =
          c(
            0.075,
            0.975
          ),
        legend.justification =
          c(
            0,
            1
          ),
        legend.direction =
          "horizontal",
        legend.box =
          "vertical",
        legend.background =
          element_rect(
            fill =
              scales::alpha(
                "white",
                0.0
              ),
            colour =
              NA
          )
      ) +
      common_legend_theme
  } else {
    p <- p +
      theme(
        legend.position =
          if (
            show_legend
          ) legend_position
          else "none",
        legend.justification =
          "left",
        legend.box.just =
          "left",
        legend.direction =
          "horizontal",
        legend.box =
          "vertical"
      ) +
      common_legend_theme
  }

  p
}

terminal_branch_qc_v210 <- function(
    tr,
    family,
    trait = data.table(),
    top_fraction = 0.01
) {
  e <- as.data.table(
    tr$edge
  )

  setnames(
    e,
    c(
      "parent",
      "child"
    )
  )

  e[, edge_length :=
      tr$edge.length]

  ntip <- length(
    tr$tip.label
  )

  tip_edges <- e[
    child <= ntip
  ]

  tip_edges[, .tip_id :=
      tr$tip.label[
        child
      ]]

  if (!nrow(tip_edges)) {
    return(
      data.table()
    )
  }

  cutoff <- stats::quantile(
    tip_edges$edge_length,
    probs =
      1 - top_fraction,
    na.rm = TRUE,
    names = FALSE
  )

  out <- tip_edges[
    edge_length >= cutoff
  ]

  out[, virus_family :=
      family]

  out[, top_fraction :=
      top_fraction]

  out[, branch_length_cutoff :=
      cutoff]

  if (
    nrow(trait)
  ) {
    keep_cols <- intersect(
      c(
        ".tip_id",
        "virus_genus",
        "host_order",
        "host_family_breadth"
      ),
      names(trait)
    )

    out <- merge(
      out,
      unique(
        trait[
          ,
          ..keep_cols
        ]
      ),
      by = ".tip_id",
      all.x = TRUE
    )
  }

  setorder(
    out,
    -edge_length
  )

  out
}


# ============================================================
# v2.11 final-layout helpers
# ============================================================

balanced_breadth_breaks_v211 <- function(
    vmax
) {
  if (
    !is.finite(vmax) ||
      vmax <= 0
  ) {
    return(c(1))
  }

  vmax <- as.numeric(vmax)

  if (vmax <= 5) {
    return(
      seq_len(
        max(
          1L,
          as.integer(
            floor(vmax)
          )
        )
      )
    )
  }

  # Publication-facing regular landmarks.
  # We deliberately do NOT append irregular maxima such as 13, 14 or 17.
  regular <- c(
    1,
    5,
    10,
    15,
    20,
    25,
    30
  )

  regular[
    regular <= vmax
  ]
}

scale_y_offset_v211 <- function(
    tip_n,
    minimum = 9
) {
  max(
    minimum,
    ceiling(
      tip_n * 0.020
    )
  )
}

tag_plot_inset_v211 <- function(
    plot,
    label,
    left = 0.006,
    bottom = 0.925,
    right = 0.065,
    top = 0.995
) {
  patchwork::wrap_plots(
    plot
  ) +
    patchwork::inset_element(
      tag_strip(
        label
      ),
      left = left,
      bottom = bottom,
      right = right,
      top = top,
      align_to = "full",
      clip = FALSE
    )
}

terminal_branch_qc_v211 <- function(
    tr,
    family,
    trait = data.table(),
    top_fraction = 0.01
) {
  e <- as.data.table(
    tr$edge
  )

  setnames(
    e,
    c(
      "parent",
      "child"
    )
  )

  e[, edge_length :=
      tr$edge.length]

  ntip <- length(
    tr$tip.label
  )

  tip_edges <- e[
    child <= ntip
  ]

  tip_edges[, .tip_id :=
      tr$tip.label[
        child
      ]]

  if (!nrow(tip_edges)) {
    return(
      data.table()
    )
  }

  med_len <- stats::median(
    tip_edges$edge_length,
    na.rm = TRUE
  )

  p95_len <- stats::quantile(
    tip_edges$edge_length,
    probs = 0.95,
    na.rm = TRUE,
    names = FALSE
  )

  cutoff <- stats::quantile(
    tip_edges$edge_length,
    probs = 1 - top_fraction,
    na.rm = TRUE,
    names = FALSE
  )

  out <- tip_edges[
    edge_length >= cutoff
  ]

  out[, virus_family :=
      family]

  out[, top_fraction :=
      top_fraction]

  out[, terminal_branch_median :=
      med_len]

  out[, terminal_branch_p95 :=
      p95_len]

  out[, branch_length_cutoff :=
      cutoff]

  if (
    is.finite(med_len) &&
      med_len > 0
  ) {
    out[, branch_to_median_ratio :=
        edge_length / med_len]
  } else {
    out[, branch_to_median_ratio :=
        NA_real_]
  }

  if (
    nrow(trait)
  ) {
    keep_cols <- intersect(
      c(
        ".tip_id",
        "virus_genus",
        "host_order",
        "host_family_breadth"
      ),
      names(trait)
    )

    out <- merge(
      out,
      unique(
        trait[
          ,
          ..keep_cols
        ]
      ),
      by = ".tip_id",
      all.x = TRUE
    )
  }

  setorder(
    out,
    -edge_length
  )

  out
}
