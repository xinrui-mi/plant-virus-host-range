suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

snake_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

standardize_names <- function(x) {
  setnames(x, names(x), make.unique(snake_name(names(x)), sep = "_"))
  x
}

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
    stop("Could not identify required column. Available columns: ", paste(nms, collapse = ", "))
  }
  NULL
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  x <- fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  standardize_names(x)
}

locate_table <- function(exact_name, fallback_regex = NULL, required = TRUE) {
  p <- file.path(TABLE_DIR, exact_name)
  if (file.exists(p)) return(p)
  if (!is.null(fallback_regex)) {
    hits <- list.files(TABLE_DIR, pattern = fallback_regex, full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) return(sort(hits)[1])
  }
  if (required) stop("Cannot locate required frozen table: ", exact_name)
  NA_character_
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  z <- tolower(trimws(as.character(x)))
  z %in% c("1", "true", "t", "yes", "y", "candidate", "selected", "pass", "retained")
}

mode_or_na <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

first_finite <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y <- y[is.finite(y)]
  if (!length(y)) return(NA_real_)
  y[1]
}

save_triplet <- function(plot, stem, outdir, width_mm, height_mm, dpi = 600) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(outdir, paste0(stem, ".pdf")), plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(outdir, paste0(stem, ".png")), plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(outdir, paste0(stem, ".tiff")), plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, compression = "lzw", bg = "white", limitsize = FALSE
  )
  invisible(TRUE)
}

extract_official_thresholds <- function(thr) {
  nms <- names(thr)
  get_wide <- function(exact, regex) {
    c <- find_col(nms, exact = exact, regex = regex, required = FALSE)
    if (is.null(c)) NA_real_ else first_finite(thr[[c]])
  }

  breadth <- get_wide(
    c("official_host_breadth_cutoff", "host_breadth_cutoff", "host_family_breadth_cutoff"),
    c("official.*host.*breadth.*cut", "host.*family.*breadth.*cut", "host.*breadth.*cut")
  )
  bridge <- get_wide(
    c("official_bridge_cutoff", "bridge_cutoff", "network_bridge_cutoff"),
    c("official.*bridge.*cut", "bridge.*cut")
  )
  support <- get_wide(
    c("official_min_vc80_support", "support_cutoff", "vc80_support_cutoff", "cluster80_support_cutoff"),
    c("official.*support", "min.*vc80.*support", "support.*cut")
  )

  if (!all(is.finite(c(breadth, bridge, support)))) {
    stop(
      "Could not recover all three official Fig.5 cutoffs from 05_candidate_thresholds.tsv. ",
      "Detected breadth=", breadth, ", bridge=", bridge, ", support=", support
    )
  }

  data.table(
    threshold = c("Host-family breadth", "Network bridge score", "VC80 support"),
    value = c(breadth, bridge, support)
  )
}

standardize_metrics <- function(x, thresholds) {
  nms <- names(x)
  genus_col <- find_col(nms,
    exact = c("virus_node", "virus_genus", "current_ictv_genus", "ictv_genus", "genus"),
    regex = c("virus.*node", "virus.*genus", "ictv.*genus", "^genus$"))
  family_col <- find_col(nms,
    exact = c("virus_family", "current_ictv_family", "ictv_family", "family"),
    regex = c("virus.*family", "ictv.*family"), required = FALSE)
  breadth_col <- find_col(nms,
    exact = c("host_family_n", "host_family_breadth", "host_breadth"),
    regex = c("host.*family.*_n$", "host.*family.*breadth", "host.*breadth"))
  bridge_col <- find_col(nms,
    exact = c("bridge_score", "network_bridge_score"), regex = c("bridge.*score"))
  support_col <- find_col(nms,
    exact = c("qc_pass_vc80_support", "vc80_support", "cluster80_support"),
    regex = c("qc.*vc80.*support", "vc80.*support", "cluster80.*support"))
  cand_col <- find_col(nms,
    exact = c("candidate", "is_candidate", "priority_candidate", "main_candidate"),
    regex = c("(^|_)candidate$", "is.*candidate", "priority.*candidate"), required = FALSE)
  module_col <- find_col(nms,
    exact = c("module_id", "module", "network_module"), regex = c("module"), required = FALSE)

  out <- data.table(
    virus_genus = as.character(x[[genus_col]]),
    virus_family = if (is.null(family_col)) NA_character_ else as.character(x[[family_col]]),
    host_family_breadth = suppressWarnings(as.numeric(x[[breadth_col]])),
    bridge_score = suppressWarnings(as.numeric(x[[bridge_col]])),
    vc80_support = suppressWarnings(as.numeric(x[[support_col]])),
    module_id = if (is.null(module_col)) NA_character_ else as.character(x[[module_col]])
  )

  b  <- thresholds[threshold == "Host-family breadth", value]
  br <- thresholds[threshold == "Network bridge score", value]
  s  <- thresholds[threshold == "VC80 support", value]
  derived <- with(out, host_family_breadth >= b & bridge_score >= br & vc80_support >= s)

  if (!is.null(cand_col)) {
    official <- truthy(x[[cand_col]])
    mismatch <- sum(official != derived, na.rm = TRUE)
    if (mismatch > 0) {
      warning(mismatch, " official candidate flags differ from cutoff-derived flags; official flags retained.")
    }
    out[, candidate := official]
  } else {
    out[, candidate := derived]
  }

  out <- out[!is.na(virus_genus) & nzchar(virus_genus)]
  unique(out, by = "virus_genus")
}

standardize_edge_table <- function(edge) {
  nms <- names(edge)
  genus_col <- find_col(nms,
    exact = c("virus_node", "virus_genus", "genus", "from"),
    regex = c("virus.*node", "virus.*genus", "^genus$", "^from$"))
  host_col <- find_col(nms,
    exact = c("host_family", "to"), regex = c("host.*family", "^to$"))
  weight_col <- find_col(nms,
    exact = c("weight", "edge_weight", "vc80_n", "vc80_support"),
    regex = c("^weight$", "edge.*weight", "vc80.*_n$", "vc80.*support"), required = FALSE)

  out <- data.table(
    virus_genus = as.character(edge[[genus_col]]),
    host_family = as.character(edge[[host_col]]),
    edge_weight = if (is.null(weight_col)) 1 else suppressWarnings(as.numeric(edge[[weight_col]]))
  )
  out[!is.finite(edge_weight), edge_weight := 1]
  out[!is.na(virus_genus) & nzchar(virus_genus) & !is.na(host_family) & nzchar(host_family)]
}

standardize_modules <- function(x) {
  nms <- names(x)
  module_col <- find_col(nms, exact = c("module_id", "module", "network_module"), regex = c("module"))
  node_col <- find_col(nms,
    exact = c("name", "node", "node_name"), regex = c("^name$", "^node$", "node.*name"))

  raw_node <- as.character(x[[node_col]])
  node_type <- rep("unknown", length(raw_node))
  node_type[grepl("^[Vv]:", raw_node)] <- "virus_genus"
  node_type[grepl("^[Hh]:", raw_node)] <- "host_family"
  clean_node <- sub("^[VvHh]:", "", raw_node)

  data.table(
    node = clean_node,
    node_raw = raw_node,
    node_type = node_type,
    module_id = as.character(x[[module_col]])
  )[!is.na(node) & nzchar(node) & !is.na(module_id)]
}

make_taxonomy_maps <- function(master) {
  nms <- names(master)
  genus_col <- find_col(nms,
    exact = c("virus_genus", "current_ictv_genus", "ictv_genus"),
    regex = c("virus.*genus", "ictv.*genus"))
  family_col <- find_col(nms,
    exact = c("virus_family", "current_ictv_family", "ictv_family"),
    regex = c("virus.*family", "ictv.*family"))
  hfamily_col <- find_col(nms, exact = c("host_family"), regex = c("host.*family"))
  horder_col <- find_col(nms, exact = c("host_order"), regex = c("host.*order"))

  genus_family <- master[
    !is.na(get(genus_col)) & nzchar(trimws(as.character(get(genus_col)))) &
    !is.na(get(family_col)) & nzchar(trimws(as.character(get(family_col)))),
    .(virus_family = mode_or_na(get(family_col))),
    by = .(virus_genus = as.character(get(genus_col)))
  ]

  host_order <- master[
    !is.na(get(hfamily_col)) & nzchar(trimws(as.character(get(hfamily_col)))) &
    !is.na(get(horder_col)) & nzchar(trimws(as.character(get(horder_col)))),
    .(host_order = mode_or_na(get(horder_col))),
    by = .(host_family = as.character(get(hfamily_col)))
  ]

  list(genus_family = genus_family, host_order = host_order)
}

candidate_family_palette <- function(metrics) {
  fam_order <- metrics[candidate == TRUE & !is.na(virus_family) & nzchar(virus_family),
                       .N, by = virus_family][order(-N, virus_family), virus_family]
  if (!length(fam_order)) return(c("Other" = "#BDBDBD"))

  # Fixed family-colour mapping used consistently in Panels C and D.
  # High-separation colours are chosen so biologically different families do not look identical.
  fixed_map <- c(
    "Secoviridae" = "#E76F92",
    "Bromoviridae" = "#D98C28",
    "Alphaflexiviridae" = "#A79C00",
    "Betaflexiviridae" = "#58B82D",
    "Caulimoviridae" = "#1FB989",
    "Geminiviridae" = "#18B8D6",
    "Potyviridae" = "#2097E6",
    "Solemoviridae" = "#9A6AD9",
    "Tospoviridae" = "#D968C3",
    "Other" = "#BDBDBD"
  )

  cols <- fixed_map[intersect(names(fixed_map), c(fam_order, "Other"))]
  missing_fams <- setdiff(fam_order, names(cols))
  if (length(missing_fams)) {
    extra_cols <- grDevices::hcl.colors(length(missing_fams), palette = "Dynamic")
    cols <- c(cols, setNames(extra_cols, missing_fams))
  }
  cols <- cols[c(fam_order, "Other")]
  cols
}


# Display-only transform: compress values below the official breadth cutoff so
# the candidate-rich upper-breadth region receives more horizontal space.
# This changes only display geometry; it does not alter any analytical value.
focus_above_cutoff_trans <- function(cutoff, lower_slope = 0.50) {
  if (!is.finite(cutoff) || cutoff <= 0) stop("focus_above_cutoff_trans requires a positive finite cutoff")
  if (!is.finite(lower_slope) || lower_slope <= 0 || lower_slope >= 1) stop("lower_slope must be between 0 and 1")
  forward <- function(x) {
    x <- as.numeric(x)
    ifelse(x <= cutoff, x * lower_slope, cutoff * lower_slope + (x - cutoff))
  }
  inverse <- function(y) {
    y <- as.numeric(y)
    hinge <- cutoff * lower_slope
    ifelse(y <= hinge, y / lower_slope, cutoff + (y - hinge))
  }
  scales::trans_new(
    name = paste0("focus_above_", format(cutoff, trim = TRUE)),
    transform = forward,
    inverse = inverse,
    domain = c(0, Inf)
  )
}

# Three readable, data-derived legend landmarks. No analytical cutoff is
# encoded here; these are presentation-only representative values.
auto_three_legend_breaks <- function(x) {
  x <- sort(unique(as.numeric(x[is.finite(x) & x > 0])))
  if (!length(x)) return(numeric())
  if (length(x) <= 3) return(x)
  q <- as.numeric(stats::quantile(x, probs = c(0.20, 0.55, 0.90), names = FALSE, type = 7))
  br <- unique(signif(q, 1))
  br <- br[br >= min(x) & br <= max(x)]
  if (length(br) < 3) {
    q2 <- unique(round(as.numeric(stats::quantile(x, probs = c(0.15, 0.50, 0.95), names = FALSE, type = 7))))
    br <- unique(c(br, q2[q2 >= min(x) & q2 <= max(x)]))
  }
  br <- sort(unique(br))
  if (length(br) > 3) br <- br[round(seq(1, length(br), length.out = 3))]
  br
}

# Data-derived axis breaks for the focused breadth display. The official
# cutoff itself is always included and comes from 05_candidate_thresholds.tsv.
focus_breadth_breaks <- function(x, cutoff) {
  x <- as.numeric(x[is.finite(x)])
  if (!length(x)) return(cutoff)
  lo <- min(x); hi <- max(x)
  nice1 <- function(v) {
    if (!is.finite(v) || v == 0) return(v)
    signif(v, 1)
  }
  raw <- c(
    lo,
    nice1(cutoff / 3),
    nice1(2 * cutoff / 3),
    cutoff,
    nice1(cutoff + (hi - cutoff) / 3),
    nice1(cutoff + 2 * (hi - cutoff) / 3),
    hi
  )
  raw <- sort(unique(raw[is.finite(raw) & raw >= lo & raw <= hi]))
  raw
}

# Helpers imported from the frozen candidate-contextualization plotting package.
# They only parse already-fitted model outputs; no model is re-fit here.
flag_from_any <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x > 0)
  z <- tolower(trimws(as.character(x)))
  false_values <- c("", "0", "false", "f", "no", "n", "none", "absent",
                    "no record", "no_record", "not found", "not_found",
                    "not assigned", "not_assigned", "unknown", "na", "n/a")
  !is.na(z) & !(z %in% false_values)
}

log_tick_labels <- function(x) {
  x <- as.numeric(x)
  vapply(x, function(z) {
    if (!is.finite(z)) return("")
    if (z < 1) formatC(z, format = "f", digits = 1)
    else format(round(z), scientific = FALSE, trim = TRUE, big.mark = ",")
  }, character(1))
}
