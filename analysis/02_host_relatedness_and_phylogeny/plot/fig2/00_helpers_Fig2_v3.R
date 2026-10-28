suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ape)
  library(patchwork)
})

FINAL12_FAMILIES <- c(
  "Alphaflexiviridae",
  "Betaflexiviridae",
  "Caulimoviridae",
  "Closteroviridae",
  "Geminiviridae",
  "Mitoviridae",
  "Partitiviridae",
  "Potyviridae",
  "Rhabdoviridae",
  "Secoviridae",
  "Solemoviridae",
  "Tombusviridae"
)

resolve_final12_tree_file <- function(fig3_root, family) {
  final_root <- file.path(
    fig3_root,
    "pipeline_results",
    "03K_family_specific_markers",
    "08_final_12family_trees"
  )

  # The release ships the twelve final treefiles under data/, which is where the
  # analysis modules read them. Keep the 03K location first so the cluster
  # provenance is unchanged whenever it is populated.
  final_repo_root <- dirname(dirname(normalizePath(fig3_root, winslash = "/", mustWork = TRUE)))
  released_root <- file.path(
    final_repo_root, "data", "analysis_data", "phylogeny", "twelve_family", "trees"
  )

  for (root in c(final_root, released_root)) {
    hits <- character()
    if (dir.exists(root)) {
      hits <- list.files(
        root,
        pattern = paste0("^", family, "\\.treefile$"),
        recursive = TRUE,
        full.names = TRUE
      )
      hits <- hits[file.exists(hits)]
    }
    if (length(hits)) {
      return(normalizePath(hits[1], winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Final 03K treefile not found for ", family, ". Looked for:\n  ",
    final_root, "/**/", family, ".treefile\n  ",
    released_root, "/", family, "/", family, ".treefile"
  )
}

final12_candidate_trait_files <- function(fig3_root, family) {
  roots <- c(
    file.path(fig3_root, "pipeline_results", "03K_family_specific_markers"),
    file.path(fig3_root, "pipeline_results", "03L_FINAL_12family_host_spectrum")
  )
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(character())

  ff <- unique(unlist(lapply(
    roots,
    function(rr) list.files(rr, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  ), use.names = FALSE))

  ff <- ff[file.exists(ff)]
  ff <- ff[grepl(family, ff, fixed = TRUE)]

  preferred <- ff[grepl("tree_tip_host_spectrum|annotation|host|membership", basename(ff), ignore.case = TRUE)]
  unique(c(preferred, ff))
}

candidate_id_columns <- function(nms) {
  low <- tolower(nms)
  exact <- c(
    ".tip_id", "tip_id", "tree_tip", "tree_tip_id",
    "representative_sequence_id", "representative_id", "representative",
    "sequence_id", "source_sequence_id", "marker_sequence_id",
    "marker_id", "orf_id", "candidate_id", "aa95_rep_id", "rep_id",
    "label", "rep_accession", "accession"
  )
  idx <- integer()
  for (x in exact) idx <- c(idx, which(low == tolower(x)))
  idx <- c(
    idx,
    grep("(^|[_.])(tip|representative|marker|orf|candidate|rep).*id$|sequence.*id$|rep.*accession$|^accession$", low, perl = TRUE)
  )
  unique(idx)
}

best_id_col_for_tree <- function(x, tree_labels) {
  idx <- candidate_id_columns(names(x))
  if (!length(idx)) return(list(col = NA_character_, overlap_n = 0L, overlap_fraction = 0))

  scores <- lapply(idx, function(i) {
    vals <- unique(as.character(x[[i]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    n <- length(intersect(vals, tree_labels))
    list(col = names(x)[i], overlap_n = n, overlap_fraction = n / max(1L, length(tree_labels)))
  })

  best <- scores[[which.max(vapply(scores, function(z) z$overlap_n, numeric(1)))]]
  best
}

find_breadth_col_v3 <- function(nms) {
  find_col(
    nms,
    exact = c(
      "host_family_n", "host_family_breadth", "host_family_coverage_n",
      "host_family_count", "host_family_count_n"
    ),
    regex = c(
      "^host.*family.*n$", "host.*family.*breadth", "host.*family.*coverage", "host.*family.*count"
    ),
    required = FALSE
  )
}

find_order_col_v3 <- function(nms) {
  # Only scalar order annotations are accepted. We deliberately do not use
  # 'host_orders' (plural set/list), because it is not a dominant order.
  find_col(
    nms,
    exact = c("dominant_host_order", "host_order", "top_host_order"),
    regex = c("^dominant.*host.*order$", "^host_order$", "^top.*host.*order$"),
    required = FALSE
  )
}


find_order_set_col_v31 <- function(nms) {
  find_col(
    nms,
    exact = c(
      "host_orders", "host_order_set", "host_orders_set",
      "host_order_list", "host_orders_list"
    ),
    regex = c("^host.*orders$", "host.*order.*set", "host.*order.*list"),
    required = FALSE
  )
}

collapse_order_value_v31 <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) return(NA_character_)
  z <- unlist(strsplit(as.character(x), "[;|]", perl = TRUE), use.names = FALSE)
  z <- trimws(z)
  z <- unique(z[!is.na(z) & nzchar(z) & !z %in% c("NA", "N/A", ".")])
  if (!length(z)) return(NA_character_)
  if (length(z) == 1L) return(z[1])
  "Multiple host orders"
}

read_best_order_component_v31 <- function(files, tree_labels, family) {
  candidates <- list()

  for (f in files) {
    h <- tryCatch(fread(f, sep = "\t", nrows = 0, showProgress = FALSE), error = function(e) NULL)
    if (is.null(h)) next

    scalar_col <- find_order_col_v3(names(h))
    set_col <- find_order_set_col_v31(names(h))
    if (is.null(scalar_col) && is.null(set_col)) next

    x <- tryCatch(fread(f, sep = "\t", showProgress = FALSE), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) next

    id <- best_id_col_for_tree(x, tree_labels)
    if (!is.finite(id$overlap_n) || id$overlap_n <= 0) next

    use_col <- if (!is.null(scalar_col)) scalar_col else set_col
    method <- if (!is.null(scalar_col)) "existing scalar host-order field" else "host-order set: singleton retained; multiple orders -> Other"

    vals <- vapply(x[[use_col]], collapse_order_value_v31, character(1))
    valid_n <- sum(!is.na(vals) & nzchar(vals))
    single_n <- sum(!is.na(vals) & nzchar(vals) & vals != "Multiple host orders")

    candidates[[length(candidates) + 1L]] <- list(
      path = f,
      table = x,
      id_col = id$col,
      trait_col = use_col,
      overlap_n = id$overlap_n,
      overlap_fraction = id$overlap_fraction,
      processed_order = vals,
      valid_n = valid_n,
      single_n = single_n,
      method = method,
      scalar_priority = if (!is.null(scalar_col)) 1L else 0L
    )
  }

  if (!length(candidates)) {
    stop(
      "Could not resolve a final per-tip host-order annotation for ", family, ".\n",
      "No scalar host-order field or host-order set with tree-tip overlap was found.\n",
      "Candidate files inspected:\n", paste(files, collapse = "\n")
    )
  }

  # Prefer actual scalar annotations when they overlap the final tree; otherwise
  # use the final host-order set without inventing a dominant category.
  score <- vapply(
    candidates,
    function(z) z$overlap_n * 100000 + z$scalar_priority * 10000 + z$single_n,
    numeric(1)
  )
  candidates[[which.max(score)]]
}

read_best_trait_component <- function(files, tree_labels, family, kind = c("breadth", "order")) {
  kind <- match.arg(kind)
  candidates <- list()

  for (f in files) {
    h <- tryCatch(fread(f, sep = "\t", nrows = 0, showProgress = FALSE), error = function(e) NULL)
    if (is.null(h)) next

    tcol <- if (kind == "breadth") find_breadth_col_v3(names(h)) else find_order_col_v3(names(h))
    if (is.null(tcol)) next

    x <- tryCatch(fread(f, sep = "\t", showProgress = FALSE), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) next

    id <- best_id_col_for_tree(x, tree_labels)
    if (!is.finite(id$overlap_n) || id$overlap_n <= 0) next

    vals <- x[[tcol]]
    if (kind == "order") {
      # Reject set/list-like values masquerading as a scalar host_order field.
      vv <- as.character(vals)
      vv <- vv[!is.na(vv) & nzchar(vv)]
      multi_fraction <- if (length(vv)) mean(grepl("[;|]", vv)) else 0
      if (multi_fraction > 0.10) next
    }

    candidates[[length(candidates) + 1L]] <- list(
      path = f,
      table = x,
      id_col = id$col,
      trait_col = tcol,
      overlap_n = id$overlap_n,
      overlap_fraction = id$overlap_fraction
    )
  }

  if (!length(candidates)) {
    stop(
      "Could not resolve a final per-tip ", kind, " annotation for ", family, ".\n",
      "The plotting script refuses to invent/fallback silently.\n",
      "Candidate files inspected:\n", paste(files, collapse = "\n")
    )
  }

  score <- vapply(candidates, function(z) z$overlap_n, numeric(1))
  candidates[[which.max(score)]]
}

prepare_final12_family_traits_v3 <- function(fig3_root, family, tree_file = NULL) {
  if (is.null(tree_file)) tree_file <- resolve_final12_tree_file(fig3_root, family)
  tr <- ape::read.tree(tree_file)
  tips <- as.character(tr$tip.label)
  files <- final12_candidate_trait_files(fig3_root, family)

  if (!length(files)) stop("No final 03K/03L TSV candidates found for ", family)

  br <- read_best_trait_component(files, tips, family, "breadth")
  od <- read_best_order_component_v31(files, tips, family)

  xb <- br$table[, .(
    .tip_id = as.character(get(br$id_col)),
    host_family_breadth = suppressWarnings(as.numeric(get(br$trait_col)))
  )]
  xb <- xb[!is.na(.tip_id) & nzchar(.tip_id)]
  xb <- xb[order(-is.finite(host_family_breadth))]
  xb <- unique(xb, by = ".tip_id")

  xo <- data.table(
    .tip_id = as.character(od$table[[od$id_col]]),
    host_order = as.character(od$processed_order)
  )
  xo[host_order %in% c("", "NA", "N/A", "."), host_order := NA_character_]
  xo <- xo[!is.na(.tip_id) & nzchar(.tip_id)]
  xo <- xo[order(-(!is.na(host_order) & nzchar(host_order)))]
  xo <- unique(xo, by = ".tip_id")

  out <- data.table(.tip_id = tips)
  out <- merge(out, xb, by = ".tip_id", all.x = TRUE)
  out <- merge(out, xo, by = ".tip_id", all.x = TRUE)
  out[, virus_family := family]
  out[, collapsed := FALSE]
  out[, collapsed_n := 1L]

  breadth_cov <- mean(is.finite(out$host_family_breadth) & out$host_family_breadth > 0)
  order_cov <- mean(!is.na(out$host_order) & nzchar(out$host_order))

  if (!is.finite(breadth_cov) || breadth_cov <= 0) {
    stop("Host-family breadth mapped to zero tree tips for ", family)
  }
  if (!is.finite(order_cov) || order_cov <= 0) {
    stop("Host-order annotation mapped to zero tree tips for ", family)
  }
  if (breadth_cov < 0.50) warning(family, ": host-family breadth coverage is only ", round(100 * breadth_cov, 1), "%")
  if (order_cov < 0.50) warning(family, ": host-order annotation coverage is only ", round(100 * order_cov, 1), "%")

  audit <- data.table(
    virus_family = family,
    tree_tip_n = length(tips),
    breadth_source = br$path,
    breadth_id_col = br$id_col,
    breadth_col = br$trait_col,
    breadth_overlap_n = br$overlap_n,
    breadth_coverage = breadth_cov,
    order_source = od$path,
    order_id_col = od$id_col,
    order_col = od$trait_col,
    order_method = od$method,
    order_overlap_n = od$overlap_n,
    order_coverage = order_cov
  )

  list(trait = out, audit = audit)
}

find_pair_low_col <- function(nms) {
  find_col(
    nms,
    exact = c("pair_sampling_q025", "rho_q025", "pair_q025", "pair_sampling_rho_q025", "stability_q025"),
    regex = c("pair.*q025", "rho.*q025", "stability.*q025")
  )
}

find_pair_high_col <- function(nms) {
  find_col(
    nms,
    exact = c("pair_sampling_q975", "rho_q975", "pair_q975", "pair_sampling_rho_q975", "stability_q975"),
    regex = c("pair.*q975", "rho.*q975", "stability.*q975")
  )
}

find_signif_col <- function(nms) {
  find_col(
    nms,
    exact = c("support_fdr_significant", "support_stratified_fdr_significant", "fdr_significant"),
    regex = c("support.*fdr.*sign", "fdr.*sign"),
    required = FALSE
  )
}

find_support_q_col <- function(nms) {
  find_col(
    nms,
    exact = c("support_stratified_bh_q", "support_stratified_fdr", "support_fdr", "bh_q_support", "support_q_value"),
    regex = c("support.*(bh.*q|fdr|q)", "bh.*q"),
    required = FALSE
  )
}

logicalize_signif <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y", "significant", "fdr significant")
}

make_clean_hist_density_plot <- function(null_dt, obs_beta) {
  dens_obj <- count_scaled_density_explicit(null_dt$.beta, bins = 30)
  xr <- range(c(null_dt$.beta, obs_beta), finite = TRUE)
  xs <- diff(xr)
  ggplot(null_dt, aes(x = .beta)) +
    geom_histogram(
      bins = 30,
      fill = scales::alpha(COL_GREY, 0.70),
      colour = "white",
      linewidth = 0.22
    ) +
    geom_line(
      data = dens_obj$df,
      aes(x = x, y = y_count_scaled),
      inherit.aes = FALSE,
      colour = COL_MID,
      linewidth = 0.72
    ) +
    geom_vline(xintercept = obs_beta, colour = COL_ORANGE, linewidth = 0.92) +
    annotate(
      "text",
      x = obs_beta + 0.02 * xs,
      y = Inf,
      label = paste0("Observed β = ", formatC(obs_beta, digits = 3, format = "f")),
      hjust = 0,
      vjust = 1.35,
      family = FONT_FAMILY,
      size = 2.45,
      colour = COL_DARK
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.09))) +
    labs(x = "Logistic decay coefficient (β)", y = "Number of permutations") +
    PUB_THEME()
}


# ============================================================
# v3.3 exact final-tree trait mapping
# ============================================================
# The final 03K tree-tip host-spectrum tables contain:
# representative_sequence_id, host_taxids, host_family_n, host_families.
# They do NOT contain host_order. For visual annotation only, v3.3 maps the
# frozen per-tip host_taxids to host_order using the frozen Stage-20 association
# master. No model, permutation, clustering, or tree inference is rerun.

resolve_stage20_association_master_v32 <- function(fig3_root) {
  final_root <- dirname(dirname(normalizePath(fig3_root, winslash = "/", mustWork = TRUE)))

  # On the cluster the Stage-20 master lived under the pipeline tree. The release
  # ships the same taxid -> host_order mapping inside the analysis input, which is
  # the association table every analysis module reads. Only two columns are used
  # here, for tip colour annotation; no model value comes from this file.
  candidates <- c(
    file.path(final_root, "final_three_source_pipeline_v1", "results",
              "20_three_source_association_master_FINAL.tsv"),
    file.path(final_root, "data", "analysis_data",
              "20_three_source_association_master_FINAL.tsv"),
    file.path(final_root, "data", "analysis_data",
              "20_three_source_association_master_FINAL.tsv.gz"),
    file.path(final_root, "data", "analysis_input",
              "association_analysis_input.tsv.gz")
  )

  for (f in candidates) {
    if (file.exists(f) && file.info(f)$size > 0) {
      return(normalizePath(f, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Frozen Stage-20 association master not found. Looked for:\n  ",
    paste(candidates, collapse = "\n  ")
  )
}

read_stage20_host_taxid_order_lookup_v32 <- function(path) {
  x <- fread(
    path,
    sep = "\t",
    select = c("host_taxid", "host_order"),
    colClasses = "character",
    showProgress = FALSE
  )
  x[, host_taxid := trimws(host_taxid)]
  x[, host_order := trimws(host_order)]
  x[host_taxid %in% c("", "NA", "N/A", "."), host_taxid := NA_character_]
  x[host_order %in% c("", "NA", "N/A", "."), host_order := NA_character_]
  x <- x[!is.na(host_taxid) & !is.na(host_order)]

  audit <- x[, .(
    order_n = uniqueN(host_order),
    order_set = paste(sort(unique(host_order)), collapse = ";")
  ), by = host_taxid]

  conflict <- audit[order_n > 1L]
  if (nrow(conflict)) {
    warning(
      "Frozen Stage-20 master contains ", nrow(conflict),
      " host TaxIDs mapping to >1 host order. These TaxIDs will be excluded from the colour lookup."
    )
  }

  good <- audit[order_n == 1L, .(host_taxid, host_order = order_set)]
  list(lookup = good, conflict = conflict)
}

resolve_final_tree_host_spectrum_v32 <- function(fig3_root, family) {
  root <- file.path(fig3_root, "pipeline_results", "03K_family_specific_markers")

  # Frozen 03K-D manifest uses the Rep-conserved Geminiviridae branch, whereas
  # the other 11 families use the final 06_AA95_alignment_QC branch. Resolve
  # these paths explicitly so superseded intermediate files cannot be selected.
  if (identical(family, "Geminiviridae")) {
    cluster_path <- file.path(
      root,
      "07_Geminiviridae_Rep_conserved",
      "Geminiviridae.tree_tip_host_spectrum.tsv"
    )
  } else {
    cluster_path <- file.path(
      root,
      "06_AA95_alignment_QC",
      family,
      paste0(family, ".tree_tip_host_spectrum.tsv")
    )
  }

  # The release ships the twelve-family tree-tip tables under data/, which is
  # where the analysis scripts read them from; the 03K pipeline directory keeps
  # the alignment QC but not these derived tip tables. Try the cluster layout
  # first so provenance is unchanged when it is present, then the released one.
  final_root <- dirname(dirname(normalizePath(fig3_root, winslash = "/", mustWork = TRUE)))
  released_path <- file.path(
    final_root, "data", "analysis_data", "phylogeny",
    "twelve_family", "tree_tip_host_spectrum",
    paste0(family, ".tree_tip_host_spectrum.tsv")
  )

  for (target in c(cluster_path, released_path)) {
    if (file.exists(target) && file.info(target)$size > 0) {
      return(normalizePath(target, winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Frozen final tree-tip host-spectrum file is missing for ", family, ". Looked for:\n  ",
    cluster_path, "\n  ", released_path
  )
}

split_taxid_set_v32 <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) return(character())
  z <- unlist(strsplit(as.character(x), "[;|,]", perl = TRUE), use.names = FALSE)
  z <- trimws(z)
  unique(z[nzchar(z) & !z %in% c("NA", "N/A", ".")])
}

prepare_final12_family_traits_v32 <- function(
    fig3_root,
    family,
    host_taxid_order_lookup,
    tree_file = NULL
) {
  if (is.null(tree_file)) tree_file <- resolve_final12_tree_file(fig3_root, family)
  tr <- ape::read.tree(tree_file)
  tips <- as.character(tr$tip.label)

  hs_file <- resolve_final_tree_host_spectrum_v32(fig3_root, family)
  hs <- fread(hs_file, sep = "\t", showProgress = FALSE)

  required <- c(
    "representative_sequence_id",
    "host_annotated",
    "host_taxid_n",
    "host_taxids",
    "host_family_n",
    "host_families"
  )
  miss <- setdiff(required, names(hs))
  if (length(miss)) {
    stop(
      family, " final tree-tip host spectrum lacks required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  hs[, .tip_id := as.character(representative_sequence_id)]
  hs[, host_family_breadth := suppressWarnings(as.numeric(host_family_n))]

  # Exact final-tree overlap: tree tips are AA95 representative sequence IDs.
  overlap_n <- length(intersect(tips, hs$.tip_id))
  if (overlap_n == 0L) {
    stop(family, ": final host-spectrum representative_sequence_id has zero overlap with final tree tips.")
  }

  ord_lookup <- host_taxid_order_lookup
  setkey(ord_lookup, host_taxid)

  hs[, host_order := vapply(host_taxids, function(z) {
    ids <- split_taxid_set_v32(z)
    if (!length(ids)) return(NA_character_)
    ord <- ord_lookup[J(ids), host_order, nomatch = 0L]
    ord <- sort(unique(ord[!is.na(ord) & nzchar(ord)]))
    if (!length(ord)) return(NA_character_)
    if (length(ord) == 1L) return(ord[1])
    "Multiple host orders"
  }, character(1))]

  genus_col <- if ("virus_genus_representative" %in% names(hs)) {
    "virus_genus_representative"
  } else {
    NA_character_
  }

  out <- data.table(.tip_id = tips)
  keep_cols <- c(".tip_id", "host_family_breadth", "host_order")
  if (!is.na(genus_col)) {
    hs[, virus_genus := as.character(get(genus_col))]
    keep_cols <- c(keep_cols, "virus_genus")
  }
  out <- merge(out, unique(hs[, ..keep_cols], by = ".tip_id"), by = ".tip_id", all.x = TRUE)
  if (!"virus_genus" %in% names(out)) out[, virus_genus := NA_character_]
  out[, virus_family := family]
  out[, collapsed := FALSE]
  out[, collapsed_n := 1L]

  breadth_cov <- mean(is.finite(out$host_family_breadth) & out$host_family_breadth > 0)
  order_cov <- mean(!is.na(out$host_order) & nzchar(out$host_order))
  multi_order_n <- sum(out$host_order == "Multiple host orders", na.rm = TRUE)

  if (!is.finite(breadth_cov) || breadth_cov <= 0) stop(family, ": zero host-family-breadth coverage.")
  if (!is.finite(order_cov) || order_cov <= 0) stop(family, ": zero host-order coverage after frozen TaxID mapping.")

  audit <- data.table(
    virus_family = family,
    tree_tip_n = length(tips),
    host_spectrum_file = hs_file,
    spectrum_row_n = nrow(hs),
    tree_tip_overlap_n = overlap_n,
    tree_tip_overlap_fraction = overlap_n / max(1L, length(tips)),
    breadth_source_col = "host_family_n",
    breadth_coverage = breadth_cov,
    order_source = "host_taxids -> frozen Stage20 host_taxid:host_order lookup",
    order_coverage = order_cov,
    multiple_host_orders_tip_n = multi_order_n,
    single_order_tip_n = sum(!is.na(out$host_order) & out$host_order != "Multiple host orders"),
    unannotated_order_tip_n = sum(is.na(out$host_order) | !nzchar(out$host_order))
  )

  list(trait = out, audit = audit)
}
