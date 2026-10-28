`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

snake_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

read_tsv_locked <- function(filename) {
  p <- file.path(TABLE_DIR, filename)
  if (!file.exists(p)) stop("Missing frozen Fig.2 input: ", p)
  x <- readr::read_tsv(p, show_col_types = FALSE, progress = FALSE)
  names(x) <- make.unique(snake_name(names(x)), sep = "_")
  x
}

read_master_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required master table: ", path)
  x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  names(x) <- make.unique(snake_name(names(x)), sep = "_")
  x
}


find_col <- function(nms, exact = character(), regex = character(), required = TRUE) {
  low <- tolower(nms)
  for (z in exact) {
    hit <- which(low == tolower(z))
    if (length(hit)) return(nms[hit[1]])
  }
  for (r in regex) {
    hit <- grep(r, low, perl = TRUE)
    if (length(hit)) return(nms[hit[1]])
  }
  if (required) stop("Could not identify required column. Available: ", paste(nms, collapse = ", "))
  NULL
}

num <- function(x) suppressWarnings(as.numeric(x))

mode_or_na <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(NA_character_)
  tt <- sort(table(x), decreasing = TRUE)
  names(tt)[1]
}

is_unclassified_taxon <- function(x) {
  z <- tolower(trimws(as.character(x)))
  is.na(z) | !nzchar(z) | grepl("unclassified|unassigned|unknown|uncultured|^na$|^n/a$", z)
}

normalise_virus_type <- function(x) {
  z <- toupper(as.character(x))
  out <- ifelse(grepl("DNA", z), "DNA", ifelse(grepl("RNA", z), "RNA", NA_character_))
  out
}

standardize_virus_family <- function(x) {
  nms <- names(x)
  fam <- find_col(nms, c("virus_family", "family"), c("virus.*family"))
  cnt <- find_col(nms,
    c("qc_pass_vc80_n", "qc_pass_vc80_count", "vc80_n", "vc80_count", "unique_vc80_n", "final_vc80_clusters", "vc80_clusters"),
    c("qc.*vc80.*(n|count)", "vc80.*(n|count|cluster)", "final.*vc80"))
  typ <- find_col(nms,
    c("virus_type", "nucleic_acid_type", "nucleic_acid", "genome_type", "genome_group"),
    c("virus.*type", "nucleic.*acid", "genome.*type"), required = FALSE)

  dplyr::tibble(
    virus_family = as.character(x[[fam]]),
    vc80_n = num(x[[cnt]]),
    virus_type = if (is.null(typ)) NA_character_ else normalise_virus_type(x[[typ]])
  ) |>
    dplyr::filter(!is_unclassified_taxon(virus_family), is.finite(vc80_n), vc80_n >= 0) |>
    dplyr::group_by(virus_family) |>
    dplyr::summarise(
      vc80_n = sum(vc80_n, na.rm = TRUE),
      virus_type = mode_or_na(virus_type),
      .groups = "drop"
    )
}

derive_family_type_lookup <- function(seq_master) {
  nms <- names(seq_master)
  fam <- find_col(nms, c("virus_family", "family"), c("virus.*family"))
  typ <- find_col(nms,
    c("virus_type", "nucleic_acid_type", "nucleic_acid", "genome_type", "genome_group"),
    c("virus.*type", "nucleic.*acid", "genome.*type"), required = FALSE)
  gen <- find_col(nms,
    c("genome", "genome_type", "genome_group"),
    c("^genome$", "genome.*(type|group)"), required = FALSE)
  if (is.null(typ) && is.null(gen)) stop("Could not determine genome/virus type columns from sequence master.")
  dplyr::tibble(
    virus_family = as.character(seq_master[[fam]]),
    virus_type = normalise_virus_type(if (!is.null(typ)) seq_master[[typ]] else seq_master[[gen]])
  ) |>
    dplyr::filter(!is_unclassified_taxon(virus_family), !is.na(virus_type)) |>
    dplyr::group_by(virus_family) |>
    dplyr::summarise(virus_type = mode_or_na(virus_type), .groups = "drop")
}

standardize_host_family <- function(x) {
  nms <- names(x)
  fam <- find_col(nms, c("host_family"), c("host.*family"))
  cnt <- find_col(nms,
    c("association_records", "association_record_n", "association_n", "record_n", "records", "n"),
    c("association.*record", "record.*n$", "^n$"))
  dplyr::tibble(host_family = as.character(x[[fam]]), association_records = num(x[[cnt]])) |>
    dplyr::filter(!is_unclassified_taxon(host_family), is.finite(association_records), association_records >= 0) |>
    dplyr::group_by(host_family) |>
    dplyr::summarise(association_records = sum(association_records, na.rm = TRUE), .groups = "drop")
}

standardize_country <- function(x) {
  nms <- names(x)
  country <- find_col(nms, c("country", "country_name", "standardized_country"), c("^country$", "country.*name", "standard.*country"))
  iso <- find_col(nms, c("iso3c", "iso_a3", "country_iso3", "iso3"), c("iso.*3"), required = FALSE)
  cont <- find_col(nms, c("continent", "continent_final"), c("continent"), required = FALSE)
  cnt <- find_col(nms,
    c("association_records", "association_record_n", "association_n", "record_n", "records", "n"),
    c("association.*record", "record.*n$", "^n$"))
  dplyr::tibble(
    country = as.character(x[[country]]),
    iso3c = if (is.null(iso)) NA_character_ else toupper(as.character(x[[iso]])),
    continent = if (is.null(cont)) NA_character_ else as.character(x[[cont]]),
    association_records = num(x[[cnt]])
  ) |>
    dplyr::filter(!is.na(country), nzchar(country), is.finite(association_records), association_records > 0)
}

standardize_continent <- function(x) {
  nms <- names(x)
  cont <- find_col(nms, c("continent", "continent_final"), c("continent"))
  cnt <- find_col(nms,
    c("association_records", "association_record_n", "association_n", "record_n", "records", "n"),
    c("association.*record", "record.*n$", "^n$"))
  dplyr::tibble(continent = as.character(x[[cont]]), association_records = num(x[[cnt]])) |>
    dplyr::filter(!is.na(continent), nzchar(continent), is.finite(association_records), association_records > 0) |>
    dplyr::group_by(continent) |>
    dplyr::summarise(association_records = sum(association_records, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(percent = association_records / sum(association_records))
}

standardize_rarefaction <- function(x) {
  nms <- names(x)
  if (all(c("records", "virus_family_median", "host_family_median", "country_median") %in% nms)) {
    out <- dplyr::bind_rows(
      dplyr::tibble(
        dimension = "Virus family",
        sampled_records = num(x[["records"]]),
        richness_mean = num(x[["virus_family_median"]]),
        richness_low = num(x[[find_col(nms, c("virus_family_q025"), c("virus_family.*q0?25"), required = FALSE) %||% "virus_family_median"]]),
        richness_high = num(x[[find_col(nms, c("virus_family_q975"), c("virus_family.*q9?75"), required = FALSE) %||% "virus_family_median"]])
      ),
      dplyr::tibble(
        dimension = "Host family",
        sampled_records = num(x[["records"]]),
        richness_mean = num(x[["host_family_median"]]),
        richness_low = num(x[[find_col(nms, c("host_family_q025"), c("host_family.*q0?25"), required = FALSE) %||% "host_family_median"]]),
        richness_high = num(x[[find_col(nms, c("host_family_q975"), c("host_family.*q9?75"), required = FALSE) %||% "host_family_median"]])
      ),
      dplyr::tibble(
        dimension = "Country",
        sampled_records = num(x[["records"]]),
        richness_mean = num(x[["country_median"]]),
        richness_low = num(x[[find_col(nms, c("country_q025"), c("country.*q0?25"), required = FALSE) %||% "country_median"]]),
        richness_high = num(x[[find_col(nms, c("country_q975"), c("country.*q9?75"), required = FALSE) %||% "country_median"]])
      )
    )
  } else {
    dimc <- find_col(nms, c("dimension", "evidence_dimension", "category_dimension"), c("dimension"))
    sc <- find_col(nms, c("sampled_records", "sample_size", "n_records", "records_sampled", "records"), c("sample.*record", "sample.*size", "n.*record", "^records$"))
    mc <- find_col(nms, c("richness_mean", "mean_richness", "observed_richness_mean", "richness"), c("richness.*mean", "mean.*richness", "^richness$"))
    lc <- find_col(nms, c("richness_low", "richness_lower", "richness_q025", "richness_2_5", "lower"), c("richness.*(low|lower|025|2_5)", "^lower$"), required = FALSE)
    hc <- find_col(nms, c("richness_high", "richness_upper", "richness_q975", "richness_97_5", "upper"), c("richness.*(high|upper|975|97_5)", "^upper$"), required = FALSE)
    out <- dplyr::tibble(
      dimension = as.character(x[[dimc]]),
      sampled_records = num(x[[sc]]),
      richness_mean = num(x[[mc]]),
      richness_low = if (is.null(lc)) num(x[[mc]]) else num(x[[lc]]),
      richness_high = if (is.null(hc)) num(x[[mc]]) else num(x[[hc]])
    )
  }
  out$dimension <- dplyr::case_when(
    grepl("virus", out$dimension, ignore.case = TRUE) ~ "Virus family",
    grepl("host", out$dimension, ignore.case = TRUE) ~ "Host family",
    grepl("countr", out$dimension, ignore.case = TRUE) ~ "Country",
    TRUE ~ out$dimension
  )
  out |>
    dplyr::filter(dimension %in% c("Virus family", "Host family", "Country"), is.finite(sampled_records), is.finite(richness_mean))
}

make_concentration <- function(name, category, evidence) {
  d <- dplyr::tibble(category = as.character(category), evidence = num(evidence)) |>
    dplyr::filter(!is.na(category), nzchar(category), is.finite(evidence), evidence > 0) |>
    dplyr::arrange(dplyr::desc(evidence)) |>
    dplyr::mutate(
      category_fraction = dplyr::row_number() / dplyr::n(),
      evidence_fraction = cumsum(evidence) / sum(evidence),
      dimension = name
    ) |>
    dplyr::select(dimension, category_fraction, evidence_fraction)
  dplyr::bind_rows(
    dplyr::tibble(dimension = name, category_fraction = 0, evidence_fraction = 0), d
  )
}

pretty_log_breaks <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(1)
  lo <- floor(log10(min(x)))
  hi <- ceiling(log10(max(x)))
  b <- 10^(lo:hi)
  b[b >= min(x) & b <= max(x) * 1.05]
}

save_triplet <- function(plot, stem, outdir, width_mm, height_mm, dpi = DPI_OUT) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(outdir, paste0(stem, ".pdf")), plot = plot,
                  width = width_mm, height = height_mm, units = "mm",
                  device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(outdir, paste0(stem, ".png")), plot = plot,
                  width = width_mm, height = height_mm, units = "mm",
                  dpi = dpi, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(outdir, paste0(stem, ".tiff")), plot = plot,
                  width = width_mm, height = height_mm, units = "mm",
                  dpi = dpi, compression = "lzw", bg = "white", limitsize = FALSE)
}

read_tsv_optional <- function(filename) {
  p <- file.path(TABLE_DIR, filename)
  if (!file.exists(p)) return(NULL)
  x <- readr::read_tsv(p, show_col_types = FALSE, progress = FALSE)
  names(x) <- make.unique(snake_name(names(x)), sep = "_")
  x
}

normalize_source_label <- function(x) {
  z <- trimws(as.character(x))
  low <- tolower(z)
  dplyr::case_when(
    grepl("land.?flora|landflora", low) ~ "Land Flora",
    grepl("ncbi", low) & grepl("vhdb|virus.?host", low) ~ "NCBI + Virus–Host DB",
    grepl("vhdb|virus.?host", low) ~ "Virus–Host DB",
    grepl("ncbi", low) ~ "NCBI Virus",
    TRUE ~ z
  )
}

standardize_source_virus <- function(x) {
  nms <- names(x)
  fam <- find_col(nms, c("virus_family", "family"), c("virus.*family"))
  src <- find_col(nms,
    c("source", "source_dataset", "source_group", "source_class", "source_category", "source_stratum", "evidence_source", "source_final"),
    c("^source$", "source.*(group|class|category|stratum|final)", "evidence.*source"))
  cnt <- find_col(nms,
    c("qc_pass_vc80_n", "vc80_n", "vc80_count", "sequence_n", "unique_sequence_n", "association_n", "association_records"),
    c("qc.*vc80", "vc80.*(n|count)", "sequence.*n", "association.*(n|record)"))
  dplyr::tibble(
    virus_family = as.character(x[[fam]]),
    source = normalize_source_label(x[[src]]),
    evidence = num(x[[cnt]])
  ) |>
    dplyr::filter(!is_unclassified_taxon(virus_family), !is.na(source), nzchar(source), is.finite(evidence), evidence >= 0) |>
    dplyr::group_by(virus_family, source) |>
    dplyr::summarise(evidence = sum(evidence, na.rm = TRUE), .groups = "drop")
}

standardize_source_host <- function(x) {
  nms <- names(x)
  fam <- find_col(nms, c("host_family"), c("host.*family"))
  src <- find_col(nms,
    c("source", "source_dataset", "source_group", "source_class", "source_category", "source_stratum", "evidence_source", "source_final"),
    c("^source$", "source.*(group|class|category|stratum|final)", "evidence.*source"))
  cnt <- find_col(nms,
    c("association_n", "association_records", "association_record_n", "record_n", "unique_sequence_n", "sequence_n"),
    c("association.*(n|record)", "record.*n$", "sequence.*n"))
  dplyr::tibble(
    host_family = as.character(x[[fam]]),
    source = normalize_source_label(x[[src]]),
    evidence = num(x[[cnt]])
  ) |>
    dplyr::filter(!is_unclassified_taxon(host_family), !is.na(source), nzchar(source), is.finite(evidence), evidence >= 0) |>
    dplyr::group_by(host_family, source) |>
    dplyr::summarise(evidence = sum(evidence, na.rm = TRUE), .groups = "drop")
}

preferred_source_levels <- function(x) {
  pref <- c("NCBI Virus", "Virus–Host DB", "NCBI + Virus–Host DB", "Land Flora")
  c(pref[pref %in% x], sort(setdiff(unique(x), pref)))
}
