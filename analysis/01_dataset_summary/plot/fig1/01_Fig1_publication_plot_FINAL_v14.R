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
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
  library(ggrepel)
  library(grid)
})

log_file <- file.path(OUT_LOG, "Fig1_publication_plot_FINAL_v14.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Fig.1 plotting from frozen Fig.2 started: ", format(Sys.time()), "\n", sep = "")
cat("FIG2_ROOT: ", FIG2_ROOT, "\n", sep = "")
cat("Resolved figure font: ", FONT_FAMILY, "\n", sep = "")
cat("R version: ", R.version.string, "\n\n", sep = "")

# -----------------------------------------------------------------------------
# Inputs: frozen Fig.2 tables only. Plotting-only extension.
# -----------------------------------------------------------------------------
raw01 <- read_tsv_locked("01_virus_family_evidence.tsv")
raw02 <- read_tsv_locked("02_host_family_evidence.tsv")
raw03 <- read_tsv_locked("03_country_evidence.tsv")
raw04 <- read_tsv_locked("04_continent_evidence.tsv")
raw07 <- read_tsv_locked("07_evidence_rarefaction_summary.tsv")
seq_master <- read_master_tsv(SEQ_MASTER)

virus_all <- standardize_virus_family(raw01)
family_type_lookup <- derive_family_type_lookup(seq_master)
virus_all <- virus_all |>
  dplyr::left_join(family_type_lookup, by = "virus_family", suffix = c("", "_master")) |>
  dplyr::mutate(virus_type = dplyr::coalesce(virus_type, virus_type_master)) |>
  dplyr::select(-virus_type_master)
if (any(is.na(virus_all$virus_type))) {
  missing_fam <- paste(sort(virus_all$virus_family[is.na(virus_all$virus_type)]), collapse = ", ")
  stop("Could not infer DNA/RNA type for virus families: ", missing_fam)
}
host_all <- standardize_host_family(raw02)
country_all <- standardize_country(raw03)
continent_all <- standardize_continent(raw04)
rare_all <- standardize_rarefaction(raw07)

computed_sequence_records <- nrow(seq_master)
computed_association_records <- sum(host_all$association_records, na.rm = TRUE)
computed_classified_vc80 <- sum(virus_all$vc80_n, na.rm = TRUE)

n_sequence_records <- FIG1_SEQUENCE_RECORDS_N
n_association_records <- FIG1_ASSOCIATION_RECORDS_N
n_vc80 <- FIG1_VC80_N

cat("Summary totals for new top strip (frozen manuscript totals):\n")
cat("  sequence records: ", scales::comma(n_sequence_records), "\n", sep = "")
cat("  virus–host association records: ", scales::comma(n_association_records), "\n", sep = "")
cat("  final VC80s: ", scales::comma(n_vc80), "\n", sep = "")
cat("Input cross-checks:\n")
cat("  sequence master rows: ", scales::comma(computed_sequence_records), "\n", sep = "")
cat("  association records summed across classified host families: ", scales::comma(computed_association_records), "\n", sep = "")
cat("  VC80s summed across classified virus families: ", scales::comma(computed_classified_vc80), "\n\n", sep = "")

if (computed_sequence_records != FIG1_SEQUENCE_RECORDS_N) {
  warning("Sequence-master row count differs from frozen Fig.1A total: ", computed_sequence_records, " vs ", FIG1_SEQUENCE_RECORDS_N)
}
if (computed_association_records != FIG1_ASSOCIATION_RECORDS_N) {
  # This difference is expected, not a regression. The frozen Fig. 1A total
  # (83,034) counts every virus-host association record; the per-host-family
  # bars can only sum records whose host resolves to a Viridiplantae family.
  # Records whose host TaxID has no family assignment are therefore absent from
  # the panel sum but present in the headline total. On the released data the
  # gap is 173 records. Report the decomposition so a reader can verify it
  # rather than having to guess whether the tables are stale.
  message(
    "NOTE: classified host-family association sum (", scales::comma(computed_association_records),
    ") is below the frozen Fig.1A total (", scales::comma(FIG1_ASSOCIATION_RECORDS_N),
    ") by ", scales::comma(FIG1_ASSOCIATION_RECORDS_N - computed_association_records),
    " records whose host has no assigned family; this is expected and does not affect the panel."
  )
}

# -----------------------------------------------------------------------------
# New compact top strip: three-source data integration.
# v14: retain the established manuscript typography while reducing the
# processing description to three concise layers and keeping the title on one line.
# -----------------------------------------------------------------------------
source_boxes <- tibble::tribble(
  ~label,                  ~xmin, ~xmax, ~ymin, ~ymax,
  "NCBI Virus",             2.0,  15.2, 14.6, 19.4,
  "Virus–Host DB",          2.0,  15.2,  8.8, 13.6,
  "Land Flora\nRNA virome", 2.0,  15.2,  3.0,  7.8
)

# Three blocks retain a matched overall vertical footprint. The processing
# block remains wide enough for three concise method layers without reducing
# the frozen TEXT_TINY size. The one-line title is separated vertically from them.
harm_box <- tibble::tibble(xmin = 22.2, xmax = 59.0, ymin = 3.3, ymax = 18.7)
integrated_box <- tibble::tibble(xmin = 62.0, xmax = 98.5, ymin = 3.3, ymax = 18.7)

stats_df <- tibble::tribble(
  ~x,   ~n,                                   ~label,
  68.2, scales::comma(n_sequence_records),    "Sequence\nrecords",
  80.5, scales::comma(n_association_records), "Virus–host\nassociation records",
  93.0, scales::comma(n_vc80),                "VC80s"
)

# Arrow helper coordinates follow the vertically enlarged source boxes.
src_mid_y <- c(17.0, 11.2, 5.4)

pA <- ggplot() +
  annotate("rect", xmin = source_boxes$xmin, xmax = source_boxes$xmax,
           ymin = source_boxes$ymin, ymax = source_boxes$ymax,
           fill = COL_LIGHT, colour = COL_MID, linewidth = 0.35) +
  annotate("rect", xmin = harm_box$xmin, xmax = harm_box$xmax,
           ymin = harm_box$ymin, ymax = harm_box$ymax,
           fill = scales::alpha(COL_GREY, 0.22), colour = COL_MID, linewidth = 0.35) +
  annotate("rect", xmin = integrated_box$xmin, xmax = integrated_box$xmax,
           ymin = integrated_box$ymin, ymax = integrated_box$ymax,
           fill = scales::alpha(COL_BLUE, 0.08), colour = COL_BLUE, linewidth = 0.42) +
  annotate("text", x = (source_boxes$xmin + source_boxes$xmax) / 2,
           y = (source_boxes$ymin + source_boxes$ymax) / 2,
           label = source_boxes$label,
           family = FONT_FAMILY, size = TEXT_SMALL / ggplot2::.pt, lineheight = 0.95) +
  annotate("text", x = 40.6, y = 15.4,
           label = "Filtering, harmonization and quality control",
           family = FONT_FAMILY, fontface = "bold",
           size = TEXT_SMALL / ggplot2::.pt) +
  annotate("text", x = 40.6, y = 8.8,
           label = paste(
             "Plant host and sequence filtering · accession deduplication",
             "Taxonomy, genome type and geographic metadata harmonization",
             "VC80 clustering and quality control",
             sep = "\n"
           ),
           family = FONT_FAMILY, size = TEXT_TINY / ggplot2::.pt, lineheight = 1.05) +
  annotate("text", x = 80.3, y = 16.3, label = "Integrated evidence base",
           family = FONT_FAMILY, fontface = "bold", size = TEXT_SMALL / ggplot2::.pt) +
  annotate("text", x = stats_df$x, y = 12.2, label = stats_df$n,
           family = FONT_FAMILY, fontface = "bold", size = TEXT_NORMAL / ggplot2::.pt) +
  annotate("text", x = stats_df$x, y = 7.8, label = stats_df$label,
           family = FONT_FAMILY, size = TEXT_SMALL / ggplot2::.pt, lineheight = 0.92) +
  annotate("segment", x = 15.2, xend = 20.3, y = src_mid_y[1], yend = src_mid_y[1], colour = COL_MID, linewidth = 0.34) +
  annotate("segment", x = 15.2, xend = 20.3, y = src_mid_y[2], yend = src_mid_y[2], colour = COL_MID, linewidth = 0.34) +
  annotate("segment", x = 15.2, xend = 20.3, y = src_mid_y[3], yend = src_mid_y[3], colour = COL_MID, linewidth = 0.34) +
  annotate("segment", x = 20.3, xend = 20.3, y = src_mid_y[3], yend = src_mid_y[1], colour = COL_MID, linewidth = 0.34) +
  annotate("segment", x = 20.3, xend = 22.2, y = 11.2, yend = 11.2, colour = COL_MID, linewidth = 0.34,
           arrow = grid::arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("segment", x = 59.0, xend = 62.0, y = 11.2, yend = 11.2, colour = COL_MID, linewidth = 0.34,
           arrow = grid::arrow(length = unit(1.6, "mm"), type = "closed")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 20), expand = FALSE, clip = "off") +
  theme_void(base_family = FONT_FAMILY, base_size = TEXT_NORMAL) +
  theme(plot.margin = margin(3, 4, 2, 4))

# -----------------------------------------------------------------------------
# Inherited former Fig.2 Panel A -> new Fig.1 virus-family evidence panel.
# -----------------------------------------------------------------------------
virus_top <- virus_all |>
  arrange(desc(vc80_n), virus_family) |>
  slice_head(n = TOP_VIRUS_FAMILY_N) |>
  mutate(virus_family = factor(virus_family, levels = rev(virus_family)))

pB <- ggplot(virus_top, aes(x = virus_family, y = vc80_n, fill = virus_type)) +
  geom_col(width = 0.66) +
  geom_text(
    aes(label = comma(vc80_n)), hjust = -0.08,
    family = FONT_FAMILY, size = TEXT_TINY / ggplot2::.pt, colour = "black"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = DNA_RNA_COLS, breaks = c("DNA", "RNA"), drop = FALSE, name = "Genome type") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Number of VC80s") +
  PUB_THEME() +
  theme(
    axis.text.y = element_text(size = TEXT_SMALL, face = "italic", margin = margin(r = 1.5)),
    panel.grid.major.x = element_line(linewidth = 0.22, colour = "grey90"),
    legend.position = c(0.82, 0.18),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.84), colour = NA),
    legend.key.width = grid::unit(4.8, "mm"),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.spacing.y = grid::unit(0.6, "mm"),
    plot.margin = margin(4, 18, 4, 4)
  )

# -----------------------------------------------------------------------------
# Inherited former Fig.2 Panel B -> new Fig.1 host-family evidence panel.
# -----------------------------------------------------------------------------
host_top <- host_all |>
  arrange(desc(association_records), host_family) |>
  slice_head(n = TOP_HOST_FAMILY_N) |>
  mutate(host_family = factor(host_family, levels = rev(host_family)))

pC <- ggplot(host_top, aes(x = host_family, y = association_records)) +
  geom_col(width = 0.66, fill = "grey38") +
  geom_text(
    aes(label = comma(association_records)), hjust = -0.07,
    family = FONT_FAMILY, size = TEXT_TINY / ggplot2::.pt, colour = "black"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "Virus–host association records") +
  PUB_THEME() +
  theme(
    axis.text.y = element_text(size = TEXT_SMALL, face = "plain"),
    panel.grid.major.x = element_line(linewidth = 0.22, colour = "grey90"),
    plot.margin = margin(4, 28, 4, 4)
  )

# -----------------------------------------------------------------------------
# Inherited former Fig.2 map panel.
# -----------------------------------------------------------------------------
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  dplyr::mutate(
    join_iso3 = dplyr::case_when(
      !is.na(iso_a3) & iso_a3 != "-99" & nzchar(iso_a3) ~ iso_a3,
      !is.na(adm0_a3) & adm0_a3 != "-99" & nzchar(adm0_a3) ~ adm0_a3,
      TRUE ~ NA_character_
    )
  )

normalize_country_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(trimws(x))
  x <- gsub("&", " and ", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

world_attr <- world_sf |> sf::st_drop_geometry()
match_fields <- intersect(c("admin", "name", "name_long", "geounit", "subunit", "formal_en"), names(world_attr))
if (!length(match_fields)) stop("No Natural Earth country-name fields available for Fig.1 map matching.")

country_map <- country_all |>
  dplyr::mutate(
    name_key = normalize_country_key(country),
    map_iso3 = NA_character_,
    matched_world_name = NA_character_,
    match_field = NA_character_,
    match_method = NA_character_
  )

valid_world_iso3 <- unique(world_sf$join_iso3[!is.na(world_sf$join_iso3)])
has_valid_iso <- !is.na(country_map$iso3c) & nzchar(country_map$iso3c) & country_map$iso3c %in% valid_world_iso3
if (any(has_valid_iso)) {
  country_map$map_iso3[has_valid_iso] <- country_map$iso3c[has_valid_iso]
  country_map$match_method[has_valid_iso] <- "ISO3"
  country_map$match_field[has_valid_iso] <- "ISO3"
  country_map$matched_world_name[has_valid_iso] <- country_map$country[has_valid_iso]
}

for (fld in match_fields) {
  lookup <- dplyr::tibble(
    matched_world_name = as.character(world_attr[[fld]]),
    join_iso3 = as.character(world_attr$join_iso3)
  ) |>
    dplyr::mutate(name_key = normalize_country_key(matched_world_name)) |>
    dplyr::filter(
      !is.na(join_iso3), nzchar(join_iso3),
      !is.na(name_key), nzchar(name_key),
      !is.na(matched_world_name), nzchar(matched_world_name)
    ) |>
    dplyr::distinct(name_key, join_iso3, .keep_all = TRUE)

  unique_keys <- lookup |>
    dplyr::group_by(name_key) |>
    dplyr::summarise(n_iso = dplyr::n_distinct(join_iso3), .groups = "drop") |>
    dplyr::filter(n_iso == 1) |>
    dplyr::pull(name_key)

  lookup <- lookup |>
    dplyr::filter(name_key %in% unique_keys) |>
    dplyr::arrange(name_key) |>
    dplyr::group_by(name_key) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  mm <- match(country_map$name_key, lookup$name_key)
  fill <- is.na(country_map$map_iso3) & !is.na(mm)
  if (any(fill)) {
    country_map$map_iso3[fill] <- lookup$join_iso3[mm[fill]]
    country_map$matched_world_name[fill] <- lookup$matched_world_name[mm[fill]]
    country_map$match_field[fill] <- fld
    country_map$match_method[fill] <- paste0("NaturalEarth:", fld)
  }
}

map_join_qc <- country_map |>
  dplyr::mutate(matched_world = !is.na(map_iso3) & nzchar(map_iso3)) |>
  dplyr::select(country, continent, association_records, iso3c,
                map_iso3, matched_world_name, match_field, match_method, matched_world) |>
  dplyr::arrange(dplyr::desc(association_records))
readr::write_tsv(map_join_qc, file.path(OUT_LOG, "Fig1_map_join_QC.tsv"))

unmatched_map <- map_join_qc |>
  dplyr::filter(!matched_world) |>
  dplyr::arrange(dplyr::desc(association_records))
readr::write_tsv(unmatched_map, file.path(OUT_LOG, "Fig1_map_unmatched.tsv"))

high_unmatched <- unmatched_map |> dplyr::filter(association_records >= 500)
if (nrow(high_unmatched) > 0) {
  stop(
    "High-evidence countries remain unmatched in Natural Earth map join: ",
    paste0(high_unmatched$country, " (", scales::comma(high_unmatched$association_records), ")", collapse = "; "),
    ". See Fig1_map_unmatched.tsv."
  )
}

country_map_for_join <- country_map |>
  dplyr::filter(!is.na(map_iso3), nzchar(map_iso3)) |>
  dplyr::group_by(map_iso3) |>
  dplyr::summarise(
    country = paste(unique(country), collapse = "; "),
    association_records = sum(association_records, na.rm = TRUE),
    .groups = "drop"
  )

world_plot <- world_sf |>
  dplyr::left_join(country_map_for_join, by = c("join_iso3" = "map_iso3"))

country_label_poly <- world_plot |>
  dplyr::filter(!is.na(association_records), association_records > 0) |>
  dplyr::arrange(dplyr::desc(association_records)) |>
  dplyr::slice_head(n = TOP_COUNTRY_LABEL_N)

country_label_pts <- suppressWarnings(
  country_label_poly |>
    sf::st_transform(3857) |>
    sf::st_point_on_surface() |>
    sf::st_transform(4326)
)

country_label_coords <- sf::st_coordinates(country_label_pts)
country_label_xy <- country_label_pts |>
  sf::st_drop_geometry() |>
  dplyr::mutate(
    lon = country_label_coords[, 1],
    lat = country_label_coords[, 2],
    country_label = country,
    label = paste0(country_label, "\n", scales::comma(association_records))
  )

map_breaks <- pretty_log_breaks(country_all$association_records)

pD_map <- ggplot() +
  geom_sf(data = world_plot, aes(fill = association_records), colour = "grey88", linewidth = 0.10) +
  ggrepel::geom_text_repel(
    data = country_label_xy,
    aes(x = lon, y = lat, label = label),
    family = FONT_FAMILY, size = TEXT_TINY / ggplot2::.pt, fontface = "bold",
    colour = "white", bg.color = scales::alpha("black", 0.55), bg.r = 0.12,
    min.segment.length = 0, segment.size = 0.18, segment.color = "grey35",
    box.padding = 0.26, point.padding = 0.10, max.overlaps = Inf, seed = 20260818
  ) +
  scale_fill_gradientn(
    colours = MAP_COLS, trans = "log10", breaks = map_breaks,
    labels = comma, na.value = "grey96", name = "Association records"
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  theme_void(base_family = FONT_FAMILY, base_size = TEXT_NORMAL) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = TEXT_SMALL, face = "bold"),
    legend.text = element_text(size = TEXT_SMALL),
    legend.key.width = grid::unit(9.2, "mm"),
    legend.key.height = grid::unit(2.4, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  )

continent_plot <- continent_all |>
  dplyr::filter(
    !grepl("antarctica", continent, ignore.case = TRUE),
    percent > 0,
    continent %in% c("Asia", "Americas", "Europe", "Africa", "Oceania")
  ) |>
  dplyr::group_by(continent) |>
  dplyr::summarise(
    association_records = sum(association_records),
    percent = sum(percent),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    continent = factor(continent, levels = c("Asia", "Americas", "Europe", "Africa", "Oceania"))
  )
cont_cols <- CONTINENT_COLS[levels(continent_plot$continent)]
cont_cols[is.na(cont_cols)] <- COL_GREY

pD_inset <- ggplot(continent_plot, aes(x = continent, y = percent, fill = continent)) +
  geom_col(width = 0.60) +
  geom_text(
    aes(y = percent + 0.018, label = percent(percent, accuracy = 0.1)),
    family = FONT_FAMILY, size = TEXT_TINY / ggplot2::.pt, vjust = 0
  ) +
  scale_fill_manual(values = cont_cols, guide = "none", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 0.65), breaks = c(0, 0.20, 0.40, 0.60),
    labels = percent_format(accuracy = 10), expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = FONT_FAMILY, base_size = TEXT_TINY) +
  theme(
    axis.text.x = element_text(size = TEXT_TINY, angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = TEXT_TINY),
    axis.line = element_line(linewidth = 0.24),
    axis.ticks = element_line(linewidth = 0.24),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(2, 0, 0, 0)
  )

pD <- pD_map + patchwork::inset_element(
  pD_inset, left = -0.010, bottom = 0.012, right = 0.225, top = 0.355, align_to = "panel"
)

# -----------------------------------------------------------------------------
# Inherited evidence concentration panel.
# -----------------------------------------------------------------------------
concentration <- bind_rows(
  make_concentration("Virus family", virus_all$virus_family, virus_all$vc80_n),
  make_concentration("Host family", host_all$host_family, host_all$association_records),
  make_concentration("Country", country_all$country, country_all$association_records)
) |>
  mutate(dimension = factor(dimension, levels = c("Virus family", "Host family", "Country")))

pE <- ggplot(concentration, aes(x = category_fraction, y = evidence_fraction, colour = dimension)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey78", linewidth = 0.45) +
  geom_line(linewidth = 0.90, lineend = "round") +
  scale_colour_manual(values = DIM_COLS, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 25), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.035))) +
  scale_y_continuous(labels = percent_format(accuracy = 25), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.025))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(x = "Ranked category share (%)", y = "Cumulative evidence share (%)") +
  PUB_THEME() +
  theme(
    panel.grid.major = element_line(linewidth = 0.22, colour = "grey90"),
    legend.position = c(0.72, 0.18),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.82), colour = NA),
    legend.key.width = grid::unit(3.8, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.25, "mm"),
    legend.box.spacing = grid::unit(0, "mm"),
    plot.margin = margin(4, 8, 4, 5)
  )

# -----------------------------------------------------------------------------
# Inherited evidence accumulation / rarefaction panel.
# -----------------------------------------------------------------------------
rare_all <- rare_all |>
  mutate(dimension = factor(dimension, levels = c("Virus family", "Host family", "Country")))

pF <- ggplot(rare_all, aes(x = sampled_records, y = richness_mean)) +
  geom_ribbon(aes(ymin = richness_low, ymax = richness_high), fill = COL_GREY, alpha = 0.40, colour = NA) +
  geom_line(linewidth = 0.78, colour = COL_BLUE) +
  facet_wrap(~dimension, scales = "free_y", nrow = 1) +
  scale_x_continuous(labels = comma, breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Subsampled association records", y = "Observed richness") +
  PUB_THEME() +
  theme(
    panel.grid.major = element_line(linewidth = 0.22, colour = "grey90"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = TEXT_SMALL, face = "bold"),
    plot.margin = margin(3, 4, 3, 4)
  )

# -----------------------------------------------------------------------------
# Composite layout. Existing B-F panels keep their inherited physical sizes.
# -----------------------------------------------------------------------------
row_bc <- patchwork::wrap_plots(
  list(patchwork::wrap_elements(full = pB), patchwork::wrap_elements(full = pC)),
  ncol = 2, widths = BC_WIDTHS_MM
)
row_ef <- patchwork::wrap_plots(
  list(patchwork::wrap_elements(full = pE), patchwork::wrap_elements(full = pF)),
  ncol = 2, widths = BOTTOM_WIDTHS_MM
)

main_fig <- patchwork::wrap_plots(
  list(
    patchwork::wrap_elements(full = pA),
    patchwork::wrap_elements(full = row_bc),
    patchwork::wrap_elements(full = pD),
    patchwork::wrap_elements(full = row_ef)
  ),
  ncol = 1,
  heights = c(A_HEIGHT_MM, BC_HEIGHT_MM, MAP_HEIGHT_MM, BOTTOM_HEIGHT_MM)
)

cat(
  "Main/panel geometry (mm): main=", MAIN_WIDTH_MM, "x", MAIN_HEIGHT_MM,
  "; top strip=", MAIN_WIDTH_MM, "x", A_HEIGHT_MM,
  "; B/C=", paste(BC_WIDTHS_MM, collapse = "/"), "x", BC_HEIGHT_MM,
  "; D=", MAP_WIDTH_MM, "x", MAP_HEIGHT_MM,
  "; E/F=", paste(BOTTOM_WIDTHS_MM, collapse = "/"), "x", BOTTOM_HEIGHT_MM, "\n", sep = ""
)

save_triplet(main_fig, "Fig1_main_publication_FINAL_v14", OUT_MAIN, MAIN_WIDTH_MM, MAIN_HEIGHT_MM)
save_triplet(pA, "Fig1_data_integration_strip", OUT_PANEL, MAIN_WIDTH_MM, A_HEIGHT_MM)
save_triplet(pB, "Fig1_virus_family_evidence", OUT_PANEL, BC_WIDTHS_MM[1], BC_HEIGHT_MM)
save_triplet(pC, "Fig1_host_family_evidence", OUT_PANEL, BC_WIDTHS_MM[2], BC_HEIGHT_MM)
save_triplet(pD, "Fig1_global_geographic_evidence", OUT_PANEL, MAP_WIDTH_MM, MAP_HEIGHT_MM)
save_triplet(pE, "Fig1_evidence_concentration", OUT_PANEL, BOTTOM_WIDTHS_MM[1], BOTTOM_HEIGHT_MM)
save_triplet(pF, "Fig1_evidence_accumulation", OUT_PANEL, BOTTOM_WIDTHS_MM[2], BOTTOM_HEIGHT_MM)

cat("\nFig.1 plotting from frozen Fig.2 completed: ", format(Sys.time()), "\n", sep = "")
