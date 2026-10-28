suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

FIG_FONT <- Sys.getenv("FIG_FONT", unset = "Open Sans")

# EXACT typography freeze inherited from Fig4_publication_FINAL_v25/v26:
# axis titles = 8.2 pt bold
# axis tick labels = 6.8 pt plain
# legend titles = 6.8 pt bold
# legend text = 6.8 pt plain
# facet-strip text = 6.8 pt bold
# no embedded A/B/C/D panel letters
# no plot titles/subtitles/captions

theme_publication <- function() {
  theme_classic(base_size = TEXT_NORMAL, base_family = FIG_FONT) +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(size = TEXT_NORMAL, face = "bold"),
      axis.text = element_text(size = TEXT_SMALL, colour = "grey15"),
      legend.title = element_text(size = TEXT_SMALL, face = "bold"),
      legend.text = element_text(size = TEXT_SMALL),
      strip.background = element_rect(
        fill = "grey96", colour = "grey82", linewidth = LINE_THIN
      ),
      strip.text = element_text(size = TEXT_SMALL, face = "bold"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.22),
      panel.grid.major.x = element_line(colour = "grey94", linewidth = 0.20),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

save_plot_all_formats <- function(
  plot_obj, filename_base, width_mm, height_mm, out_dir, dpi = DPI_OUT
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  ggsave(
    file.path(out_dir, paste0(filename_base, ".pdf")),
    plot_obj, width = w, height = h, units = "in",
    device = grDevices::cairo_pdf,
    bg = "white", limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      file.path(out_dir, paste0(filename_base, ".png")),
      plot_obj, width = w, height = h, units = "in",
      dpi = dpi, device = ragg::agg_png,
      bg = "white", limitsize = FALSE
    )
    ggsave(
      file.path(out_dir, paste0(filename_base, ".tiff")),
      plot_obj, width = w, height = h, units = "in",
      dpi = dpi, device = ragg::agg_tiff,
      compression = "lzw", bg = "white", limitsize = FALSE
    )
  } else {
    ggsave(
      file.path(out_dir, paste0(filename_base, ".png")),
      plot_obj, width = w, height = h, units = "in",
      dpi = dpi, bg = "white", limitsize = FALSE
    )
    ggsave(
      file.path(out_dir, paste0(filename_base, ".tiff")),
      plot_obj, width = w, height = h, units = "in",
      dpi = dpi, compression = "lzw", bg = "white", limitsize = FALSE
    )
  }
}
