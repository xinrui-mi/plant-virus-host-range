suppressPackageStartupMessages(library(ggplot2))

resolve_font_family <- function(request = Sys.getenv("FIG_FONT", unset = "Open Sans")) {
  if (nzchar(Sys.which("fc-match"))) {
    cmd <- sprintf("fc-match -f '%%{family[0]}\\n' %s 2>/dev/null", shQuote(request))
    ans <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
    ans <- trimws(ans)
    if (length(ans) && nzchar(ans[1])) return(ans[1])
  }
  "sans"
}

REQUESTED_FONT <- Sys.getenv("FIG_FONT", unset = "Open Sans")
FONT_FAMILY <- resolve_font_family(REQUESTED_FONT)

if (
  identical(tolower(trimws(REQUESTED_FONT)), "open sans") &&
  !identical(tolower(trimws(FONT_FAMILY)), "open sans")
) {
  stop(
    "Open Sans was requested but fontconfig resolved '", FONT_FAMILY,
    "'. Run fc-match \"Open Sans\" before plotting."
  )
}

# Frozen Fig.3 typography scale.
TEXT_TINY   <- 5.7
TEXT_SMALL  <- 6.8
TEXT_NORMAL <- 8.2
TEXT_PANEL  <- 11.2

# Fig.3 shared colors.
COL_BLUE   <- "#3182BD"
COL_ORANGE <- "#E69F00"
COL_GREY   <- "#D0D0D0"
COL_LIGHT  <- "#F2F2F2"
COL_DARK   <- "#222222"
COL_MID    <- "#6F6F6F"

# Legacy Fig.5 A colors, intentionally retained exactly.
COL_A_VIRUS <- "#91CF60"
COL_A_HOST  <- "#80B1D3"

# Legacy Fig.5 B gradient, intentionally retained exactly.
COL_B_LOW  <- "#F7FBFF"
COL_B_HIGH <- "#08519C"

PUB_THEME <- function(base_size = TEXT_NORMAL, base_family = FONT_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(size = TEXT_NORMAL, face = "bold", colour = "black"),
      axis.text = element_text(size = TEXT_SMALL, face = "plain", colour = "black"),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.42, colour = "black"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      legend.title = element_text(size = TEXT_SMALL, face = "bold", colour = "black"),
      legend.text = element_text(size = TEXT_SMALL, face = "plain", colour = "black"),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = TEXT_NORMAL, face = "plain", colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

# Legacy Fig.5 B/C look, but with frozen Fig.3 typography.
LEGACY_FIG5_THEME <- function(base_size = TEXT_NORMAL, base_family = FONT_FAMILY) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.24, colour = "grey88"),
      panel.border = element_rect(colour = "black", linewidth = 0.38, fill = NA),
      axis.title = element_text(size = TEXT_NORMAL, face = "bold", colour = "black"),
      axis.text = element_text(size = TEXT_SMALL, face = "plain", colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_rect(fill = "grey94", colour = "grey45", linewidth = 0.28),
      strip.text = element_text(size = TEXT_SMALL, face = "bold", colour = "black"),
      legend.title = element_text(size = TEXT_SMALL, face = "bold", colour = "black"),
      legend.text = element_text(size = TEXT_SMALL, face = "plain", colour = "black"),
      legend.key.height = grid::unit(0.26, "cm"),
      legend.key.width = grid::unit(0.34, "cm"),
      legend.spacing.y = grid::unit(0.05, "cm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}
