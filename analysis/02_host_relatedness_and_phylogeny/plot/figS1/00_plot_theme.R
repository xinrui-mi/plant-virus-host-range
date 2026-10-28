suppressPackageStartupMessages(library(ggplot2))

resolve_font_family <- function(request = Sys.getenv("FIG_FONT", unset = "Open Sans")) {
  # Resolve the explicitly requested publication font through fontconfig.
  if (nzchar(Sys.which("fc-match"))) {
    cmd <- sprintf("fc-match -f '%%{family[0]}\\n' %s 2>/dev/null", shQuote(request))
    ans <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
    ans <- trimws(ans)
    if (length(ans) && nzchar(ans[1])) return(ans[1])
  }
  "sans"
}

REQUESTED_FONT <- Sys.getenv(
  "FIG_FONT",
  unset = "Open Sans"
)

FONT_FAMILY <- resolve_font_family(
  REQUESTED_FONT
)

if (
  identical(
    tolower(trimws(REQUESTED_FONT)),
    "open sans"
  ) &&
  !identical(
    tolower(trimws(FONT_FAMILY)),
    "open sans"
  )
) {
  stop(
    "Open Sans was requested but fontconfig resolved '",
    FONT_FAMILY,
    "'. Run fc-match \"Open Sans\" before plotting."
  )
}

TEXT_TINY   <- 5.7
TEXT_SMALL  <- 6.8
TEXT_NORMAL <- 8.2
TEXT_PANEL  <- 11.2

COL_BLUE   <- "#3182BD"
COL_ORANGE <- "#E69F00"
COL_GREY   <- "#D0D0D0"
COL_LIGHT  <- "#F2F2F2"
COL_DARK   <- "#222222"
COL_MID    <- "#6F6F6F"

FAMILY_COLORS <- c(
  "Geminiviridae"  = "#3182BD",
  "Caulimoviridae" = "#E69F00",
  "Potyviridae"    = "#009E73",
  "Tombusviridae"  = "#CC79A7"
)

PUB_THEME <- function(base_size = TEXT_NORMAL, base_family = FONT_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),

      axis.title = element_text(
        size = TEXT_NORMAL,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = TEXT_SMALL,
        face = "plain",
        colour = "black"
      ),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.42, colour = "black"),
      axis.ticks.length = grid::unit(1.4, "mm"),

      legend.title = element_text(
        size = TEXT_SMALL,
        face = "bold",
        colour = "black"
      ),
      legend.text = element_text(
        size = TEXT_SMALL,
        face = "plain",
        colour = "black"
      ),
      legend.key = element_blank(),

      strip.background = element_blank(),
      strip.text = element_text(
        size = TEXT_NORMAL,
        face = "plain",
        colour = "black"
      ),

      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

PANEL_TAG_THEME <- theme(
  plot.tag = element_text(
    size = TEXT_PANEL,
    face = "bold",
    family = FONT_FAMILY,
    colour = "black"
  ),
  plot.tag.position = c(0.005, 0.995)
)
