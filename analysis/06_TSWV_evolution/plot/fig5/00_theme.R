suppressPackageStartupMessages(library(ggplot2))

FONT_FAMILY <- "Open Sans"

# Frozen publication font. Fail loudly on the cluster rather than silently
# falling back to another family.
if (nzchar(Sys.which("fc-match"))) {
  resolved_font <- tryCatch(
    trimws(
      system(
        "fc-match -f '%{family[0]}\\n' 'Open Sans' 2>/dev/null",
        intern = TRUE
      )
    ),
    error = function(e) character()
  )

  if (
    !length(resolved_font) ||
    !identical(
      tolower(trimws(resolved_font[1])),
      "open sans"
    )
  ) {
    stop(
      "Open Sans is required for frozen Fig.2–7 typography, ",
      "but fontconfig resolved: ",
      if (length(resolved_font)) resolved_font[1] else "<nothing>"
    )
  }
}

# Frozen manuscript typography inherited from finalized Fig.2–6.
TEXT_TINY   <- 5.7
TEXT_SMALL  <- 6.8
TEXT_NORMAL <- 8.2
TEXT_LARGE  <- 9.2

# Compatibility aliases used in this plotting package.
BASE_FAMILY <- FONT_FAMILY
BASE_SIZE <- TEXT_NORMAL
STRIP_SIZE <- TEXT_SMALL
AXIS_TITLE_SIZE <- TEXT_NORMAL
AXIS_TEXT_SIZE <- TEXT_SMALL
LEGEND_TITLE_SIZE <- TEXT_SMALL
LEGEND_TEXT_SIZE <- TEXT_SMALL
PLOT_TITLE_SIZE <- TEXT_NORMAL
TREE_SMALL_SIZE <- TEXT_TINY

COL_BLUE <- "#3182BD"
COL_ORANGE <- "#E69F00"
COL_GREEN <- "#009E73"

COL_DARK <- "#222222"
COL_MID <- "#6F6F6F"
COL_TREE <- "#7A7A7A"
COL_SUPPORT <- "#989898"
COL_GRID <- "#ECECEC"

# GARD breakpoint red retained from the manuscript-style version.
COL_BREAKPOINT <- "#E64B35"
COL_EVENT <- "#2F2F2F"

SEGMENT_COLORS <- c(
  "L" = COL_BLUE,
  "M" = COL_ORANGE,
  "S" = COL_GREEN
)

PUB_THEME <- function(
  base_size = TEXT_SMALL,
  base_family = FONT_FAMILY
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
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

      axis.line = element_line(
        linewidth = 0.45,
        colour = "black"
      ),

      axis.ticks = element_line(
        linewidth = 0.42,
        colour = "black"
      ),

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
        size = TEXT_SMALL,
        face = "bold",
        colour = "black"
      ),

      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

TREE_THEME <- function(
  title_colour = COL_DARK
) {
  theme_void(
    base_family = FONT_FAMILY
  ) +
    theme(
      text = element_text(
        family = FONT_FAMILY,
        colour = COL_DARK
      ),

      plot.title = element_text(
        family = FONT_FAMILY,
        size = TEXT_SMALL,
        face = "bold",
        hjust = 0.5,
        colour = title_colour,
        margin = margin(0, 0, 0.4, 0)
      ),

      plot.margin = margin(1, 1, 1, 1)
    )
}
