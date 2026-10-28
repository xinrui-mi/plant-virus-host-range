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
if (identical(tolower(trimws(REQUESTED_FONT)), "open sans") &&
    !identical(tolower(trimws(FONT_FAMILY)), "open sans")) {
  stop("Open Sans was requested but fontconfig resolved '", FONT_FAMILY,
       "'. Run: fc-match \"Open Sans\"")
}

# Frozen manuscript typography inherited from finalized Fig.3/Fig.5.
TEXT_TINY   <- 5.7
TEXT_SMALL  <- 6.8
TEXT_NORMAL <- 8.2
TEXT_PANEL  <- 11.2

# Shared manuscript colours retained where appropriate.
COL_BLUE   <- "#3182BD"
COL_ORANGE <- "#E69F00"
COL_GREY   <- "#D0D0D0"
COL_LIGHT  <- "#F2F2F2"
COL_DARK   <- "#222222"
COL_MID    <- "#555555"

# Fig.2 legacy DNA/RNA palette restored exactly from the mature two-source Fig.2.
DNA_RNA_COLS <- c("DNA" = "#4DBBD5", "RNA" = "#E64B35")

# Evidence concentration uses the same semantic colours as Fig.2 A plus neutral grey.
DIM_COLS <- c(
  "Virus family" = "#4DBBD5",
  "Host family" = "#E64B35",
  "Country" = "#555555"
)

# Legacy Fig.2 continent colours.
CONTINENT_COLS <- c(
  "Asia" = "#E64B35",
  "Europe" = "#4DBBD5",
  "Americas" = "#00A087",
  "Africa" = "#3C5488",
  "Oceania" = "#F39B7F",
  "North America" = "#00A087",
  "South America" = "#00A087"
)

# Legacy Fig.2 geographic blue gradient.
MAP_COLS <- c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")

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
      strip.text = element_text(size = TEXT_SMALL, face = "bold", colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
