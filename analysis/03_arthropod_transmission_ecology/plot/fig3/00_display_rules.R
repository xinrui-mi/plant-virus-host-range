# Display-only geometry / typography. No analytical thresholds/results here.
# Typography is EXACTLY inherited from Fig4_publication_FINAL_v25/v26.

MAIN_WIDTH_MM <- 180
A_HEIGHT_MM <- 52
BC_HEIGHT_MM <- 58
D_HEIGHT_MM <- 62
MAIN_HEIGHT_MM <- A_HEIGHT_MM + BC_HEIGHT_MM + D_HEIGHT_MM

# Final main layout: A full width; B (OR) + C (distance decay) on one row; D full width.
A_WIDTH_MM <- 180
BC_WIDTHS_MM <- c(80, 100)
D_WIDTH_MM <- 180

SUPP_WIDTH_MM <- 180
SA_HEIGHT_MM <- 70
SB_HEIGHT_MM <- 94
SUPP_HEIGHT_MM <- SA_HEIGHT_MM + SB_HEIGHT_MM

DPI_OUT <- 600

TEXT_TINY <- 5.7
TEXT_SMALL <- 6.8
TEXT_NORMAL <- 8.2
TEXT_LARGE <- 9.2

LINE_THIN <- 0.28
LINE_MEDIUM <- 0.48

COL_GREY <- "#9A9A9A"
COL_GREY_DARK <- "#6F6F6F"
COL_GREY_LIGHT <- "#C9C9C9"
COL_BLUE <- "#4DBBD5"
COL_BLUE_DARK <- "#2F8FB3"
COL_BLUE_LIGHT <- "#A6DCEF"
COL_ORANGE <- "#E64B35"
COL_ORANGE_LIGHT <- "#F39B7F"

SEQ_COLS <- c(
  "Absent" = COL_GREY,
  "Present" = COL_BLUE_DARK
)

EVIDENCE_COLS <- c(
  "Arthropod-associated sequence evidence" = COL_BLUE,
  "Documented vector relationship" = COL_ORANGE
)

MODE_COLS <- c(
  "Non-persistent" = "#D95F02",
  "Semi-persistent" = "#E6AB02",
  "Persistent non-propagative" = "#7E6148",
  "Persistent propagative" = "#1B9E77"
)

ASSOC_COLS <- c(EVIDENCE_COLS, MODE_COLS)
