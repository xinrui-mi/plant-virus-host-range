# DISPLAY-ONLY rules for Fig.1 rebuilt from frozen Fig.2.
# Existing Fig.2 panels are inherited without analytical change.
TOP_VIRUS_FAMILY_N <- 20L
TOP_HOST_FAMILY_N  <- 20L
TOP_COUNTRY_LABEL_N <- 6L

# Geometry in millimetres.
A_HEIGHT_MM <- 24
BC_HEIGHT_MM <- 66
MAP_HEIGHT_MM <- 74
BOTTOM_HEIGHT_MM <- 54

MAIN_WIDTH_MM  <- 180
MAIN_HEIGHT_MM <- A_HEIGHT_MM + BC_HEIGHT_MM + MAP_HEIGHT_MM + BOTTOM_HEIGHT_MM

BC_WIDTHS_MM <- c(90, 90)      # virus-family evidence, host-family evidence
MAP_WIDTH_MM <- 180            # full-width map
BOTTOM_WIDTHS_MM <- c(70, 110) # evidence concentration, rarefaction

DPI_OUT <- 600

# Frozen manuscript totals shown in the new Fig.1A opening strip.
# These are display constants from the finalized three-source analysis.
FIG1_SEQUENCE_RECORDS_N <- 80028L
FIG1_ASSOCIATION_RECORDS_N <- 83034L
FIG1_VC80_N <- 11708L
