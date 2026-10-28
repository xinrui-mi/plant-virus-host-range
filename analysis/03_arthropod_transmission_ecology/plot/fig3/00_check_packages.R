required_pkgs <- c(
  "readr", "dplyr", "tidyr", "stringr", "data.table",
  "ggplot2", "patchwork", "scales", "tibble"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}
