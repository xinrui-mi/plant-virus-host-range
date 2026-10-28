required_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "ggplot2", "scales",
  "patchwork", "sf", "rnaturalearth", "ggrepel", "grid"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them into your R library before rerunning."
  )
}
