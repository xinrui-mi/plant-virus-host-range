required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "), "\n",
    "Install them with:\n",
    "Rscript -e 'install.packages(c(",
    paste(sprintf('\\\"%s\\\"', missing_packages), collapse = ","),
    "), repos=\\\"https://cloud.r-project.org\\\")'"
  )
}
