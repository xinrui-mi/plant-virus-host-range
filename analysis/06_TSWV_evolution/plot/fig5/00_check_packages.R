required <- c(
  "data.table", "ggplot2", "scales",
  "ape", "ggtree", "patchwork"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "Missing required R packages: ",
    paste(missing, collapse = ", ")
  )
}

cat("Fig7 publication FINAL v1.3.6 package check PASS\n")
