if (!exists("R_LIB", inherits = TRUE)) {
  R_LIB <- Sys.getenv("R_LIB", unset = "")
  if (nzchar(R_LIB) && dir.exists(R_LIB)) .libPaths(unique(c(R_LIB, .libPaths())))
}

required_pkgs <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "ape",
  "ggtree"
)

missing <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "Missing R package(s): ", paste(missing, collapse = ", "), "\n\n",
    "This publication plotting script does NOT modify your R environment automatically.\n",
    "Install the missing package(s) first, then rerun.\n",
    "Current R library: ", R_LIB
  )
}

message("Package check passed.")
