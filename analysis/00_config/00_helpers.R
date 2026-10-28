suppressPackageStartupMessages({
  library(data.table)
})

std_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  # Use fixed replacement for actual control characters. Do NOT use
  # [\r\n\t] with the default TRE regex engine, because it can treat
  # r/n/t as literal letters and corrupt taxonomy names.
  x <- gsub("\r", " ", x, fixed = TRUE)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- gsub("\t", " ", x, fixed = TRUE)
  x <- trimws(gsub("[[:space:]]+", " ", x))
  x
}

present <- function(x) std_text(x) != ""

# TRUE/FALSE flags stored as text in the released tables (source_NCBI_Virus,
# source_Virus_Host_DB, source_LandFlora).
is_true_flag <- function(x) toupper(std_text(x)) %in% c("TRUE","T","1","YES")

collapse_unique <- function(x) {
  z <- sort(unique(std_text(x)[present(x)]))
  if (length(z) == 0L) "" else paste(z, collapse = ";")
}

mode_value <- function(x) {
  z <- std_text(x)
  z <- z[z != ""]
  if (!length(z)) return("")
  tt <- sort(table(z), decreasing = TRUE)
  cand <- names(tt)[tt == max(tt)]
  sort(cand)[1]
}

scale01 <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(rep(0, length(x)))
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(rep(0, length(x)))
  y <- (x - lo) / (hi - lo)
  y[!is.finite(y)] <- 0
  y
}

safe_quantile <- function(x, p) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, names = FALSE, na.rm = TRUE, type = 7))
}

derive_nucleic_acid <- function(virus_type, genome, rdrp_cluster = NULL) {
  vt <- toupper(std_text(virus_type)); gm <- toupper(std_text(genome))
  txt <- paste(vt, gm)
  out <- rep(NA_character_, length(vt))
  has_rna <- grepl("RNA", txt)
  has_dna <- grepl("DNA", txt)
  out[has_rna & !has_dna] <- "RNA"
  out[has_dna & !has_rna] <- "DNA"
  out[has_rna & has_dna] <- "AMBIGUOUS_DNA_RNA"
  if (!is.null(rdrp_cluster)) {
    rr <- present(rdrp_cluster)
    out[is.na(out) & rr] <- "RNA"
  }
  out
}

current_ictv_family_eligible <- function(match_rank, family) {
  tolower(std_text(match_rank)) %in% c("species", "genus", "family") & present(family)
}

current_ictv_genus_eligible <- function(match_rank, genus) {
  tolower(std_text(match_rank)) %in% c("species", "genus") & present(genus)
}

current_ictv_species_eligible <- function(match_rank, species) {
  tolower(std_text(match_rank)) == "species" & present(species)
}

complete_host_taxonomy <- function(dt) {
  req <- c("host_kingdom","host_phylum","host_class","host_order","host_family","host_genus","host_species","host_taxid")
  Reduce(`&`, lapply(req, function(z) present(dt[[z]])))
}

host_distance_vectors <- function(host, i, j) {
  # host must have one row per host_taxid and complete standard ranks.
  same <- function(col) host[[col]][i] == host[[col]][j]
  d <- rep(NA_integer_, length(i))
  idx <- same("host_kingdom"); d[idx] <- 6L
  idx <- idx & same("host_phylum"); d[idx] <- 5L
  idx <- idx & same("host_class"); d[idx] <- 4L
  idx <- idx & same("host_order"); d[idx] <- 3L
  idx <- idx & same("host_family"); d[idx] <- 2L
  idx <- idx & same("host_genus"); d[idx] <- 1L
  idx <- idx & same("host_species"); d[idx] <- 0L
  d
}

fit_binomial_aggregated <- function(summary_dt, label = "analysis") {
  x <- as.data.frame(summary_dt)
  x <- x[is.finite(x$taxonomy_distance) & x$total_pairs > 0, , drop=FALSE]
  if (nrow(x) < 2 || sum(x$positive_pairs) == 0 || sum(x$positive_pairs) == sum(x$total_pairs)) {
    return(data.table(analysis=label,beta=NA_real_,se=NA_real_,beta_ci_low=NA_real_,beta_ci_high=NA_real_,odds_ratio=NA_real_,or_ci_low=NA_real_,or_ci_high=NA_real_,n_pairs=sum(x$total_pairs),positive_pairs=sum(x$positive_pairs)))
  }
  fit <- suppressWarnings(glm(cbind(positive_pairs, total_pairs-positive_pairs) ~ taxonomy_distance, family=binomial(), data=x))
  cc <- summary(fit)$coefficients
  if (!("taxonomy_distance" %in% rownames(cc))) {
    return(data.table(analysis=label,beta=NA_real_,se=NA_real_,beta_ci_low=NA_real_,beta_ci_high=NA_real_,odds_ratio=NA_real_,or_ci_low=NA_real_,or_ci_high=NA_real_,n_pairs=sum(x$total_pairs),positive_pairs=sum(x$positive_pairs)))
  }
  b <- cc["taxonomy_distance","Estimate"]; se <- cc["taxonomy_distance","Std. Error"]
  data.table(
    analysis=label, beta=b, se=se,
    beta_ci_low=b-1.96*se, beta_ci_high=b+1.96*se,
    odds_ratio=exp(b), or_ci_low=exp(b-1.96*se), or_ci_high=exp(b+1.96*se),
    n_pairs=sum(x$total_pairs), positive_pairs=sum(x$positive_pairs)
  )
}

coef_table <- function(fit, model_name) {
  if (is.null(fit)) return(data.table())
  cc <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
  if (is.null(cc)) return(data.table())
  out <- data.table(term=rownames(cc), estimate=cc[,1], std_error=cc[,2])
  if (ncol(cc) >= 4) { out[, statistic := cc[,3]]; out[, p_value := cc[,4]] }
  out[, model := model_name]
  out[]
}

# Released inputs are gzip-compressed. data.table::fread() cannot decompress them
# without the optional R.utils package, so every read of a released input goes
# through this wrapper. Decompression is transparent: the bytes handed to fread()
# are identical to the uncompressed file, so nothing downstream changes.
fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    if (!file.exists(path)) stop("Released input missing: ", path)
    data.table::fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  } else {
    data.table::fread(path, ...)
  }
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  data.table::fwrite(as.data.table(x), path, sep="\t", na="", quote=FALSE)
}

make_dirs <- function(base) {
  for (z in c("tables","qc","logs","work")) dir.create(file.path(base,z), recursive=TRUE, showWarnings=FALSE)
}

max_host_distance_one_group <- function(h) {
  h <- unique(as.data.table(h), by="host_taxid")
  h <- h[complete_host_taxonomy(h)]
  if (nrow(h) < 2) return(0L)
  ij <- utils::combn(seq_len(nrow(h)), 2L)
  max(host_distance_vectors(h, ij[1,], ij[2,]), na.rm=TRUE)
}
