#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

# ---------------------------------------------------------------------------
# Released version.
#
# On the original cluster the input was the third-party workbook
# VIPDB_relation_raw.xlsx, which is not redistributed. The released repository
# ships `data/analysis_data/VIPDB_relation_normalized.tsv.gz`, a normalized
# projection of that same workbook that already carries the four mode classes
# (`fine4_mode`) and the standardized lineage / vector fields this control needs.
# The script therefore accepts EITHER:
#
#   Rscript 04D_transmission_specificity_control.R <VIPDB_relation_raw.xlsx>
#   Rscript 04D_transmission_specificity_control.R <VIPDB_relation_normalized.tsv.gz>
#   Rscript 04D_transmission_specificity_control.R <analysis_root> <input file>
#
# Output goes to <analysis_root>/additional_validation/transmission_mode_consistency
# (previously <analysis_root>/Fig4_transmission_lineage_feasibility_TEST_v2).
# ---------------------------------------------------------------------------
DEFAULT_ROOT <- Sys.getenv("ANALYSIS_ROOT", unset = "")
DEFAULT_TSV <- if (nzchar(DEFAULT_ROOT)) {
  file.path(dirname(DEFAULT_ROOT), "data", "analysis_data", "VIPDB_relation_normalized.tsv.gz")
} else ""

is_tabular_input <- function(p) grepl("\\.(xlsx|xls|tsv|tsv\\.gz|txt|txt\\.gz)$", p, ignore.case = TRUE)

INPUT_FILE <- ""
if (length(args) >= 1 && is_tabular_input(args[1])) {
  ANALYSIS_ROOT <- DEFAULT_ROOT
  INPUT_FILE <- args[1]
} else {
  ANALYSIS_ROOT <- if (length(args) >= 1) args[1] else DEFAULT_ROOT
  INPUT_FILE <- if (length(args) >= 2) args[2] else DEFAULT_TSV
}
if (!nzchar(ANALYSIS_ROOT)) {
  stop("Pass the analysis root as the first argument, or set ANALYSIS_ROOT.\n",
       "Example: Rscript 04D_transmission_specificity_control.R ../.. <input file>")
}
ANALYSIS_ROOT <- normalizePath(ANALYSIS_ROOT, mustWork = TRUE)

outdir <- Sys.getenv("SPECIFICITY_OUTDIR", unset = "")
if (!nzchar(outdir)) {
  # Default: alongside this script, so the released frozen tables are overwritten
  # in place rather than scattered into a new tree.
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  outdir <- if (!is.na(script_path) && nzchar(script_path)) dirname(normalizePath(script_path, mustWork = FALSE))
            else file.path(ANALYSIS_ROOT, "additional_validation", "transmission_mode_consistency")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

present <- function(z) !is.na(z) & nzchar(trimws(as.character(z)))
std_text <- function(z){
  z <- as.character(z); z[is.na(z)] <- ""
  z <- gsub("\u00A0"," ",z,fixed=TRUE)
  trimws(gsub("[[:space:]]+"," ",z))
}
normalize_mode_key <- function(z){
  z <- as.character(z); z[is.na(z)] <- ""
  z <- gsub("\u00A0"," ",z,fixed=TRUE)
  z <- gsub("[‐‑‒–—−]","-",z,perl=TRUE)
  z <- gsub("[[:punct:]]+"," ",z)
  trimws(gsub("[[:space:]]+"," ",tolower(z)))
}

if (!nzchar(INPUT_FILE) || !file.exists(INPUT_FILE)) {
  stop("Provide the VIP-DB input path. Example:\n",
       "Rscript 04D_transmission_specificity_control.R ../.. ",
       "data/analysis_data/VIPDB_relation_normalized.tsv.gz")
}

cat("ANALYSIS_ROOT:", ANALYSIS_ROOT, "\n")
cat("INPUT_FILE:", INPUT_FILE, "\n")
cat("OUTDIR:", outdir, "\n")

# Two accepted input layouts:
#  (a) the original VIPDB workbook  -> raw column names, mode strings need parsing
#  (b) the released normalized TSV   -> already standardised columns
read_vipdb <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    as.data.table(read_excel(path, sheet = "Sheet1", col_types = "text"))
  } else {
    as.data.table(fread(path, sep = "\t", quote = "", na.strings = c("", "NA")))
  }
}
x <- read_vipdb(INPUT_FILE)

# Map the normalized TSV onto the four fields the control consumes.
if (!all(c("Virus Transmission Mode") %in% names(x)) && "fine4_mode" %in% names(x)) {
  setnames(x,
    c("virus_family", "virus_genus", "virus_taxid",
      "vector_genus", "vector_name", "vector_taxid", "raw_transmission_mode"),
    c("Virus Family", "Virus Genus", "Virus TaxID",
      "Vector Genus", "Vector", "Vector TaxID", "Virus Transmission Mode"),
    skip_absent = TRUE)
  # The normalized table already resolved the four-class mode; reuse it but
  # lower-case it so labels match the fcase() branches below exactly (the frozen
  # cluster run also emitted lower-case labels).
  x[, fine4_mode_pre := tolower(std_text(fine4_mode))]
} else {
  x[, fine4_mode_pre := NA_character_]
}

need <- c("Virus Family", "Virus Genus", "Virus TaxID",
          "Vector Genus", "Vector", "Vector TaxID", "Virus Transmission Mode")
miss <- setdiff(need, names(x))
if (length(miss)) stop("VIPDB input missing columns: ", paste(miss, collapse = ", "))

x[, mode_key := normalize_mode_key(`Virus Transmission Mode`)]
x[, fine4_mode := fcase(
  !is.na(fine4_mode_pre) & nzchar(fine4_mode_pre), fine4_mode_pre,
  mode_key == "non persistent transmission", "non-persistent",
  mode_key == "non circulative semi persistent transmission", "semi-persistent",
  mode_key == "circulative persistent non propagative transmission", "persistent non-propagative",
  mode_key == "circulative persistent propagative transmission", "persistent propagative",
  default = NA_character_
)]
x[, virus_taxid := std_text(`Virus TaxID`)]
x[, virus_genus := std_text(`Virus Genus`)]
x[, virus_family := std_text(`Virus Family`)]
x[, vector_taxid := std_text(`Vector TaxID`)]
x[, vector_name := std_text(`Vector`)]
x[, vector_genus := std_text(`Vector Genus`)]
x[, vector_species_id := fifelse(
  present(vector_taxid), paste0("taxid:",vector_taxid),
  fifelse(present(vector_name),paste0("name:",tolower(vector_name)),NA_character_)
)]

d <- unique(x[
  present(virus_taxid) & present(vector_species_id) & present(fine4_mode),
  .(virus_taxid,virus_genus,virus_family,
    vector_species_id,vector_genus,vector_name,vector_taxid,fine4_mode)
])

# One true label per virus × vector species is required.
pairqc <- d[,.(n_modes=uniqueN(fine4_mode),
              modes=paste(sort(unique(fine4_mode)),collapse="|")),
            by=.(virus_taxid,virus_genus,virus_family,
                 vector_species_id,vector_genus,vector_name,vector_taxid)]
fwrite(pairqc[n_modes>1],file.path(outdir,"00_ambiguous_pairs.tsv"),
       sep="\t",na="NA",quote=FALSE)
pairs <- pairqc[n_modes==1]
pairs[,fine4_mode:=modes]
pairs[,c("n_modes","modes"):=NULL]
pairs[,target_id:=.I]

cat("Eligible virus-vector pairs:",nrow(pairs),"\n")
cat("Unique virus TaxIDs:",uniqueN(pairs$virus_taxid),"\n")
cat("Unique vector species:",uniqueN(pairs$vector_species_id),"\n\n")

schemes <- list(
  virus_genus_only = c("virus_genus"),
  virus_family_only = c("virus_family"),
  vector_species_only = c("vector_species_id"),
  vector_genus_only = c("vector_genus"),
  vector_species_x_virus_genus = c("vector_species_id","virus_genus"),
  vector_genus_x_virus_genus = c("vector_genus","virus_genus"),
  vector_species_x_virus_family = c("vector_species_id","virus_family"),
  vector_genus_x_virus_family = c("vector_genus","virus_family")
)

predict_scheme <- function(dat, cols, scheme, min_train=2L, min_purity=1.0){
  out <- vector("list",nrow(dat))
  for(i in seq_len(nrow(dat))){
    target <- dat[i]
    valid_target <- all(vapply(cols,function(cc) present(target[[cc]]),logical(1)))
    z <- data.table(
      target_id=target$target_id,
      virus_taxid=target$virus_taxid,
      true_mode=target$fine4_mode,
      scheme=scheme,
      min_train=min_train,
      min_purity=min_purity,
      train_virus_n=0L,
      train_purity=NA_real_,
      predicted_mode=NA_character_,
      predicted=FALSE,
      correct=FALSE
    )
    if(!valid_target){out[[i]]<-z;next}

    mask <- rep(TRUE,nrow(dat))
    for(cc in cols) mask <- mask & dat[[cc]]==target[[cc]]
    tr <- dat[mask & virus_taxid!=target$virus_taxid,
              .(virus_taxid,fine4_mode)]
    # Resolve each training virus within this grouping; conflicting viruses are dropped.
    trv <- tr[,.(n_modes=uniqueN(fine4_mode),
                 mode=if(uniqueN(fine4_mode)==1) fine4_mode[1] else NA_character_),
              by=virus_taxid]
    trv <- trv[n_modes==1 & present(mode)]
    ntr <- nrow(trv)
    z[,train_virus_n:=ntr]
    if(ntr>=min_train){
      tt <- table(trv$mode); mx <- max(tt); winners <- names(tt)[tt==mx]
      pur <- as.numeric(mx/ntr)
      z[,train_purity:=pur]
      if(length(winners)==1 && pur>=min_purity){
        z[,predicted_mode:=winners[1]]
        z[,predicted:=TRUE]
        z[,correct:=predicted_mode==true_mode]
      }
    }
    out[[i]]<-z
  }
  rbindlist(out)
}

settings <- data.table(
  min_train=c(2L,2L,2L),
  min_purity=c(0.0,0.8,1.0)
)

allpred <- list(); k<-0L
for(nm in names(schemes)){
  for(j in seq_len(nrow(settings))){
    k<-k+1L
    allpred[[k]] <- predict_scheme(
      pairs,schemes[[nm]],nm,
      min_train=settings$min_train[j],
      min_purity=settings$min_purity[j]
    )
  }
}
pred <- rbindlist(allpred)
fwrite(pred,file.path(outdir,"01_predictions_all_schemes.tsv"),sep="\t",na="NA",quote=FALSE)

perf_one <- function(dd){
  pp <- dd[predicted==TRUE]
  recalls <- pp[,.(recall=mean(correct),n=.N),by=true_mode]
  data.table(
    total_pairs=nrow(dd),
    predicted_pairs=nrow(pp),
    unique_predicted_virus_taxid=uniqueN(pp$virus_taxid),
    coverage=if(nrow(dd)) nrow(pp)/nrow(dd) else NA_real_,
    accuracy=if(nrow(pp)) mean(pp$correct) else NA_real_,
    balanced_accuracy=if(nrow(recalls)) mean(recalls$recall) else NA_real_
  )
}
summary <- pred[,perf_one(.SD),by=.(scheme,min_train,min_purity)]
fwrite(summary,file.path(outdir,"02_specificity_control_summary.tsv"),sep="\t",na="NA",quote=FALSE)

recall <- pred[predicted==TRUE,.(
  n=.N,correct_n=sum(correct),recall=mean(correct)
),by=.(scheme,min_train,min_purity,true_mode)]
fwrite(recall,file.path(outdir,"03_recall_by_mode.tsv"),sep="\t",na="NA",quote=FALSE)

# Exact same-target comparison for the strict (>=2 training viruses, unanimous) rule.
strict <- pred[min_train==2 & min_purity==1]
wide <- dcast(strict,target_id+virus_taxid+true_mode ~ scheme,
              value.var="correct",fill=NA)
wide_pred <- dcast(strict,target_id ~ scheme,value.var="predicted",fill=FALSE)

comparisons <- list(
  c("virus_genus_only","vector_species_x_virus_genus"),
  c("virus_genus_only","vector_genus_x_virus_genus"),
  c("virus_family_only","vector_species_x_virus_family"),
  c("virus_family_only","vector_genus_x_virus_family")
)

cmpout <- list()
for(i in seq_along(comparisons)){
  a<-comparisons[[i]][1]; b<-comparisons[[i]][2]
  A <- strict[scheme==a,.(target_id,a_pred=predicted,a_correct=correct)]
  B <- strict[scheme==b,.(target_id,b_pred=predicted,b_correct=correct)]
  M <- merge(A,B,by="target_id")
  both <- M[a_pred & b_pred]
  cmpout[[i]] <- data.table(
    baseline=a,conditioned=b,
    baseline_coverage=mean(M$a_pred),
    conditioned_coverage=mean(M$b_pred),
    common_predicted_n=nrow(both),
    baseline_accuracy_on_common=if(nrow(both)) mean(both$a_correct) else NA_real_,
    conditioned_accuracy_on_common=if(nrow(both)) mean(both$b_correct) else NA_real_
  )
}
cmp <- rbindlist(cmpout)
fwrite(cmp,file.path(outdir,"04_incremental_vector_value.tsv"),sep="\t",na="NA",quote=FALSE)

# Global lineage mode diversity: does family/genus alone already determine the mode?
lineage_genus <- unique(pairs[,.(virus_taxid,virus_genus,fine4_mode)])[present(virus_genus),
  .(virus_n=uniqueN(virus_taxid),
    mode_n=uniqueN(fine4_mode),
    dominant_mode=names(which.max(table(fine4_mode))),
    consistency=max(table(fine4_mode))/.N,
    mode_counts=paste(names(table(fine4_mode)),as.integer(table(fine4_mode)),
                      sep="=",collapse=";")),
  by=virus_genus]
lineage_family <- unique(pairs[,.(virus_taxid,virus_family,fine4_mode)])[present(virus_family),
  .(virus_n=uniqueN(virus_taxid),
    mode_n=uniqueN(fine4_mode),
    dominant_mode=names(which.max(table(fine4_mode))),
    consistency=max(table(fine4_mode))/.N,
    mode_counts=paste(names(table(fine4_mode)),as.integer(table(fine4_mode)),
                      sep="=",collapse=";")),
  by=virus_family]
fwrite(lineage_genus[order(-virus_n)],file.path(outdir,"05_genus_global_mode_consistency.tsv"),sep="\t",na="NA",quote=FALSE)
fwrite(lineage_family[order(-virus_n)],file.path(outdir,"06_family_global_mode_consistency.tsv"),sep="\t",na="NA",quote=FALSE)

cat("============================================================\n")
cat("SPECIFICITY CONTROL: strict rule (min_train=2, purity=1.0)\n")
cat("============================================================\n")
print(summary[min_train==2 & min_purity==1][order(scheme)])

cat("\n============================================================\n")
cat("INCREMENTAL VALUE OF VECTOR CONDITIONING\n")
cat("============================================================\n")
print(cmp)

cat("\nOutputs:\n",outdir,"\n",sep="")
cat("Please send terminal output or:\n")
cat("  02_specificity_control_summary.tsv\n")
cat("  04_incremental_vector_value.tsv\n")
cat("  06_family_global_mode_consistency.tsv\n")
