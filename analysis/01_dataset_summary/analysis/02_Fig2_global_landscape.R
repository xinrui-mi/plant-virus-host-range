source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_config.R"))
source(file.path(ANALYSIS_ROOT, "00_config", "00_helpers.R"))
# Shipped under 01_dataset_summary: this script produces the global dataset
# description (virus families, host families, countries, RNA evidence space)
# reported in Results 3.1, not the host-relatedness inference of Results 3.2-3.3.
base <- file.path(ANALYSIS_ROOT, "01_dataset_summary"); make_dirs(base)
seq <- fread_any(SEQ_FINAL, sep="\t", showProgress=FALSE)
assoc <- fread_any(ASSOC_FINAL, sep="\t", showProgress=FALSE)
vc <- fread(file.path(ANALYSIS_ROOT,"01_dataset_summary","tables","01_VC80_internal_consistency_QC.tsv"), sep="\t")
seq <- merge(seq, vc[,.(vc80_final_id,vc80_support_qc_pass)], by="vc80_final_id", all.x=TRUE)
seq[, current_family := fifelse(current_ictv_family_eligible(virus_ictv_match_rank,virus_family),std_text(virus_family),"")]
seq[, current_genus := fifelse(current_ictv_genus_eligible(virus_ictv_match_rank,virus_genus),std_text(virus_genus),"")]
seq[, nucleic_acid_class := derive_nucleic_acid(virus_type, genome, rdrp_operational_cluster_id)]
family_seq <- seq[present(current_family), .(
  sequence_n=uniqueN(sequence_id),
  qc_pass_VC80_n=uniqueN(vc80_final_id[vc80_support_qc_pass %in% TRUE]),
  current_ICTV_genus_n=uniqueN(current_genus[present(current_genus)]),
  virus_taxid_n=uniqueN(virus_taxid[present(virus_taxid)])
), by=.(virus_family=current_family)]
family_assoc <- assoc[current_ictv_family_eligible(virus_ictv_match_rank,virus_family),.(association_n=.N,host_family_n=uniqueN(host_family[present(host_family)])),by=virus_family]
family <- merge(family_seq,family_assoc,by="virus_family",all=TRUE); setorder(family,-qc_pass_VC80_n,-sequence_n)
write_tsv(family,file.path(base,"tables","01_virus_family_evidence.tsv"))
hostfam <- assoc[present(host_family),.(
  association_n=.N,
  unique_sequence_n=uniqueN(sequence_id),
  virus_taxid_n=uniqueN(virus_taxid[present(virus_taxid)]),
  current_ICTV_genus_n=uniqueN(virus_genus[current_ictv_genus_eligible(virus_ictv_match_rank,virus_genus)])
),by=host_family][order(-association_n)]
write_tsv(hostfam,file.path(base,"tables","02_host_family_evidence.tsv"))
country <- assoc[present(country),.(association_n=.N,unique_sequence_n=uniqueN(sequence_id),host_family_n=uniqueN(host_family[present(host_family)])),by=.(country,continent)][order(-association_n)]
continent <- assoc[,.(association_n=.N,unique_sequence_n=uniqueN(sequence_id),country_n=uniqueN(country[present(country)])),by=.(continent=fifelse(present(continent),continent,"NA"))][order(-association_n)]
write_tsv(country,file.path(base,"tables","03_country_evidence.tsv")); write_tsv(continent,file.path(base,"tables","04_continent_evidence.tsv"))
rna <- seq[nucleic_acid_class=="RNA"]
rna[, formal_current_genus := present(current_genus)]
rna[, rdrp65 := present(rdrp_operational_cluster_id)]
rna[, coverage_class := fifelse(formal_current_genus & rdrp65,"formal_genus+RdRP65",fifelse(formal_current_genus,"formal_genus_only",fifelse(rdrp65,"RdRP65_only","neither")))]
rna_space <- rna[,.(sequence_n=.N,VC80_n=uniqueN(vc80_final_id)),by=coverage_class]; rna_space[,fraction:=sequence_n/sum(sequence_n)]
write_tsv(rna_space,file.path(base,"tables","05_RNA_classified_novel_evidence_space.tsv"))
# Public-evidence rarefaction over association records; 100 fixed permutations and 101 checkpoints.
set.seed(SEED)
n <- nrow(assoc); checkpoints <- unique(round(seq(1,n,length.out=101)))
rar <- vector("list", FIG2_RAREFACTION_ITER)
for (it in seq_len(FIG2_RAREFACTION_ITER)) {
  ord <- sample.int(n)
  d <- assoc[ord, .(virus_family=std_text(virus_family),host_family=std_text(host_family),country=std_text(country))]
  firstpos <- function(v) {
    ok <- v!=""; z <- data.table(value=v[ok],pos=which(ok)); z[,.(first_pos=min(pos)),by=value]$first_pos
  }
  fv <- firstpos(d$virus_family); fh <- firstpos(d$host_family); fc <- firstpos(d$country)
  rar[[it]] <- data.table(iteration=it,records=checkpoints,
    virus_family_n=vapply(checkpoints,function(k) sum(fv<=k),integer(1)),
    host_family_n=vapply(checkpoints,function(k) sum(fh<=k),integer(1)),
    country_n=vapply(checkpoints,function(k) sum(fc<=k),integer(1)))
}
rar <- rbindlist(rar); write_tsv(rar,file.path(base,"tables","06_evidence_rarefaction_100_permutations.tsv"))
rar_sum <- rar[,.(virus_family_median=median(virus_family_n),virus_family_q025=quantile(virus_family_n,.025),virus_family_q975=quantile(virus_family_n,.975),host_family_median=median(host_family_n),host_family_q025=quantile(host_family_n,.025),host_family_q975=quantile(host_family_n,.975),country_median=median(country_n),country_q025=quantile(country_n,.025),country_q975=quantile(country_n,.975)),by=records]
write_tsv(rar_sum,file.path(base,"tables","07_evidence_rarefaction_summary.tsv"))
tri <- assoc[present(virus_family)&present(host_family)&present(continent),.(association_n=.N,unique_sequence_n=uniqueN(sequence_id)),by=.(virus_family,host_family,continent)]
write_tsv(tri,file.path(base,"tables","08_family_host_family_continent_long.tsv"))

# Source-stratified descriptive landscapes (supplementary/QC; not a prevalence comparison).
src_family <- assoc[current_ictv_family_eligible(virus_ictv_match_rank,virus_family),.(association_n=.N,unique_sequence_n=uniqueN(sequence_id),host_family_n=uniqueN(host_family[present(host_family)])),by=.(source_label_4way,virus_family=std_text(virus_family))]
src_host <- assoc[present(host_family),.(association_n=.N,unique_sequence_n=uniqueN(sequence_id),virus_family_n=uniqueN(virus_family[present(virus_family)])),by=.(source_label_4way,host_family)]
write_tsv(src_family,file.path(base,"tables","09_source_stratified_virus_family_evidence.tsv"))
write_tsv(src_host,file.path(base,"tables","10_source_stratified_host_family_evidence.tsv"))
qc <- data.table(metric=c("virus_family_n","host_family_n","country_n","continent_labels_n","RNA_sequence_n","rarefaction_iterations"),value=c(nrow(family),nrow(hostfam),nrow(country),uniqueN(continent$continent),nrow(rna),FIG2_RAREFACTION_ITER))
write_tsv(qc,file.path(base,"qc","01_Fig2_QC.tsv")); cat("Fig2 analysis PASS\n"); print(qc)
