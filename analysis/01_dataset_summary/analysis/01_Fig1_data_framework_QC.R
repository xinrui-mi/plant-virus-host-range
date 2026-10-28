source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_config.R"))
source(file.path(ANALYSIS_ROOT, "00_config", "00_helpers.R"))
base <- file.path(ANALYSIS_ROOT, "01_dataset_summary"); make_dirs(base)
seq <- fread_any(SEQ_FINAL, sep="\t", showProgress=FALSE)
assoc <- fread_any(ASSOC_FINAL, sep="\t", showProgress=FALSE)
seq[, nucleic_acid_class := derive_nucleic_acid(virus_type, genome, rdrp_operational_cluster_id)]
seq[, current_ictv_family := fifelse(current_ictv_family_eligible(virus_ictv_match_rank, virus_family), std_text(virus_family), "")]
seq[, current_ictv_genus := fifelse(current_ictv_genus_eligible(virus_ictv_match_rank, virus_genus), std_text(virus_genus), "")]
if (any(!present(seq$vc80_final_id))) stop("Blank vc80_final_id in the released sequence input")
vc <- seq[, .(
  member_sequence_n=.N,
  known_nucleic_acid_n=sum(present(nucleic_acid_class)),
  nucleic_acid_n=uniqueN(nucleic_acid_class[present(nucleic_acid_class)]),
  nucleic_acid_set=collapse_unique(nucleic_acid_class),
  DNA_member_n=sum(nucleic_acid_class=="DNA",na.rm=TRUE),
  RNA_member_n=sum(nucleic_acid_class=="RNA",na.rm=TRUE),
  known_family_n=sum(present(virus_family)),
  virus_family_n=uniqueN(std_text(virus_family)[present(virus_family)]),
  virus_family_set=collapse_unique(virus_family),
  current_ictv_genus_n=uniqueN(current_ictv_genus[present(current_ictv_genus)]),
  current_ictv_genus_set=collapse_unique(current_ictv_genus),
  ambiguous_nucleic_acid_member_n=sum(nucleic_acid_class=="AMBIGUOUS_DNA_RNA",na.rm=TRUE)
), by=vc80_final_id]
vc[, `:=`(
  cross_nucleic_acid = DNA_member_n>0L & RNA_member_n>0L,
  cross_virus_family = virus_family_n > 1L,
  cross_current_ictv_genus = current_ictv_genus_n > 1L
)]
vc[, vc80_support_qc_pass := !cross_nucleic_acid & !cross_virus_family & ambiguous_nucleic_acid_member_n==0L]
setorder(vc, vc80_final_id)
write_tsv(vc, file.path(base,"tables","01_VC80_internal_consistency_QC.tsv"))
map <- vc[, .(vc80_final_id, vc80_support_qc_pass, cross_nucleic_acid, cross_virus_family, cross_current_ictv_genus, ambiguous_nucleic_acid_member_n)]
map <- merge(seq[,.(sequence_id,vc80_final_id)], map, by="vc80_final_id", all.x=TRUE)
write_tsv(map, file.path(base,"tables","02_sequence_VC80_QC_map.tsv"))
# Source provenance uses the released final provenance coding only (source_label_4way).
source_comp <- seq[,.(sequence_n=.N), by=source_label_4way][order(-sequence_n)]
assoc_comp <- assoc[,.(association_n=.N,unique_sequence_n=uniqueN(sequence_id)), by=.(source_label_4way)][order(-association_n)]
# Sequence-level provenance classes are inferred only when all association rows for a sequence agree on source label.
seqprov <- assoc[,.(source_label_set=collapse_unique(source_label_4way),source_label_n=uniqueN(std_text(source_label_4way)[present(source_label_4way)])),by=sequence_id]
seqprov[,source_provenance_class:=fifelse(source_label_n==1 & present(source_label_set),source_label_set,"multiple/ambiguous")]
seqprov_sum<-seqprov[,.(sequence_n=.N),by=source_provenance_class][order(-sequence_n)]
write_tsv(source_comp, file.path(base,"tables","03_sequence_source_composition.tsv"))
write_tsv(assoc_comp, file.path(base,"tables","04_association_source_provenance_composition.tsv"))
write_tsv(seqprov_sum, file.path(base,"tables","04B_sequence_source_provenance_composition.tsv"))
fields_seq <- c("virus_realm","virus_kingdom","virus_phylum","virus_class","virus_order","virus_family","virus_genus","virus_species","virus_taxid","vc80_final_id","rdrp_operational_cluster_id","arthropod_sequence_supported")
cov1 <- rbindlist(lapply(fields_seq, function(z) data.table(entity="sequence",field=z,present_n=sum(present(seq[[z]])),total_n=nrow(seq))))
fields_assoc <- c("host_kingdom","host_phylum","host_class","host_order","host_family","host_genus","host_species","host_taxid","country","continent")
cov2 <- rbindlist(lapply(fields_assoc, function(z) data.table(entity="association",field=z,present_n=sum(present(assoc[[z]])),total_n=nrow(assoc))))
cov <- rbind(cov1,cov2); cov[, fraction := present_n/total_n]
write_tsv(cov, file.path(base,"tables","05_taxonomy_metadata_coverage.tsv"))
seq[, formal_taxid := present(virus_taxid)]
seq[, formal_current_genus := current_ictv_genus_eligible(virus_ictv_match_rank, virus_genus)]
seq[, rdrp65 := present(rdrp_operational_cluster_id)]
units <- data.table(
  unit=c("sequence","VC80_all","VC80_support_QC_pass","NCBI_Virus_TaxID","ICTV_current_family","ICTV_current_genus","RdRP65_operational_lineage"),
  n=c(nrow(seq),uniqueN(seq$vc80_final_id),sum(vc$vc80_support_qc_pass),uniqueN(seq$virus_taxid[present(seq$virus_taxid)]),uniqueN(seq$current_ictv_family[present(seq$current_ictv_family)]),uniqueN(seq$current_ictv_genus[seq$formal_current_genus]),uniqueN(seq$rdrp_operational_cluster_id[seq$rdrp65]))
)
write_tsv(units, file.path(base,"tables","06_analysis_unit_coverage.tsv"))
rna <- seq[nucleic_acid_class=="RNA"]
rna[, coverage_class := fifelse(formal_current_genus & rdrp65, "formal_genus+RdRP65", fifelse(formal_current_genus, "formal_genus_only", fifelse(rdrp65, "RdRP65_only", "neither")))]
rna_cov <- rna[,.(sequence_n=.N),by=coverage_class]; rna_cov[,fraction:=sequence_n/sum(sequence_n)]
write_tsv(rna_cov, file.path(base,"tables","07_RNA_formal_genus_RdRP65_coverage.tsv"))
qc <- data.table(
  metric=c("sequence_rows","association_rows","VC80_clusters","VC80_main_QC_pass","VC80_cross_nucleic_acid","VC80_ambiguous_nucleic_acid_annotation","VC80_cross_family","VC80_cross_current_genus_diagnostic","RNA_sequences","current_ICTV_genus_sequences","RdRP65_sequences"),
  value=c(nrow(seq),nrow(assoc),nrow(vc),sum(vc$vc80_support_qc_pass),sum(vc$cross_nucleic_acid),sum(vc$ambiguous_nucleic_acid_member_n>0),sum(vc$cross_virus_family),sum(vc$cross_current_ictv_genus),nrow(rna),sum(seq$formal_current_genus),sum(seq$rdrp65))
)
write_tsv(qc, file.path(base,"qc","01_Fig1_QC.tsv"))
cat("Fig1 analysis PASS\n"); print(qc)
