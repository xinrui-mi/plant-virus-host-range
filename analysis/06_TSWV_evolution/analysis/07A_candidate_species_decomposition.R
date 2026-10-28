suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"))
source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base <- file.path(ANALYSIS_ROOT,"06_TSWV_evolution"); make_dirs(base)
seq <- fread_any(SEQ_FINAL, sep="\t", showProgress=FALSE)
assoc <- fread_any(ASSOC_FINAL, sep="\t", showProgress=FALSE)
vc <- fread(file.path(ANALYSIS_ROOT,"01_dataset_summary","tables","01_VC80_internal_consistency_QC.tsv"),sep="\t")
cand <- fread(file.path(ANALYSIS_ROOT,"04_genus_host_network","tables","06_network_candidates.tsv"),sep="\t")
vip <- fread(file.path(ANALYSIS_ROOT,"03_arthropod_transmission_ecology","tables","02_VIPDB_virus_level_summary.tsv"),sep="\t")
# TaxIDs are identifiers; normalize join-key type explicitly.
seq[,virus_taxid:=std_text(virus_taxid)];assoc[,virus_taxid:=std_text(virus_taxid)];vip[,virus_taxid:=std_text(virus_taxid)]
seq <- merge(seq, vc[,.(vc80_final_id,vc80_support_qc_pass)], by="vc80_final_id", all.x=TRUE)
cg <- sort(unique(std_text(cand$virus_node)[present(cand$virus_node)]))
if (!length(cg)) stop("No Fig5 candidate genera available")
# Current formal species only: exact current ICTV species anchor, not source-retained species under higher anchor.
sq <- seq[tolower(std_text(virus_ictv_match_rank))=="species" & present(virus_species) & present(virus_genus) & std_text(virus_genus)%in%cg]
as <- assoc[tolower(std_text(virus_ictv_match_rank))=="species" & present(virus_species) & present(virus_genus) & std_text(virus_genus)%in%cg]
if (!nrow(sq) || !nrow(as)) stop("No current formal species found within Fig5 candidate genera")
# Parent-genus breadth is defined from all current-genus eligible associations, independent of species availability.
genus_breadth <- assoc[current_ictv_genus_eligible(virus_ictv_match_rank,virus_genus) & std_text(virus_genus)%in%cg,
  .(parent_genus_host_family_n=uniqueN(host_family[present(host_family)]),
    parent_genus_host_order_n=uniqueN(host_order[present(host_order)])),
  by=.(virus_genus=std_text(virus_genus))]
sm <- sq[,.(sequence_n=uniqueN(sequence_id),
            virus_taxid_n=uniqueN(virus_taxid[present(virus_taxid)]),
            virus_taxid_set=collapse_unique(virus_taxid),
            qc_pass_VC80_support=uniqueN(vc80_final_id[vc80_support_qc_pass%in%TRUE]),
            total_VC80_support=uniqueN(vc80_final_id),
            arthropod_supported_sequence_n=uniqueN(sequence_id[toupper(std_text(arthropod_sequence_supported))=="TRUE"])),
         by=.(virus_genus=std_text(virus_genus),virus_species=std_text(virus_species))]
am <- as[,.(association_n=.N,
            host_taxid_n=uniqueN(host_taxid[present(host_taxid)]),
            host_family_n=uniqueN(host_family[present(host_family)]),
            host_order_n=uniqueN(host_order[present(host_order)]),
            country_n=uniqueN(country[present(country)]),
            continent_n=uniqueN(continent[present(continent)])),
         by=.(virus_genus=std_text(virus_genus),virus_species=std_text(virus_species))]
out <- merge(sm,am,by=c("virus_genus","virus_species"),all=TRUE)
out <- merge(out,genus_breadth,by="virus_genus",all.x=TRUE)
out[,host_family_contribution_to_parent_genus := host_family_n/pmax(parent_genus_host_family_n,1)]
out[,host_order_contribution_to_parent_genus := host_order_n/pmax(parent_genus_host_order_n,1)]
# Exact TaxID VIPDB context summarized across all NCBI TaxIDs carried by the current species.
sp_tax <- unique(sq[present(virus_taxid),.(virus_genus=std_text(virus_genus),virus_species=std_text(virus_species),virus_taxid=std_text(virus_taxid))])
spvip <- merge(sp_tax,vip,by="virus_taxid",all.x=TRUE)
spvip <- spvip[,.(VIPDB_documented_vector=any(!is.na(relation_row_n)&relation_row_n>0),
                    VIPDB_matched_taxid_n=sum(!is.na(relation_row_n)&relation_row_n>0),
                    VIPDB_fine_mode_set=collapse_unique(fine_mode_set),
                    VIPDB_four_class_mode_set=collapse_unique(fine4_mode),
                    VIPDB_retention_group_set=collapse_unique(retention_group)),
                 by=.(virus_genus,virus_species)]
out <- merge(out,spvip,by=c("virus_genus","virus_species"),all.x=TRUE)
out[is.na(VIPDB_documented_vector),VIPDB_documented_vector:=FALSE]
# Parent-genus network metrics for context, not a new species risk score.
cm <- cand[,.(virus_genus=virus_node,parent_genus_candidate_score=candidate_score,parent_genus_bridge_score=bridge_score,parent_genus_QCpass_VC80_support=qc_pass_VC80_support,parent_genus_host_family_network_n=host_family_n)]
out <- merge(out,cm,by="virus_genus",all.x=TRUE)
setorder(out,-host_family_contribution_to_parent_genus,-host_family_n,-qc_pass_VC80_support)
write_tsv(out,file.path(base,"tables","01_candidate_genus_all_current_species_metrics.tsv"))
# Historical six-pathogen preselection is deliberately not reused in the new analysis.
qc <- data.table(
  metric=c("Fig5_candidate_genus_n","current_species_n_within_candidates","candidate_species_with_VIPDB_vector_n"),
  value=c(length(cg),nrow(out),sum(out$VIPDB_documented_vector))
)
write_tsv(qc,file.path(base,"qc","01_Fig7_species_decomposition_QC.tsv"))
cat("Fig7 candidate-species decomposition PASS\n")
print(qc)
