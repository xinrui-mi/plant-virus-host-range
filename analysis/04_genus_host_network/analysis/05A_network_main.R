suppressPackageStartupMessages({library(data.table);library(igraph)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"));source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"));source(file.path(ANALYSIS_ROOT,"04_genus_host_network","analysis","05_network_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"04_genus_host_network");make_dirs(base)
# Frozen-output guard (U26, 2026-10-27): `tables/01–06` are the FROZEN results
# the paper's Fig. 4 renders from. A re-run must NOT overwrite them (Louvain is
# a randomised heuristic — a re-run will not land on the published 12-module
# partition). If the frozen tables exist, this run writes to tables/rerun/
# instead; the plotting code only ever reads the frozen tables/.
OUT_TABLES <- file.path(base, "tables")
FROZEN_SENTINEL <- file.path(OUT_TABLES, "03_network_module_membership.tsv")
OUT_QC <- file.path(base, "qc")
if (file.exists(FROZEN_SENTINEL)) {
  OUT_TABLES <- file.path(base, "tables", "rerun")
  OUT_QC <- file.path(base, "qc", "rerun")
  dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_QC, recursive = TRUE, showWarnings = FALSE)
  message("Frozen network tables detected; this re-run writes to ", OUT_TABLES,
          " and ", OUT_QC,
          " and will not overwrite the published partition or its QC file.")
}
assoc<-fread_any(ASSOC_FINAL,sep="\t",showProgress=FALSE);vc<-fread(file.path(ANALYSIS_ROOT,"01_dataset_summary","tables","01_VC80_internal_consistency_QC.tsv"),sep="\t")
assoc<-merge(assoc,vc[,.(vc80_final_id,vc80_support_qc_pass,cross_current_ictv_genus)],by="vc80_final_id",all.x=TRUE)
# Main network uses only current ICTV genus explicitly supported by a current species/genus anchor.
a<-assoc[current_ictv_genus_eligible(virus_ictv_match_rank,virus_genus) & present(host_family) & vc80_support_qc_pass%in%TRUE]
long<-unique(a[,.(virus_node=std_text(virus_genus),host_family=std_text(host_family),vc80_final_id)])
write_tsv(long,file.path(OUT_TABLES,"01_network_input_current_ICTV_genus_host_family_VC80.tsv"))
res<-network_metrics(long,"virus_node",SEED,TRUE)
metric<-res$metrics
# add descriptive family/geography from all associations belonging to eligible current genus
meta<-assoc[current_ictv_genus_eligible(virus_ictv_match_rank,virus_genus),.(virus_family=mode_value(virus_family),host_order_n=uniqueN(host_order[present(host_order)]),country_n=uniqueN(country[present(country)]),continent_n=uniqueN(continent[present(continent)]),sequence_n=uniqueN(sequence_id)),by=.(virus_node=std_text(virus_genus))]
metric<-merge(metric,meta,by="virus_node",all.x=TRUE);setorder(metric,-candidate,-candidate_score)
write_tsv(res$edge_weights,file.path(OUT_TABLES,"02_network_edge_weights.tsv"));write_tsv(res$module_map,file.path(OUT_TABLES,"03_network_module_membership.tsv"));write_tsv(metric,file.path(OUT_TABLES,"04_current_ICTV_genus_network_metrics.tsv"));write_tsv(res$thresholds,file.path(OUT_TABLES,"05_candidate_thresholds.tsv"));write_tsv(metric[candidate==TRUE],file.path(OUT_TABLES,"06_network_candidates.tsv"))
# Cross-genus VC80 is diagnostic only, not a main QC failure. Summarize candidate sensitivity exposure.
cross_diag<-a[,.(cross_genus_VC80_n=uniqueN(vc80_final_id[cross_current_ictv_genus%in%TRUE]),all_VC80_n=uniqueN(vc80_final_id)),by=.(virus_node=std_text(virus_genus))];cross_diag[,cross_genus_fraction:=cross_genus_VC80_n/pmax(all_VC80_n,1)];write_tsv(cross_diag,file.path(OUT_TABLES,"07_cross_genus_VC80_diagnostic_by_genus.tsv"))
qc<-data.table(metric=c("network_genus_n","host_family_n","weighted_edges_n","modules_n","modularity","candidate_n","host_breadth_cutoff","bridge_cutoff","support_cutoff"),value=c(nrow(metric),uniqueN(long$host_family),nrow(res$edge_weights),uniqueN(res$module_map$module_id),res$modularity,sum(metric$candidate),res$thresholds$host_breadth_cutoff,res$thresholds$bridge_cutoff,res$thresholds$support_cutoff));write_tsv(qc,file.path(OUT_QC,"01_Fig5_main_QC.tsv"));cat("Fig5 main network PASS\n");print(qc)
