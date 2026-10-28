suppressPackageStartupMessages({library(data.table);library(igraph)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"));source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"));source(file.path(ANALYSIS_ROOT,"04_genus_host_network","analysis","05_network_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"04_genus_host_network");make_dirs(base)
long<-fread(file.path(base,"tables","01_network_input_current_ICTV_genus_host_family_VC80.tsv"),sep="\t")
vc<-fread(file.path(ANALYSIS_ROOT,"01_dataset_summary","tables","01_VC80_internal_consistency_QC.tsv"),sep="\t")
full<-fread(file.path(base,"tables","04_current_ICTV_genus_network_metrics.tsv"),sep="\t")
pure_ids<-vc[vc80_support_qc_pass%in%TRUE & !cross_current_ictv_genus,vc80_final_id]
pure<-long[vc80_final_id%in%pure_ids]
r<-network_metrics(pure,"virus_node",SEED+650000L,TRUE);m<-r$metrics
write_tsv(m,file.path(base,"tables","19_genus_pure_VC80_network_metrics_sensitivity.tsv"))
cmp<-merge(full[,.(virus_node,main_host_family_n=host_family_n,main_bridge_score=bridge_score,main_candidate=candidate,main_candidate_score=candidate_score)],m[,.(virus_node,pure_host_family_n=host_family_n,pure_bridge_score=bridge_score,pure_candidate=candidate,pure_candidate_score=candidate_score)],by="virus_node",all=TRUE)
cmp[,candidate_agreement:=main_candidate==pure_candidate]
write_tsv(cmp,file.path(base,"tables","cross_genus_VC80_network_sensitivity_comparison.tsv"))
main_set<-sort(unique(full[candidate%in%TRUE,virus_node]));pure_set<-sort(unique(m[candidate%in%TRUE,virus_node]));shared_set<-intersect(main_set,pure_set);union_set<-union(main_set,pure_set)
change<-data.table(virus_node=union_set);change[,status:=fifelse(virus_node%in%shared_set,"shared_candidate",fifelse(virus_node%in%main_set,"main_only","genus_pure_only"))]
write_tsv(change,file.path(base,"tables","cross_genus_VC80_candidate_set_changes.tsv"))
jaccard<-if(length(union_set))length(shared_set)/length(union_set) else NA_real_;retention<-if(length(main_set))length(shared_set)/length(main_set) else NA_real_;pure_shared<-if(length(pure_set))length(shared_set)/length(pure_set) else NA_real_
qc<-data.table(metric=c("main_VC80_n","genus_pure_VC80_n","main_candidate_n","pure_candidate_n","shared_candidate_n","candidate_union_n","candidate_set_Jaccard","main_candidate_retention_fraction","genus_pure_candidate_shared_fraction","all_node_classification_agreement"),value=c(uniqueN(long$vc80_final_id),uniqueN(pure$vc80_final_id),length(main_set),length(pure_set),length(shared_set),length(union_set),jaccard,retention,pure_shared,mean(cmp$candidate_agreement,na.rm=TRUE)))
write_tsv(qc,file.path(base,"qc","03_Fig5_cross_genus_VC80_sensitivity_QC.tsv"));cat("Fig5 cross-genus VC80 sensitivity PASS\n");print(qc)
