suppressPackageStartupMessages({library(data.table);library(Matrix);library(MASS)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"))
source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base <- file.path(ANALYSIS_ROOT,"03_arthropod_transmission_ecology"); make_dirs(base)
seq <- fread_any(SEQ_FINAL,sep="\t",showProgress=FALSE)
assoc <- fread_any(ASSOC_FINAL,sep="\t",showProgress=FALSE)
vc <- fread(file.path(ANALYSIS_ROOT,"01_dataset_summary","tables","01_VC80_internal_consistency_QC.tsv"),sep="\t")
vip <- fread(file.path(base,"tables","02_VIPDB_virus_level_summary.tsv"),sep="\t")
seq <- merge(seq,vc[,.(vc80_final_id,vc80_support_qc_pass)],by="vc80_final_id",all.x=TRUE)
seq[, arth_supported := toupper(std_text(arthropod_sequence_supported))=="TRUE"]
seq[, current_family := fifelse(current_ictv_family_eligible(virus_ictv_match_rank,virus_family),std_text(virus_family),"")]
# sequence-level evidence aggregated to formal NCBI Virus TaxID
sm <- seq[present(virus_taxid),.(
  sequence_n=uniqueN(sequence_id),
  arthropod_supported_sequence_n=uniqueN(sequence_id[arth_supported]),
  arthropod_support_fraction=mean(arth_supported),
  any_arthropod_sequence_support=any(arth_supported),
  qc_pass_VC80_support=uniqueN(vc80_final_id[vc80_support_qc_pass %in% TRUE]),
  virus_family=mode_value(current_family),
  virus_genus=mode_value(virus_genus)
),by=virus_taxid]
# Association breadth per formal TaxID.
am <- assoc[present(virus_taxid),.(
  association_n=.N,
  host_taxid_n=uniqueN(host_taxid[present(host_taxid)]),
  host_family_n=uniqueN(host_family[present(host_family)]),
  host_order_n=uniqueN(host_order[present(host_order)]),
  country_n=uniqueN(country[present(country)]),
  continent_n=uniqueN(continent[present(continent)])
),by=virus_taxid]
metric <- merge(sm,am,by="virus_taxid",all.x=TRUE)
# Maximum host taxonomic distance. Only current complete host lineages can contribute.
complete <- assoc[present(virus_taxid) & complete_host_taxonomy(assoc)]
span <- complete[, .(max_host_taxonomic_distance=max_host_distance_one_group(.SD)), by=virus_taxid,
  .SDcols=c("host_taxid","host_species","host_genus","host_family","host_order","host_class","host_phylum","host_kingdom")]
metric <- merge(metric,span,by="virus_taxid",all.x=TRUE); metric[is.na(max_host_taxonomic_distance),max_host_taxonomic_distance:=0]
metric <- merge(metric,vip,by="virus_taxid",all.x=TRUE)
metric[, documented_vector := !is.na(relation_row_n) & relation_row_n>0]
metric[is.na(documented_vector),documented_vector:=FALSE]
metric[, seq_support_group:=ifelse(any_arthropod_sequence_support,"sequence-supported","no sequence support")]
write_tsv(metric,file.path(base,"tables","04_formal_taxid_arthropod_vector_metrics.tsv"))
# Descriptive summaries.
desc_seq <- metric[,.(virus_taxid_n=.N,host_family_median=as.numeric(median(host_family_n,na.rm=TRUE)),host_family_q25=as.numeric(quantile(host_family_n,.25,na.rm=TRUE,names=FALSE)),host_family_q75=as.numeric(quantile(host_family_n,.75,na.rm=TRUE,names=FALSE)),host_order_median=as.numeric(median(host_order_n,na.rm=TRUE)),max_distance_median=as.numeric(median(max_host_taxonomic_distance,na.rm=TRUE)),VC80_support_median=as.numeric(median(qc_pass_VC80_support,na.rm=TRUE))),by=seq_support_group]
desc_vec <- metric[,.(virus_taxid_n=.N,host_family_median=as.numeric(median(host_family_n,na.rm=TRUE)),host_order_median=as.numeric(median(host_order_n,na.rm=TRUE)),max_distance_median=as.numeric(median(max_host_taxonomic_distance,na.rm=TRUE))),by=documented_vector]
write_tsv(desc_seq,file.path(base,"tables","05_sequence_evidence_breadth_summary.tsv"));write_tsv(desc_vec,file.path(base,"tables","06_documented_vector_breadth_summary.tsv"))
# Regression helpers. Family-adjusted models are sensitivity models, not causal claims.
fit_count <- function(outcome, exposure, dat, tag) {
  d <- copy(dat)[is.finite(get(outcome)) & get(outcome)>=0 & present(virus_family)]
  if (nrow(d)<30 || uniqueN(d[[exposure]])<2) return(data.table())
  form <- as.formula(paste(outcome,"~",exposure,"+ log1p(qc_pass_VC80_support) + factor(virus_family)"))
  # Negative-binomial theta estimation was unstable; use the
  # pre-specified quasi-Poisson model for overdispersed count outcomes.
  qp <- tryCatch(glm(form,family=quasipoisson(link="log"),data=d,control=glm.control(maxit=200)),error=function(e)NULL)
  coef_table(qp,paste0(tag,"_quasipoisson"))
}
model_terms <- rbindlist(list(
 fit_count("host_family_n","seq_support_group",metric,"host_family_by_sequence_evidence"),
 fit_count("host_order_n","seq_support_group",metric,"host_order_by_sequence_evidence"),
 fit_count("host_family_n","documented_vector",metric,"host_family_by_documented_vector"),
 fit_count("host_order_n","documented_vector",metric,"host_order_by_documented_vector")
),fill=TRUE)
# Continuous max-distance sensitivity.
for (ex in c("seq_support_group","documented_vector")) {
 d<-metric[present(virus_family) & is.finite(max_host_taxonomic_distance)]
 if (uniqueN(d[[ex]])>=2) {
  f<-lm(as.formula(paste("max_host_taxonomic_distance ~",ex,"+ log1p(qc_pass_VC80_support) + factor(virus_family)")),data=d)
  model_terms<-rbind(model_terms,coef_table(f,paste0("max_distance_by_",ex)),fill=TRUE)
 }
}
write_tsv(model_terms,file.path(base,"tables","07_support_and_family_adjusted_models.tsv"))
# Fine four-mode comparison: only virus TaxIDs with exactly one documented four-class mode.
mode4 <- metric[fine4_status=="single_four_class_mode" & present(fine4_mode)]
mode4_desc <- mode4[,.(virus_taxid_n=.N,host_family_median=as.numeric(median(host_family_n,na.rm=TRUE)),host_order_median=as.numeric(median(host_order_n,na.rm=TRUE)),max_distance_median=as.numeric(median(max_host_taxonomic_distance,na.rm=TRUE)),VC80_support_median=as.numeric(median(qc_pass_VC80_support,na.rm=TRUE))),by=fine4_mode]
write_tsv(mode4_desc,file.path(base,"tables","08_transmission_four_mode_summary.tsv"))
mode_tests <- list(); k<-0L
for(outcome in c("host_family_n","host_order_n","max_host_taxonomic_distance")) {
 if(nrow(mode4)>10 && uniqueN(mode4$fine4_mode)>1){
  kr<-kruskal.test(mode4[[outcome]]~factor(mode4$fine4_mode));k<-k+1L;mode_tests[[k]]<-data.table(test="Kruskal-Wallis",outcome=outcome,statistic=unname(kr$statistic),df=unname(kr$parameter),p_value=kr$p.value)
  lmfit<-lm(as.formula(paste0("log1p(",outcome,") ~ factor(fine4_mode) + log1p(qc_pass_VC80_support) + factor(virus_family)")),data=mode4[present(virus_family)])
  model_terms<-rbind(model_terms,coef_table(lmfit,paste0(outcome,"_four_mode_loglinear")),fill=TRUE)
 }
}
write_tsv(rbindlist(mode_tests,fill=TRUE),file.path(base,"tables","09_transmission_four_mode_global_tests.tsv"))
# Pairwise non-parametric contrasts are descriptive and BH-adjusted within each outcome.
pw <- list(); pk <- 0L
for(outcome in c("host_family_n","host_order_n","max_host_taxonomic_distance")){
  lev <- sort(unique(mode4$fine4_mode))
  if(length(lev)>1){
    cmb <- combn(lev,2,simplify=FALSE)
    for(z in cmb){
      d <- mode4[fine4_mode%in%z]
      wt <- tryCatch(wilcox.test(d[[outcome]] ~ factor(d$fine4_mode), exact=FALSE), error=function(e)NULL)
      if(!is.null(wt)){pk<-pk+1L;pw[[pk]]<-data.table(outcome=outcome,group1=z[1],group2=z[2],W=unname(wt$statistic),p_value=wt$p.value)}
    }
  }
}
pw <- rbindlist(pw,fill=TRUE);if(nrow(pw))pw[,p_BH:=p.adjust(p_value,"BH"),by=outcome]
write_tsv(pw,file.path(base,"tables","10_transmission_four_mode_pairwise_Wilcoxon.tsv"))
# Count-model sensitivity for four documented transmission classes.
if(nrow(mode4)>30 && uniqueN(mode4$fine4_mode)>1){
  model_terms<-rbind(model_terms,fit_count("host_family_n","fine4_mode",mode4,"host_family_by_four_mode"),fill=TRUE)
  model_terms<-rbind(model_terms,fit_count("host_order_n","fine4_mode",mode4,"host_order_by_four_mode"),fill=TRUE)
}
# Two retention groups: persistent circulative / propagation unspecified can enter here if its virus-level group is unambiguous.
ret <- metric[retention_status=="single_group" & present(retention_group)]
ret_desc <- ret[,.(virus_taxid_n=.N,host_family_median=as.numeric(median(host_family_n,na.rm=TRUE)),host_order_median=as.numeric(median(host_order_n,na.rm=TRUE)),max_distance_median=as.numeric(median(max_host_taxonomic_distance,na.rm=TRUE))),by=retention_group]
write_tsv(ret_desc,file.path(base,"tables","11_retention_group_summary.tsv"))
ret_tests <- list();rk<-0L
if(uniqueN(ret$retention_group)==2){
 for(outcome in c("host_family_n","host_order_n","max_host_taxonomic_distance")){
  wt<-wilcox.test(ret[[outcome]]~factor(ret$retention_group),exact=FALSE);rk<-rk+1L;ret_tests[[rk]]<-data.table(outcome=outcome,W=unname(wt$statistic),p_value=wt$p.value)
 }
 model_terms<-rbind(model_terms,fit_count("host_family_n","retention_group",ret,"host_family_by_retention_group"),fill=TRUE)
 model_terms<-rbind(model_terms,fit_count("host_order_n","retention_group",ret,"host_order_by_retention_group"),fill=TRUE)
}
write_tsv(rbindlist(ret_tests,fill=TRUE),file.path(base,"tables","12_retention_group_tests.tsv"))
write_tsv(model_terms,file.path(base,"tables","07_support_and_family_adjusted_models.tsv"))
# Host-distance sharing slopes for sequence-supported vs unsupported formal TaxIDs.
# This is a descriptive interaction; inference is interpreted with the Fig3 null/sensitivity framework.
hostcols<-c("host_taxid","host_species","host_genus","host_family","host_order","host_class","host_phylum","host_kingdom")
host<-unique(assoc[complete_host_taxonomy(assoc),..hostcols]);if(host[,uniqueN(host_taxid)]!=nrow(host))stop("host taxonomy conflict")
host<-host[order(host_taxid)];host[,hi:=.I];hm<-host[,.(host_taxid,hi)];ij<-combn(seq_len(nrow(host)),2L);pi<-ij[1,];pj<-ij[2,];dd<-host_distance_vectors(host,pi,pj)
formal_groups<-metric[,.(virus_taxid,seq_support_group)]
e<-unique(merge(assoc[present(virus_taxid)&complete_host_taxonomy(assoc),.(host_taxid,virus_taxid)],formal_groups,by="virus_taxid"));e<-merge(e,hm,by="host_taxid")
interaction_rows<-list();z<-0L
for(gp in c("no sequence support","sequence-supported")){
 ee<-e[seq_support_group==gp];hh<-sort(unique(ee$hi));mask<-pi%in%hh&pj%in%hh;vl<-sort(unique(ee$virus_taxid));X<-sparseMatrix(i=ee$hi,j=match(ee$virus_taxid,vl),x=1,dims=c(nrow(host),length(vl)));SS<-tcrossprod(X);sh<-as.integer(as.numeric(SS[cbind(pi[mask],pj[mask])])>0);z<-z+1L;interaction_rows[[z]]<-data.table(group=gp,taxonomy_distance=dd[mask],share=sh)[,.(total_pairs=.N,positive_pairs=sum(share),sharing_probability=mean(share)),by=.(group,taxonomy_distance)]
}
inter<-rbindlist(interaction_rows);write_tsv(inter,file.path(base,"tables","13_host_distance_by_sequence_evidence_group.tsv"))
inter_fit<-glm(cbind(positive_pairs,total_pairs-positive_pairs)~taxonomy_distance*factor(group),family=binomial(),data=inter);write_tsv(coef_table(inter_fit,"descriptive_host_distance_x_sequence_evidence"),file.path(base,"tables","14_host_distance_interaction_model.tsv"))
# RdRP65 extension, with main A+B LF evidence and strict A-only LF sensitivity.
seq_r <- seq[present(rdrp_operational_cluster_id), .(any_seq_support=any(arth_supported)), by=rdrp_operational_cluster_id]
calc_rdrp <- function(strict=FALSE){
 a<-copy(assoc)[present(rdrp_operational_cluster_id)]
 if(strict){islf<-is_true_flag(a$source_LandFlora)|present(a$landflora_record_id);okstrict<-toupper(std_text(a$landflora_strict_plant_infection))%in%c("TRUE","1","YES");a<-a[!islf|okstrict]}
 m<-a[,.(host_family_n=uniqueN(host_family[present(host_family)]),host_order_n=uniqueN(host_order[present(host_order)])),by=rdrp_operational_cluster_id];merge(m,seq_r,by="rdrp_operational_cluster_id",all.x=TRUE)[,dataset:=ifelse(strict,"strict_LF_A_only","main_LF_AplusB")]
}
rr<-rbind(calc_rdrp(FALSE),calc_rdrp(TRUE));write_tsv(rr,file.path(base,"tables","15_RdRP65_arthropod_sequence_evidence_sensitivity.tsv"));rrsum<-rr[,.(lineage_n=.N,host_family_median=as.numeric(median(host_family_n,na.rm=TRUE)),host_order_median=as.numeric(median(host_order_n,na.rm=TRUE))),by=.(dataset,any_seq_support)];write_tsv(rrsum,file.path(base,"tables","16_RdRP65_sensitivity_summary.tsv"))
qc<-data.table(metric=c("formal_virus_taxid_n","sequence_supported_taxid_n","VIPDB_documented_vector_taxid_n_in_dataset","four_mode_single_taxid_n","retention_single_group_taxid_n","four_mode_multi_taxid_n","RdRP65_lineage_rows_main_and_strict"),value=c(nrow(metric),sum(metric$any_arthropod_sequence_support),sum(metric$documented_vector),nrow(mode4),nrow(ret),sum(metric$fine4_status=="multiple_four_class_modes",na.rm=TRUE),nrow(rr)))
write_tsv(qc,file.path(base,"qc","01_Fig4_QC.tsv"));cat("Fig4 analysis PASS\n");print(qc)
