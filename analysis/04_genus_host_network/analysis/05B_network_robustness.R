suppressPackageStartupMessages({library(data.table);library(igraph);library(parallel)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"));source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"));source(file.path(ANALYSIS_ROOT,"04_genus_host_network","analysis","05_network_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"04_genus_host_network");make_dirs(base)
long<-fread(file.path(base,"tables","01_network_input_current_ICTV_genus_host_family_VC80.tsv"),sep="\t");full<-fread(file.path(base,"tables","04_current_ICTV_genus_network_metrics.tsv"),sep="\t")
clusters<-sort(unique(long$vc80_final_id))
long_path<-file.path(base,"tables","08_network_80pct_VC80_resampling_long.tsv")
required_cols<-c("virus_node","host_family_n","qc_pass_VC80_support","bridge_score","candidate_score","candidate","iteration")
reuse_existing<-FALSE
if (file.exists(long_path) && file.info(long_path)$size>0) {
  old<-tryCatch(fread(long_path,sep="\t"),error=function(e) NULL)
  if (!is.null(old) && all(required_cols%in%names(old))) {
    got<-sort(unique(as.integer(old$iteration)))
    reuse_existing<-identical(got,seq_len(FIG5_ROBUSTNESS_ITER))
    if (reuse_existing) {
      longr<-old
      cat("Fig5 robustness: reusing complete existing resampling table:",long_path,"\n")
    }
  }
}
if (!reuse_existing) {
  set.seed(SEED+500000L)
  worker<-function(it){set.seed(SEED+500000L+it);keep<-sample(clusters,max(2L,floor(length(clusters)*FIG5_ROBUSTNESS_FRACTION)),replace=FALSE);r<-network_metrics(long[vc80_final_id%in%keep],"virus_node",SEED+it,TRUE);m<-r$metrics[,.(virus_node,host_family_n,qc_pass_VC80_support,bridge_score,candidate_score,candidate)];m[,iteration:=it];m}
  rr<-mclapply(seq_len(FIG5_ROBUSTNESS_ITER),worker,mc.cores=min(NCORES,FIG5_ROBUSTNESS_ITER),mc.preschedule=TRUE,mc.cleanup=TRUE)
  longr<-rbindlist(rr,fill=TRUE);write_tsv(longr,long_path)
}
# Explicit numeric coercion avoids data.table group-type instability when median(integer)
# returns integer for some groups and double for others.
rec<-longr[,.(observed_iteration_n=as.integer(uniqueN(iteration)),candidate_recurrence=as.numeric(mean(candidate,na.rm=TRUE)),median_host_family_n=as.numeric(median(host_family_n,na.rm=TRUE)),median_bridge_score=as.numeric(median(bridge_score,na.rm=TRUE)),median_candidate_score=as.numeric(median(candidate_score,na.rm=TRUE))),by=virus_node]
rec<-merge(full[,.(virus_node,full_candidate=candidate,full_candidate_score=candidate_score,full_host_family_n=host_family_n,full_bridge_score=bridge_score)],rec,by="virus_node",all.x=TRUE);rec[is.na(candidate_recurrence),candidate_recurrence:=0];setorder(rec,-full_candidate,-candidate_recurrence,-full_candidate_score);write_tsv(rec,file.path(base,"tables","09_candidate_recurrence.tsv"))
# rank/metric stability per iteration
stab<-longr[,{z<-merge(.SD,full[,.(virus_node,full_candidate_score=candidate_score,full_bridge=bridge_score,full_breadth=host_family_n)],by="virus_node");list(candidate_score_spearman=as.numeric(cor(z$candidate_score,z$full_candidate_score,method="spearman",use="complete.obs")),bridge_spearman=as.numeric(cor(z$bridge_score,z$full_bridge,method="spearman",use="complete.obs")),breadth_spearman=as.numeric(cor(z$host_family_n,z$full_breadth,method="spearman",use="complete.obs")))},by=iteration];write_tsv(stab,file.path(base,"tables","10_network_metric_stability.tsv"))
# threshold perturbation: predefined quantile/support grid, not tuned after seeing candidates.
grid<-CJ(host_q=c(.90,.95,.975),bridge_q=c(.50,.75,.90),support_cutoff=c(5L,10L,20L));out<-list();k<-0L
for(i in seq_len(nrow(grid))){g<-grid[i];hc<-ceiling(safe_quantile(full$host_family_n,g$host_q)-1e-12);bnz<-full$bridge_score[is.finite(full$bridge_score)&full$bridge_score>0];bc<-if(length(bnz))safe_quantile(bnz,g$bridge_q) else 0;z<-full[,.(virus_node)];z[,candidate:=full$host_family_n>=hc & full$bridge_score>=bc & full$qc_pass_VC80_support>=g$support_cutoff];z[,`:=`(host_q=g$host_q,bridge_q=g$bridge_q,support_cutoff=g$support_cutoff,host_cutoff=hc,bridge_cutoff=bc)];k<-k+1L;out[[k]]<-z}
th<-rbindlist(out);write_tsv(th,file.path(base,"tables","11_candidate_threshold_perturbation_long.tsv"));ths<-th[,.(candidate_scenario_rate=as.numeric(mean(candidate))),by=virus_node];write_tsv(ths,file.path(base,"tables","12_candidate_threshold_perturbation_stability.tsv"));fc<-rec[full_candidate==TRUE,candidate_recurrence];qc<-data.table(metric=c("iterations","resample_fraction","full_candidate_n","full_candidate_median_recurrence","full_candidate_min_recurrence"),value=as.numeric(c(FIG5_ROBUSTNESS_ITER,FIG5_ROBUSTNESS_FRACTION,sum(full$candidate),if(length(fc))median(fc) else NA_real_,if(length(fc))min(fc) else NA_real_)));write_tsv(qc,file.path(base,"qc","02_Fig5_robustness_QC.tsv"));cat("Fig5 robustness PASS\n");print(qc)
