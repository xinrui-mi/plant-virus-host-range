suppressPackageStartupMessages({library(data.table);library(Matrix)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"));source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"02_host_relatedness_and_phylogeny");make_dirs(base);o<-readRDS(file.path(base,"work","03_main_objects.rds"))
fit_subset<-function(host_keep=NULL,virus_keep=NULL,label="x"){
 e<-copy(o$formal_edges);if(!is.null(host_keep))e<-e[host_index%in%host_keep];if(!is.null(virus_keep))e<-e[virus_index%in%virus_keep];hk<-if(is.null(host_keep)) sort(unique(e$host_index)) else sort(host_keep);mask<-o$pair_i%in%hk&o$pair_j%in%hk;if(sum(mask)<100||nrow(e)<10)return(data.table(analysis=label,beta=NA_real_));vv<-sort(unique(e$virus_index));X<-sparseMatrix(i=e$host_index,j=match(e$virus_index,vv),x=1,dims=c(nrow(o$host),length(vv)));S<-tcrossprod(X);sh<-as.integer(as.numeric(S[cbind(o$pair_i[mask],o$pair_j[mask])])>0);sm<-data.table(taxonomy_distance=o$dist[mask],share=sh)[,.(total_pairs=.N,positive_pairs=sum(share)),by=taxonomy_distance];fit_binomial_aggregated(sm,label)
}
set.seed(SEED+300000L);hs<-sort(unique(o$formal_edges$host_index));vs<-sort(unique(o$formal_edges$virus_index));out<-list();k<-0L
for(i in seq_len(FIG3_SENSITIVITY_ITER)){k<-k+1L;out[[k]]<-fit_subset(sample(hs,max(2L,floor(length(hs)*FIG3_SUBSAMPLE_FRACTION))),NULL,paste0("host70_iter",i));k<-k+1L;out[[k]]<-fit_subset(NULL,sample(vs,max(2L,floor(length(vs)*FIG3_SUBSAMPLE_FRACTION))),paste0("virus70_iter",i))}
hdeg<-o$formal_edges[,.N,by=host_index];vdeg<-o$formal_edges[,.N,by=virus_index];hcut<-quantile(hdeg$N,.99);vcut<-quantile(vdeg$N,.99);k<-k+1L;out[[k]]<-fit_subset(hdeg[N<hcut,host_index],NULL,"exclude_top1pct_host_degree");k<-k+1L;out[[k]]<-fit_subset(NULL,vdeg[N<vcut,virus_index],"exclude_top1pct_virus_degree")
res<-rbindlist(out,fill=TRUE);write_tsv(res,file.path(base,"tables","10_global_decay_sensitivity.tsv"));sumr<-res[,.(median_beta=median(beta,na.rm=TRUE),q025=quantile(beta,.025,na.rm=TRUE),q975=quantile(beta,.975,na.rm=TRUE),valid_n=sum(is.finite(beta))),by=.(scenario=sub("_iter[0-9]+$","",analysis))];write_tsv(sumr,file.path(base,"tables","11_global_decay_sensitivity_summary.tsv"));
# Categorical-distance sensitivity: no linear-distance assumption; distance 0 is reference.
ps <- fread(file.path(base,"tables","01_formal_taxid_host_sharing_by_distance.tsv"),sep="\t")
ps[,distance_factor:=factor(taxonomy_distance,levels=0:6)]
catfit <- glm(cbind(positive_pairs,total_pairs-positive_pairs)~distance_factor,family=binomial(),data=ps)
ct <- coef_table(catfit,"categorical_taxonomy_distance_sensitivity")
write_tsv(ct,file.path(base,"tables","12_categorical_distance_sensitivity_model.tsv"))
cat("Fig3 sensitivity PASS\n")
