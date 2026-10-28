suppressPackageStartupMessages({library(data.table);library(ape)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"))
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"06_TSWV_evolution");make_dirs(base)
# Frozen L/M/S segment trees and matched isolate metadata are shipped as released
# derived analysis inputs; they are read directly from the repository data tree.
d <- TSWV_TREE_DIR
meta<-fread_any(TSWV_LMS_META,sep="\t",showProgress=FALSE)
segs<-c("L","M","S"); trees<-setNames(lapply(segs,function(s) read.tree(file.path(d,paste0(s,".treefile")))),segs)
sets<-lapply(trees,function(t) sort(t$tip.label));if(length(unique(vapply(sets,paste,collapse="|",FUN.VALUE="")))!=1) stop("L/M/S trees do not contain identical isolate labels")
tips<-sets[[1]];n<-length(tips);if(n<20) stop("Too few matched isolates for segment concordance")
perm_iter<-as.integer(Sys.getenv("TSWV_SEGMENT_CONCORDANCE_PERM","999"));seed0<-as.integer(Sys.getenv("TSWV_SEGMENT_SEED","20260813"))

parse_support<-function(x){
  if(length(x)==0||is.na(x)||!nzchar(x)) return(c(SH=NA_real_,UF=NA_real_))
  z<-regmatches(x,gregexpr("[0-9]+(?:\\.[0-9]+)?",x,perl=TRUE))[[1]];z<-suppressWarnings(as.numeric(z))
  if(length(z)>=2) c(SH=z[1],UF=z[2]) else c(SH=NA_real_,UF=if(length(z)==1)z[1] else NA_real_)
}
canon_key<-function(side,univ){
  side<-sort(side);other<-sort(setdiff(univ,side))
  if(length(side)<length(other)) use<-side else if(length(other)<length(side)) use<-other else use<-if(paste(side,collapse="|")<=paste(other,collapse="|"))side else other
  paste(use,collapse="|")
}
extract_splits<-function(tr,seg){
  nt<-length(tr$tip.label); all<-tr$tip.label
  children<-split(tr$edge[,2],tr$edge[,1])
  cache<-new.env(parent=emptyenv())
  desc<-function(node){
    k<-as.character(node);if(exists(k,cache,inherits=FALSE))return(get(k,cache))
    if(node<=nt){z<-tr$tip.label[node]} else {ch<-children[[k]];z<-sort(unique(unlist(lapply(ch,desc),use.names=FALSE)))}
    assign(k,z,cache);z
  }
  rows<-list();ii<-0L
  for(e in seq_len(nrow(tr$edge))){
    child<-tr$edge[e,2];if(child<=nt)next
    side<-desc(child);if(length(side)<2||length(side)>nt-2)next
    sup<-parse_support(tr$node.label[child-nt]);key<-canon_key(side,all)
    can_tips<-strsplit(key,"\\|",fixed=FALSE)[[1]]
    ii<-ii+1L;rows[[ii]]<-data.table(segment=seg,split_key=key,split_size=length(can_tips),split_tips=paste(can_tips,collapse=";"),SH_aLRT=sup['SH'],UFBoot=sup['UF'])
  }
  z<-rbindlist(rows,fill=TRUE)
  # The arbitrary Newick root can yield the same unrooted split twice; retain the best-supported representation.
  if(nrow(z)){ z[,sort_SH:=ifelse(is.finite(SH_aLRT),SH_aLRT,-Inf)]; z[,sort_UF:=ifelse(is.finite(UFBoot),UFBoot,-Inf)]; setorder(z,-sort_SH,-sort_UF); z<-z[!duplicated(split_key)]; z[,c("sort_SH","sort_UF"):=NULL] }
  z[,high_support:=is.finite(SH_aLRT)&is.finite(UFBoot)&SH_aLRT>=80&UFBoot>=95]
  z
}
compatible_split<-function(a,b,univ){
  A<-strsplit(a,"\\|",fixed=FALSE)[[1]];B<-strsplit(b,"\\|",fixed=FALSE)[[1]];Ac<-setdiff(univ,A);Bc<-setdiff(univ,B)
  any(c(length(intersect(A,B)),length(intersect(A,Bc)),length(intersect(Ac,B)),length(intersect(Ac,Bc)))==0)
}

splits<-lapply(segs,function(s)extract_splits(trees[[s]],s));names(splits)<-segs
all_split<-rbindlist(splits,fill=TRUE);fwrite(all_split,file.path(base,"tables","12_TSWV_segment_tree_splits.tsv"),sep="\t",na="")

pairs<-list(c("L","M"),c("L","S"),c("M","S"));sumrows<-list();confrows<-list();isolate_rows<-list()
for(pi in seq_along(pairs)){
  a<-pairs[[pi]][1];b<-pairs[[pi]][2];ta<-trees[[a]];tb<-trees[[b]]
  # reorder distances to the exact same label order
  Da<-cophenetic.phylo(ta)[tips,tips];Db<-cophenetic.phylo(tb)[tips,tips]
  lo<-lower.tri(Da);obsrho<-suppressWarnings(cor(Da[lo],Db[lo],method="spearman"))
  set.seed(seed0+pi*1000L);pr<-numeric(perm_iter)
  for(k in seq_len(perm_iter)){p<-sample.int(n);pr[k]<-suppressWarnings(cor(Da[lo],Db[p,p][lo],method="spearman"))}
  ka<-unique(splits[[a]]$split_key);kb<-unique(splits[[b]]$split_key)
  rf<-length(setdiff(ka,kb))+length(setdiff(kb,ka));rf_norm<-rf/(2*(n-3))
  ha<-splits[[a]][high_support==TRUE];hb<-splits[[b]][high_support==TRUE]
  shared_hs<-length(intersect(ha$split_key,hb$split_key));union_hs<-length(union(ha$split_key,hb$split_key))
  cc<-list();ci<-0L
  if(nrow(ha)&&nrow(hb)) for(i in seq_len(nrow(ha))) for(j in seq_len(nrow(hb))) if(!compatible_split(ha$split_key[i],hb$split_key[j],tips)){
    ci<-ci+1L;cc[[ci]]<-data.table(segment_pair=paste0(a,"_vs_",b),segment_A=a,segment_B=b,
      A_split_key=ha$split_key[i],A_split_size=ha$split_size[i],A_SH_aLRT=ha$SH_aLRT[i],A_UFBoot=ha$UFBoot[i],A_split_tips=ha$split_tips[i],
      B_split_key=hb$split_key[j],B_split_size=hb$split_size[j],B_SH_aLRT=hb$SH_aLRT[j],B_UFBoot=hb$UFBoot[j],B_split_tips=hb$split_tips[j])
  }
  if(length(cc))confrows[[pi]]<-rbindlist(cc)
  sumrows[[pi]]<-data.table(segment_pair=paste0(a,"_vs_",b),matched_isolate_n=n,
    patristic_spearman_rho=obsrho,patristic_null_median=median(pr,na.rm=TRUE),patristic_null_ci_low=quantile(pr,.025,na.rm=TRUE),patristic_null_ci_high=quantile(pr,.975,na.rm=TRUE),
    patristic_empirical_p_greater=(sum(pr>=obsrho,na.rm=TRUE)+1)/(sum(is.finite(pr))+1),permutations=perm_iter,
    RF_distance=rf,normalized_RF_distance=rf_norm,
    high_support_split_A_n=nrow(ha),high_support_split_B_n=nrow(hb),shared_high_support_split_n=shared_hs,
    high_support_split_jaccard=if(union_hs)shared_hs/union_hs else NA_real_,incompatible_high_support_split_pair_n=ci)
  # Per-isolate neighborhood concordance. This is exploratory ranking, not a reassortment test.
  for(t in tips){
    others<-setdiff(tips,t);rho_i<-suppressWarnings(cor(Da[t,others],Db[t,others],method="spearman"))
    kn<-function(D,k){z<-sort(D[t,others],index.return=TRUE);others[z$ix[seq_len(min(k,length(others)))]]}
    a5<-kn(Da,5);b5<-kn(Db,5);a10<-kn(Da,10);b10<-kn(Db,10)
    isolate_rows[[length(isolate_rows)+1L]]<-data.table(paired_id=t,segment_pair=paste0(a,"_vs_",b),distance_profile_spearman=rho_i,
      nearest5_jaccard=length(intersect(a5,b5))/length(union(a5,b5)),nearest10_jaccard=length(intersect(a10,b10))/length(union(a10,b10)))
  }
}
summary<-rbindlist(sumrows);fwrite(summary,file.path(base,"tables","13_TSWV_segment_concordance_summary.tsv"),sep="\t",na="")
conf<-if(length(confrows))rbindlist(confrows,fill=TRUE) else data.table(segment_pair=character())
fwrite(conf,file.path(base,"tables","14_TSWV_high_support_split_conflicts.tsv"),sep="\t",na="")
iso_long<-rbindlist(isolate_rows);iso<-iso_long[,.(median_distance_profile_spearman=median(distance_profile_spearman,na.rm=TRUE),min_distance_profile_spearman=min(distance_profile_spearman,na.rm=TRUE),
  median_nearest5_jaccard=median(nearest5_jaccard,na.rm=TRUE),median_nearest10_jaccard=median(nearest10_jaccard,na.rm=TRUE)),by=paired_id]
iso<-merge(iso,meta,by="paired_id",all.x=TRUE);q05<-quantile(iso$median_distance_profile_spearman,.05,na.rm=TRUE);iso[,exploratory_low_concordance_5pct:=median_distance_profile_spearman<=q05]
setorder(iso,median_distance_profile_spearman);iso[,discordance_rank:=.I]
fwrite(iso_long,file.path(base,"tables","15_TSWV_segment_isolate_concordance_long.tsv"),sep="\t",na="")
fwrite(iso,file.path(base,"tables","15_TSWV_segment_isolate_discordance_ranking.tsv"),sep="\t",na="")
qc<-rbindlist(list(
  data.table(metric="matched_isolate_n",value=n),
  data.table(metric="L_high_support_split_n",value=splits$L[high_support==TRUE,.N]),
  data.table(metric="M_high_support_split_n",value=splits$M[high_support==TRUE,.N]),
  data.table(metric="S_high_support_split_n",value=splits$S[high_support==TRUE,.N]),
  summary[,.(metric=paste0(segment_pair,"::normalized_RF_distance"),value=normalized_RF_distance)],
  summary[,.(metric=paste0(segment_pair,"::incompatible_high_support_split_pair_n"),value=incompatible_high_support_split_pair_n)]
),fill=TRUE)
fwrite(qc,file.path(base,"qc","06_TSWV_segment_concordance_QC.tsv"),sep="\t")
cat("TSWV segment concordance/reassortment analysis PASS\n")
print(summary)
cat("\nInterpretation guardrail: high-support incompatible segment splits are evidence of segment phylogenetic discordance consistent with reassortment; they are not, by themselves, proof of a specific reassortment event.\n")
