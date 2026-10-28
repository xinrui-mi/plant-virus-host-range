suppressPackageStartupMessages({library(data.table); library(ape)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_config.R"))
source(file.path(ANALYSIS_ROOT, "00_config", "00_helpers.R"))

# Fig.3G sensitivity analyses for the completed 03E-03F family phylogenies.
# This script NEVER rebuilds a tree and NEVER overwrites 03E/03F outputs.
# It evaluates:
#   A) 50k pair-sampling stability (100 repeated pair samples)
#   B) support-stratified host-spectrum permutations (999)
#   C) high-support-only subsets (cluster_accession_n >=2 and >=3)
#   D) alignment-filter bias, with special attention to Potyviridae

base <- file.path(ANALYSIS_ROOT, "02_host_relatedness_and_phylogeny")
make_dirs(base)
# 12 family-specific marker trees (see 03F and twelve_family/03K_D_tree_input_manifest.tsv).
families <- c(
  "Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
  "Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
  "Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
)
seed0 <- as.integer(Sys.getenv("FIG3_03G_SEED", as.character(SEED + 3000L)))
max_pair <- as.integer(Sys.getenv("FIG3_03G_MAX_PAIR", "50000"))
pair_reps <- as.integer(Sys.getenv("FIG3_03G_PAIR_REPS", "100"))
perm_iter <- as.integer(Sys.getenv("FIG3_03G_PERM_ITER", "999"))
min_stratum <- as.integer(Sys.getenv("FIG3_03G_MIN_STRATUM", "20"))
max_strata <- as.integer(Sys.getenv("FIG3_03G_MAX_STRATA", "8"))

ann_all <- fread(file.path(base,"tables","03E_tree_tip_annotation_master.tsv"), sep="\t", showProgress=FALSE)
long_all <- fread(file.path(base,"tables","22_03E_cluster_member_host_long.tsv"), sep="\t", showProgress=FALSE)
mem_all <- fread(file.path(base,"tables","03E_cluster95_membership_all.tsv"), sep="\t", showProgress=FALSE)
pair03f <- fread(file.path(base,"tables","23_03F_phylogeny_host_pairwise_tested.tsv"), sep="\t", showProgress=FALSE)
seqm <- fread_any(SEQ_FINAL, sep="\t", showProgress=FALSE)
assoc <- fread_any(ASSOC_FINAL, sep="\t", showProgress=FALSE)

need <- c("virus_family","representative_sequence_id","cluster_accession_n","host_family_n")
if (!all(need %in% names(ann_all))) stop("03G requires completed 03F annotation master; missing: ", paste(setdiff(need,names(ann_all)),collapse=","))
if (!("host_family" %in% names(long_all))) stop("03G requires host_family in 22_03E_cluster_member_host_long.tsv")

sample_pairs <- function(n, max_n, seed) {
  total <- n*(n-1)/2
  if (total <= max_n) {
    ij <- which(upper.tri(matrix(FALSE,n,n)), arr.ind=TRUE)
    return(data.table(i=ij[,1],j=ij[,2],mode="all_pairs",total_pair_n=total))
  }
  set.seed(seed)
  # Seed with a ring so every tip is represented, then fill with unique random pairs.
  ring <- data.table(i=seq_len(n),j=c(2:n,1L))
  ring[,`:=`(a=pmin(i,j),b=pmax(i,j))]
  out <- unique(ring[a<b,.(i=a,j=b)])
  while (nrow(out) < max_n) {
    needn <- max_n - nrow(out)
    x <- sample.int(n,max(needn*3L,1000L),replace=TRUE)
    y <- sample.int(n,length(x),replace=TRUE)
    add <- data.table(i=pmin(x,y),j=pmax(x,y))[i<j]
    out <- unique(rbind(out,add),by=c("i","j"))
  }
  out <- out[seq_len(max_n)]
  out[,`:=`(mode=paste0("random_pair_sample_",max_n),total_pair_n=total)]
  out
}

jaccard_matrix <- function(hostsets) {
  lev <- sort(unique(unlist(hostsets,use.names=FALSE)))
  n <- length(hostsets)
  if (!length(lev)) return(matrix(NA_real_,n,n))
  H <- matrix(0,nrow=n,ncol=length(lev))
  for (ii in seq_len(n)) if (length(hostsets[[ii]])) H[ii,match(hostsets[[ii]],lev)] <- 1
  inter <- tcrossprod(H)
  breadth <- rowSums(H)
  uni <- outer(breadth,breadth,"+") - inter
  z <- 1-inter/uni
  z[uni==0] <- NA_real_
  diag(z) <- 0
  z
}

make_support_strata <- function(support, min_size=20L, max_bins=8L) {
  x <- log1p(pmax(0,as.numeric(support)))
  n <- length(x)
  kmax <- min(max_bins, max(1L, floor(n/min_size)))
  if (kmax < 2L) return(factor(rep("S1",n)))
  for (k in seq(kmax,2L,by=-1L)) {
    br <- unique(as.numeric(quantile(x, probs=seq(0,1,length.out=k+1L), type=8, na.rm=TRUE)))
    if (length(br) < 3L) next
    # Extend the outer bounds slightly so every finite value is included.
    br[1] <- br[1]-1e-12; br[length(br)] <- br[length(br)]+1e-12
    g <- cut(x,breaks=br,include.lowest=TRUE,labels=FALSE)
    if (anyNA(g)) next
    tab <- table(g)
    if (length(tab)>=2L && min(tab)>=min_size) return(factor(paste0("S",g),levels=paste0("S",sort(unique(g)))))
  }
  factor(rep("S1",n))
}

# The release keeps trees and AA95 membership in subdirectories of
# FAMILY_TREE_DIR; the flat <family>.<suffix> layout still comes first so a
# cluster tree behaves exactly as before.
fam_tree_path <- function(fam) {
  cand <- c(file.path(FAMILY_TREE_DIR, paste0(fam, ".treefile")),
            file.path(FAMILY_TREE_DIR, "trees", fam, paste0(fam, ".treefile")))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Treefile not found for ", fam, ". Tried:\n  ", paste(cand, collapse = "\n  "))
  hit[1]
}
fam_membership_path <- function(fam) {
  cand <- c(
    file.path(FAMILY_TREE_DIR, paste0(fam, ".nr95.membership.tsv")),
    file.path(FAMILY_TREE_DIR, "AA95_membership", paste0(fam, ".AA95_membership.tsv")),
    file.path(FAMILY_TREE_DIR, "AA95_membership", paste0(fam, ".Rep_conserved.AA95_membership.tsv"))
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Membership table not found for ", fam, ". Tried:\n  ", paste(cand, collapse = "\n  "))
  hit[1]
}
fam_alignment_qc_path <- function(fam) {
  cand <- c(
    file.path(FAMILY_TREE_DIR, paste0(fam, ".alignment_sequence_QC.tsv")),
    file.path(ANALYSIS_ROOT, "02_host_relatedness_and_phylogeny", "pipeline_results",
              "03K_family_specific_markers", "06_AA95_alignment_QC", fam,
              paste0(fam, ".alignment_sequence_QC.tsv")),
    file.path(ANALYSIS_ROOT, "02_host_relatedness_and_phylogeny", "pipeline_results",
              "03K_family_specific_markers", "07_Geminiviridae_Rep_conserved",
              "AA95_alignment_QC", "Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv")
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Alignment QC not found for ", fam, ". Tried:\n  ", paste(cand, collapse = "\n  "))
  hit[1]
}

prepare_family <- function(fam, subset_min_support=1L) {
  tr <- read.tree(fam_tree_path(fam))
  ann <- ann_all[virus_family==fam & host_family_n>0 & cluster_accession_n>=subset_min_support]
  tips <- intersect(tr$tip.label,ann$representative_sequence_id)
  if (length(tips)<4L) stop("Too few tips in ",fam," at support threshold ",subset_min_support)
  tr <- keep.tip(tr,tips)
  tip_order <- tr$tip.label
  long <- long_all[virus_family==fam & representative_sequence_id %in% tip_order]
  hs <- lapply(tip_order,function(x) sort(unique(std_text(long[representative_sequence_id==x & present(host_family),host_family]))))
  hs <- lapply(hs,function(x) x[x!=""])
  names(hs) <- tip_order
  support <- ann$cluster_accession_n[match(tip_order,ann$representative_sequence_id)]
  if (anyNA(support)) stop("Missing cluster_accession_n for ",fam)
  list(tree=tr,tips=tip_order,Dp=cophenetic.phylo(tr),Dh=jaccard_matrix(hs),support=as.integer(support),hostsets=hs)
}

rho_for_pairs <- function(x,pairs) {
  pd <- x$Dp[cbind(pairs$i,pairs$j)]
  hd <- x$Dh[cbind(pairs$i,pairs$j)]
  ok <- is.finite(pd)&is.finite(hd)
  suppressWarnings(cor(pd[ok],hd[ok],method="spearman"))
}

# A. Repeated pair-sampling stability -----------------------------------------
pair_iter <- list(); pair_sum <- list()
for (fi in seq_along(families)) {
  fam <- families[fi]; x <- prepare_family(fam,1L)
  total <- length(x$tips)*(length(x$tips)-1)/2
  reps <- if (total<=max_pair) 1L else pair_reps
  zz <- rbindlist(lapply(seq_len(reps),function(r) {
    pp <- sample_pairs(length(x$tips),max_pair,seed0+fi*100000L+r)
    data.table(virus_family=fam,replicate=r,total_possible_pair_n=total,tested_pair_n=nrow(pp),sampling_mode=unique(pp$mode)[1],rho=rho_for_pairs(x,pp))
  }))
  pair_iter[[fam]] <- zz
  pair_sum[[fam]] <- zz[,.(
    replicate_n=.N,
    rho_median=median(rho,na.rm=TRUE),rho_mean=mean(rho,na.rm=TRUE),rho_sd=sd(rho,na.rm=TRUE),
    rho_q025=quantile(rho,.025,na.rm=TRUE),rho_q25=quantile(rho,.25,na.rm=TRUE),rho_q75=quantile(rho,.75,na.rm=TRUE),rho_q975=quantile(rho,.975,na.rm=TRUE),
    rho_min=min(rho,na.rm=TRUE),rho_max=max(rho,na.rm=TRUE)
  ),by=.(virus_family,total_possible_pair_n,tested_pair_n,sampling_mode)]
}
write_tsv(rbindlist(pair_iter),file.path(base,"tables","26_03G_pair_sampling_stability_iterations.tsv"))
write_tsv(rbindlist(pair_sum),file.path(base,"tables","26_03G_pair_sampling_stability_summary.tsv"))

# B. Support-stratified permutations ----------------------------------------
strata_rows <- list(); sp_iter <- list(); sp_sum <- list()
for (fi in seq_along(families)) {
  fam <- families[fi]; x <- prepare_family(fam,1L)
  # Hold the exact 03F tested pair set fixed; only the annotation permutation changes.
  pp0 <- pair03f[virus_family==fam]
  ii <- match(pp0$tip1,x$tips); jj <- match(pp0$tip2,x$tips)
  if(anyNA(ii)||anyNA(jj)) stop("03G cannot map the original 03F pair set for ",fam)
  pp <- data.table(i=ii,j=jj)
  pd <- as.numeric(pp0$patristic_distance); hd <- as.numeric(pp0$host_jaccard_distance)
  ok <- is.finite(pd)&is.finite(hd); pp<-pp[ok]; pd<-pd[ok]; hd<-hd[ok]
  obs <- suppressWarnings(cor(pd,hd,method="spearman"))
  strata <- make_support_strata(x$support,min_stratum,max_strata)
  ss <- data.table(virus_family=fam,tip=x$tips,cluster_accession_n=x$support,stratum=as.character(strata))
  strata_rows[[fam]] <- ss[,.(tip_n=.N,support_min=min(cluster_accession_n),support_median=median(cluster_accession_n),support_max=max(cluster_accession_n)),by=.(virus_family,stratum)]
  set.seed(seed0+fi*1000L)
  pr <- numeric(perm_iter)
  groups <- split(seq_along(x$tips),strata)
  for (b in seq_len(perm_iter)) {
    map <- seq_along(x$tips)
    for (g in groups) if (length(g)>1L) map[g] <- sample(g,length(g),replace=FALSE)
    yy <- x$Dh[cbind(map[pp$i],map[pp$j])]
    pr[b] <- suppressWarnings(cor(pd,yy,method="spearman"))
  }
  sp_iter[[fam]] <- data.table(virus_family=fam,iteration=seq_len(perm_iter),perm_rho=pr)
  sp_sum[[fam]] <- data.table(
    virus_family=fam,tree_tip_n=length(x$tips),support_stratum_n=length(unique(strata)),tested_pair_n=nrow(pp),
    observed_rho=obs,null_median_rho=median(pr,na.rm=TRUE),null_ci_low=quantile(pr,.025,na.rm=TRUE),null_ci_high=quantile(pr,.975,na.rm=TRUE),
    empirical_p_greater=(sum(pr>=obs,na.rm=TRUE)+1)/(sum(is.finite(pr))+1),
    empirical_p_two_sided_abs=(sum(abs(pr)>=abs(obs),na.rm=TRUE)+1)/(sum(is.finite(pr))+1),permutations=perm_iter
  )
}
write_tsv(rbindlist(strata_rows),file.path(base,"tables","27_03G_support_strata.tsv"))
write_tsv(rbindlist(sp_iter),file.path(base,"tables","27_03G_support_stratified_permutation_iterations.tsv"))
write_tsv(rbindlist(sp_sum),file.path(base,"tables","27_03G_support_stratified_permutation_summary.tsv"))

# C. High-support-only sensitivity ------------------------------------------
hs_iter <- list(); hs_sum <- list(); kk <- 0L
for (fi in seq_along(families)) for (thr in c(2L,3L)) {
  fam <- families[fi]
  annn <- ann_all[virus_family==fam & host_family_n>0 & cluster_accession_n>=thr,.N]
  if (annn<20L) {
    kk<-kk+1L; hs_sum[[kk]]<-data.table(virus_family=fam,min_cluster_accession_n=thr,tip_n=annn,status="SKIP_lt20_tips")
    next
  }
  x <- prepare_family(fam,thr); pp <- sample_pairs(length(x$tips),max_pair,seed0+fi*100L+thr)
  pd <- x$Dp[cbind(pp$i,pp$j)]; hd <- x$Dh[cbind(pp$i,pp$j)]; ok<-is.finite(pd)&is.finite(hd);pp<-pp[ok];pd<-pd[ok];hd<-hd[ok]
  obs <- suppressWarnings(cor(pd,hd,method="spearman"))
  set.seed(seed0+fi*10000L+thr)
  pr <- numeric(perm_iter)
  for (b in seq_len(perm_iter)) {
    map<-sample.int(length(x$tips)); yy<-x$Dh[cbind(map[pp$i],map[pp$j])]
    pr[b]<-suppressWarnings(cor(pd,yy,method="spearman"))
  }
  kk<-kk+1L
  hs_iter[[kk]]<-data.table(virus_family=fam,min_cluster_accession_n=thr,iteration=seq_len(perm_iter),perm_rho=pr)
  hs_sum[[kk]]<-data.table(
    virus_family=fam,min_cluster_accession_n=thr,tip_n=length(x$tips),tested_pair_n=nrow(pp),status="PASS",
    observed_rho=obs,null_median_rho=median(pr,na.rm=TRUE),null_ci_low=quantile(pr,.025,na.rm=TRUE),null_ci_high=quantile(pr,.975,na.rm=TRUE),
    empirical_p_greater=(sum(pr>=obs,na.rm=TRUE)+1)/(sum(is.finite(pr))+1),
    empirical_p_two_sided_abs=(sum(abs(pr)>=abs(obs),na.rm=TRUE)+1)/(sum(is.finite(pr))+1),permutations=perm_iter
  )
}
if (length(hs_iter)) write_tsv(rbindlist(hs_iter,fill=TRUE),file.path(base,"tables","28_03G_high_support_permutation_iterations.tsv"))
write_tsv(rbindlist(hs_sum,fill=TRUE),file.path(base,"tables","28_03G_high_support_sensitivity_summary.tsv"))

# D. Alignment-QC exclusion-bias audit --------------------------------------
first_col <- function(dt,candidates) { z<-candidates[candidates%in%names(dt)]; if(length(z)) z[1] else NA_character_ }
getv <- function(dt,col) if(!is.na(col)&&col%in%names(dt)) std_text(dt[[col]]) else rep("",nrow(dt))
acc_col <- first_col(seqm,c("accession","accession_id","sequence_accession"))
genus_col <- first_col(seqm,c("virus_genus","current_ictv_genus"))
source_col <- first_col(seqm,c("source_label_4way","source_provenance","data_source","source"))
seq_meta <- data.table(sequence_id=std_text(seqm$sequence_id),accession=getv(seqm,acc_col),virus_genus=getv(seqm,genus_col),source_label=getv(seqm,source_col))
seq_meta[!present(accession),accession:=sequence_id];seq_meta<-unique(seq_meta,by="sequence_id")
if (!("host_family"%in%names(assoc))) assoc[,host_family:=""]
assoc_small <- unique(assoc[,.(sequence_id=std_text(sequence_id),host_family=std_text(host_family))])

bias_detail<-list();bias_sum<-list();bias_genus<-list()
for (fam in families) {
  qf<-fam_alignment_qc_path(fam); mf<-fam_membership_path(fam)
  q<-fread(qf,sep="\t",showProgress=FALSE);setnames(q,"tip","representative_sequence_id")
  m<-fread(mf,sep="\t",showProgress=FALSE)
  z<-merge(m,seq_meta,by.x="member_sequence_id",by.y="sequence_id",all.x=TRUE)
  z<-merge(z,assoc_small,by.x="member_sequence_id",by.y="sequence_id",all.x=TRUE,allow.cartesian=TRUE)
  cl<-z[,.(
    cluster_member_n=uniqueN(member_sequence_id),cluster_accession_n=uniqueN(accession[present(accession)]),host_family_n=uniqueN(host_family[present(host_family)]),
    virus_genus_set=collapse_unique(virus_genus),source_label_set=collapse_unique(source_label),
    any_landflora=any(is_true_flag(getv(seqm,"source_LandFlora")),na.rm=TRUE)
  ),by=.(virus_family,cdhit95_cluster_id,representative_sequence_id)]
  zz<-merge(cl,q,by="representative_sequence_id",all.x=TRUE)
  zz[,drop_gap:=gap_fraction>0.50];zz[,drop_effective:=effective_aa<100]
  zz[,drop_reason:=fcase(action=="keep","keep",drop_gap&drop_effective,"gap_and_effective",drop_gap,"gap",drop_effective,"effective_aa",default="other")]
  bias_detail[[fam]]<-zz
  keep<-zz[action=="keep"];drop<-zz[action=="drop"]
  wp1<-if(nrow(keep)>0&&nrow(drop)>0) suppressWarnings(wilcox.test(keep$cluster_accession_n,drop$cluster_accession_n,exact=FALSE)$p.value) else NA_real_
  wp2<-if(nrow(keep)>0&&nrow(drop)>0) suppressWarnings(wilcox.test(keep$host_family_n,drop$host_family_n,exact=FALSE)$p.value) else NA_real_
  fp<-NA_real_
  if(nrow(keep)>0&&nrow(drop)>0) {
    tb<-table(factor(zz$action,levels=c("keep","drop")),factor(zz$any_landflora,levels=c(FALSE,TRUE)))
    if(all(dim(tb)==c(2,2))) fp<-fisher.test(tb)$p.value
  }
  bias_sum[[fam]]<-data.table(
    virus_family=fam,nr95_cluster_n=nrow(zz),kept_n=nrow(keep),dropped_n=nrow(drop),retention_fraction=nrow(keep)/nrow(zz),
    drop_gap_only_n=zz[action=="drop"&drop_gap&!drop_effective,.N],drop_effective_only_n=zz[action=="drop"&!drop_gap&drop_effective,.N],drop_both_n=zz[action=="drop"&drop_gap&drop_effective,.N],
    kept_median_cluster_accession_n=if(nrow(keep)) median(keep$cluster_accession_n) else NA_real_,dropped_median_cluster_accession_n=if(nrow(drop)) median(drop$cluster_accession_n) else NA_real_,support_wilcox_p=wp1,
    kept_median_host_family_n=if(nrow(keep)) median(keep$host_family_n) else NA_real_,dropped_median_host_family_n=if(nrow(drop)) median(drop$host_family_n) else NA_real_,host_breadth_wilcox_p=wp2,
    kept_landflora_fraction=if(nrow(keep)) mean(keep$any_landflora) else NA_real_,dropped_landflora_fraction=if(nrow(drop)) mean(drop$any_landflora) else NA_real_,landflora_fisher_p=fp
  )
  # Descriptive genus-level retention; no inferential claim, because cluster genus sets can contain >1 genus.
  gx<-zz[,.(cluster_n=.N,dropped_n=sum(action=="drop"),retention_fraction=mean(action=="keep")),by=.(virus_family,virus_genus_set)][cluster_n>=5]
  setorder(gx,retention_fraction,-cluster_n);bias_genus[[fam]]<-gx
}
write_tsv(rbindlist(bias_detail,fill=TRUE),file.path(base,"tables","29_03G_alignment_filter_bias.tsv"))
write_tsv(rbindlist(bias_sum,fill=TRUE),file.path(base,"tables","29_03G_alignment_filter_bias_summary.tsv"))
write_tsv(rbindlist(bias_genus,fill=TRUE),file.path(base,"tables","29_03G_alignment_filter_by_genus.tsv"))

# Compact QC index -----------------------------------------------------------
ps<-rbindlist(pair_sum);sp<-rbindlist(sp_sum);hs<-rbindlist(hs_sum,fill=TRUE);ab<-rbindlist(bias_sum)
qc<-rbindlist(list(
  ps[,.(module="pair_sampling_stability",virus_family,metric="rho_sd",value=as.character(rho_sd))],
  sp[,.(module="support_stratified_permutation",virus_family,metric="empirical_p_two_sided_abs",value=as.character(empirical_p_two_sided_abs))],
  hs[status=="PASS",.(module=paste0("high_support_ge",min_cluster_accession_n),virus_family,metric="observed_rho",value=as.character(observed_rho))],
  ab[,.(module="alignment_filter",virus_family,metric="retention_fraction",value=as.character(retention_fraction))]
),fill=TRUE)
write_tsv(qc,file.path(base,"qc","04_Fig3_03G_sensitivity_QC.tsv"))

cat("Fig3 03G phylogeny-host spectrum sensitivity PASS\n")
cat("\nPair-sampling stability:\n");print(ps)
cat("\nSupport-stratified permutation:\n");print(sp)
cat("\nHigh-support subsets:\n");print(hs)
cat("\nAlignment-filter bias summary:\n");print(ab)
