suppressPackageStartupMessages({library(data.table);library(Matrix)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_config.R"))
source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base <- file.path(ANALYSIS_ROOT,"02_host_relatedness_and_phylogeny"); make_dirs(base)
assoc <- fread_any(ASSOC_FINAL,sep="\t",showProgress=FALSE)
assoc[, complete_host := complete_host_taxonomy(.SD)]
a <- assoc[complete_host & present(virus_taxid)]
# one current taxonomy row per host taxid
host_cols <- c("host_taxid","host_species","host_genus","host_family","host_order","host_class","host_phylum","host_kingdom")
host <- unique(a[,..host_cols])
conf <- host[, .N, by=host_taxid][N>1]
if (nrow(conf)) stop("Host taxid taxonomy conflict in the released association input")
host <- host[order(host_taxid)]; host[,host_index:=.I]
hmap <- host[,.(host_taxid,host_index)]
# formal NCBI Virus TaxID unit only
edges <- unique(a[,.(host_taxid,virus_taxid)])
edges <- merge(edges,hmap,by="host_taxid")
vlevels <- sort(unique(edges$virus_taxid)); vmap <- data.table(virus_taxid=vlevels,virus_index=seq_along(vlevels))
edges <- merge(edges,vmap,by="virus_taxid")
if (nrow(host) < 50 || nrow(edges)<100) stop("Fig3 formal taxid network unexpectedly small")
ij <- combn(seq_len(nrow(host)),2L)
pair_i <- ij[1,]; pair_j <- ij[2,]
dist <- host_distance_vectors(host,pair_i,pair_j)
if (any(is.na(dist))) stop("NA taxonomy distance despite complete host taxonomy")
X <- sparseMatrix(i=edges$host_index,j=edges$virus_index,x=1,dims=c(nrow(host),length(vlevels)))
S <- tcrossprod(X)
shared_n <- as.numeric(S[cbind(pair_i,pair_j)])
share <- as.integer(shared_n>0)
pair_summary <- data.table(taxonomy_distance=dist,share=share,shared_n=shared_n)[,.(total_pairs=.N,positive_pairs=sum(share),sharing_probability=mean(share),mean_shared_virus_taxid=mean(shared_n)),by=taxonomy_distance][order(taxonomy_distance)]
write_tsv(pair_summary,file.path(base,"tables","01_formal_taxid_host_sharing_by_distance.tsv"))
obs <- fit_binomial_aggregated(pair_summary,"formal_NCBI_Virus_TaxID")
write_tsv(obs,file.path(base,"tables","02_formal_taxid_global_decay_model.tsv"))
# family-specific slopes. TaxIDs with ambiguous family are excluded from family-specific models only.
taxfam <- a[current_ictv_family_eligible(virus_ictv_match_rank,virus_family),.(virus_family_n=uniqueN(virus_family),virus_family=mode_value(virus_family)),by=virus_taxid]
valid_tf <- taxfam[virus_family_n==1]
family_res <- list(); fi <- 0L
for (fam in sort(unique(valid_tf$virus_family))) {
  vids <- valid_tf[virus_family==fam,virus_taxid]
  ee <- edges[virus_taxid %in% vids]
  hkeep <- sort(unique(ee$host_index)); vt_n <- uniqueN(ee$virus_taxid)
  if (vt_n < 10L || length(hkeep)<20L) next
  mask <- pair_i %in% hkeep & pair_j %in% hkeep
  if (sum(mask)<100L) next
  vv <- sort(unique(ee$virus_index)); xx <- sparseMatrix(i=ee$host_index,j=match(ee$virus_index,vv),x=1,dims=c(nrow(host),length(vv)))
  ss <- tcrossprod(xx); sh <- as.integer(as.numeric(ss[cbind(pair_i[mask],pair_j[mask])])>0)
  sm <- data.table(taxonomy_distance=dist[mask],share=sh)[,.(total_pairs=.N,positive_pairs=sum(share)),by=taxonomy_distance]
  if (sum(sm$positive_pairs)<20L) next
  fi <- fi+1L; z <- fit_binomial_aggregated(sm,paste0("family:",fam)); z[,`:=`(virus_family=fam,virus_taxid_n=vt_n,host_taxid_n=length(hkeep))]; family_res[[fi]]<-z
}
family_res <- rbindlist(family_res,fill=TRUE); write_tsv(family_res,file.path(base,"tables","03_family_specific_decay_models.tsv"))
# RdRP65 operational RNA lineage extension
r <- assoc[complete_host & present(rdrp_operational_cluster_id),.(host_taxid,rdrp_operational_cluster_id)] |> unique()
r <- merge(r,hmap,by="host_taxid")
rlevels <- sort(unique(r$rdrp_operational_cluster_id)); rmap <- data.table(rdrp_operational_cluster_id=rlevels,rdrp_index=seq_along(rlevels)); r<-merge(r,rmap,by="rdrp_operational_cluster_id")
rhosts <- sort(unique(r$host_index)); rmask <- pair_i %in% rhosts & pair_j %in% rhosts
RX <- sparseMatrix(i=r$host_index,j=r$rdrp_index,x=1,dims=c(nrow(host),length(rlevels))); RS <- tcrossprod(RX)
rshare <- as.integer(as.numeric(RS[cbind(pair_i[rmask],pair_j[rmask])])>0)
rsm <- data.table(taxonomy_distance=dist[rmask],share=rshare)[,.(total_pairs=.N,positive_pairs=sum(share),sharing_probability=mean(share)),by=taxonomy_distance][order(taxonomy_distance)]
rmod <- fit_binomial_aggregated(rsm,"RdRP65_operational_RNA_lineage")
write_tsv(rsm,file.path(base,"tables","04_RdRP65_host_sharing_by_distance.tsv"));write_tsv(rmod,file.path(base,"tables","05_RdRP65_global_decay_model.tsv"))
# save compact objects for null/sensitivity
obj <- list(host=host,pair_i=pair_i,pair_j=pair_j,dist=dist,formal_edges=edges[,.(host_index,virus_index)],formal_nvirus=length(vlevels),formal_observed=obs,rdrp_edges=r[,.(host_index,rdrp_index)],rdrp_nvirus=length(rlevels),rdrp_pair_mask=rmask,rdrp_observed=rmod)
saveRDS(obj,file.path(base,"work","03_main_objects.rds"),compress=TRUE)
filter_qc <- data.table(metric=c("association_rows_total","association_rows_complete_host_and_formal_taxid","host_taxid_n","formal_virus_taxid_n","host_pairs_n","RdRP65_lineage_n","RdRP65_host_pairs_n","family_slope_n"),value=c(nrow(assoc),nrow(a),nrow(host),length(vlevels),length(pair_i),length(rlevels),sum(rmask),nrow(family_res)))
write_tsv(filter_qc,file.path(base,"qc","01_Fig3_main_QC.tsv")); cat("Fig3 main PASS\n"); print(filter_qc)
