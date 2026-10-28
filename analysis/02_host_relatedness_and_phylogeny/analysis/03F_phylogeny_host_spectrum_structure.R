suppressPackageStartupMessages({library(data.table); library(ape)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"), "00_config", "00_config.R"))
source(file.path(ANALYSIS_ROOT, "00_config", "00_helpers.R"))

base <- file.path(ANALYSIS_ROOT, "02_host_relatedness_and_phylogeny")
make_dirs(base)
seqm <- fread_any(SEQ_FINAL, sep="\t", showProgress=FALSE)
assoc <- fread_any(ASSOC_FINAL, sep="\t", showProgress=FALSE)

# 12 family-specific marker trees shipped under data/analysis_data/phylogeny/twelve_family/.
# The marker used per family is recorded in twelve_family/03K_D_tree_input_manifest.tsv.
families <- c(
  "Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae",
  "Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae",
  "Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"
)
perm_iter <- as.integer(Sys.getenv("FIG3_TREE_PERM_ITER", as.character(FIG3_TREE_PERM_ITER)))
max_pair_test <- as.integer(Sys.getenv("FIG3_TREE_MAX_PAIR_TEST", "50000"))
seed0 <- as.integer(Sys.getenv("FIG3_TREE_SEED", as.character(SEED)))

first_col <- function(dt, candidates) {
  z <- candidates[candidates %in% names(dt)]
  if (length(z)) z[1] else NA_character_
}
getv <- function(dt, col) if (!is.na(col) && col %in% names(dt)) std_text(dt[[col]]) else rep("", nrow(dt))

acc_col <- first_col(seqm, c("accession","accession_id","sequence_accession"))
taxid_col <- first_col(seqm, c("virus_taxid","virus_tax_id"))
species_col <- first_col(seqm, c("virus_species","current_ictv_species"))
genus_col <- first_col(seqm, c("virus_genus","current_ictv_genus"))
source_col <- first_col(seqm, c("source_label_4way","source_provenance","data_source","source"))
vc80_col <- first_col(seqm, c("vc80_final_id","VC80_id","vc80_id","vc80_cluster_id","cluster80_id","cluster80"))

seq_meta <- data.table(
  sequence_id=std_text(seqm$sequence_id),
  accession=getv(seqm, acc_col), virus_taxid=getv(seqm, taxid_col),
  virus_species=getv(seqm, species_col), virus_genus=getv(seqm, genus_col),
  source_label=getv(seqm, source_col), VC80_id=getv(seqm, vc80_col)
)
seq_meta[!present(accession), accession := sequence_id]
seq_meta <- unique(seq_meta, by="sequence_id")

host_cols <- c("host_taxid","host_species","host_genus","host_family","host_order","host_class","host_phylum","host_kingdom","country","continent")
for (z in host_cols) if (!(z %in% names(assoc))) assoc[, (z) := ""]
assoc_use <- unique(assoc[, c("sequence_id",host_cols), with=FALSE])
assoc_use[, sequence_id := std_text(sequence_id)]

all_ann <- list(); all_long <- list(); all_mem <- list(); all_pair <- list(); all_perm <- list(); all_sum <- list(); qc <- list()

sample_pairs <- function(n, max_n, seed) {
  total <- n*(n-1)/2
  if (total <= max_n) {
    ij <- which(upper.tri(matrix(FALSE,n,n)), arr.ind=TRUE)
    return(data.table(i=ij[,1], j=ij[,2], sampling_mode="all_pairs", total_pair_n=total))
  }
  set.seed(seed)
  # Ensure every tip is represented at least once, then fill uniformly with unique random pairs.
  ring <- data.table(i=seq_len(n), j=c(2:n,1L))
  ring[, `:=`(a=pmin(i,j), b=pmax(i,j))]; ring <- unique(ring[a<b,.(i=a,j=b)])
  out <- ring
  while (nrow(out) < max_n) {
    need <- max_n - nrow(out)
    x <- sample.int(n, max(need*3L, 1000L), replace=TRUE)
    y <- sample.int(n, length(x), replace=TRUE)
    add <- data.table(i=pmin(x,y), j=pmax(x,y))[i<j]
    out <- unique(rbind(out, add), by=c("i","j"))
  }
  out <- out[seq_len(max_n)]
  out[, `:=`(sampling_mode=paste0("deterministic_random_pair_sample_max_",max_n), total_pair_n=total)]
  out
}

jaccard_matrix <- function(hostsets) {
  lev <- sort(unique(unlist(hostsets, use.names=FALSE)))
  n <- length(hostsets)
  if (!length(lev)) return(matrix(NA_real_, n, n))
  H <- matrix(0, nrow=n, ncol=length(lev), dimnames=list(names(hostsets),lev))
  for (ii in seq_len(n)) if (length(hostsets[[ii]])) H[ii, match(hostsets[[ii]],lev)] <- 1
  inter <- tcrossprod(H)
  breadth <- rowSums(H)
  uni <- outer(breadth,breadth,"+") - inter
  out <- 1 - inter/uni
  out[uni==0] <- NA_real_
  diag(out) <- 0
  out
}

for (fi in seq_along(families)) {
  fam <- families[fi]
  # Released layout: trees live in <FAMILY_TREE_DIR>/trees/<family>/<family>.treefile
  # and the AA95 membership tables in <FAMILY_TREE_DIR>/AA95_membership/. The older
  # flat layout (<FAMILY_TREE_DIR>/<family>.<suffix>) is still tried first so the
  # cluster layout keeps working unchanged.
  # The released AA95 membership tables take precedence over the older flat
  # <family>.nr95.membership.tsv, which is a coarser 95%-identity clustering that
  # does not correspond one-to-one with the AA95 tree tips. Both are present in
  # the release, so the order here decides correctness, not just convenience.
  mem_candidates <- c(
    file.path(FAMILY_TREE_DIR, "AA95_membership", paste0(fam, ".AA95_membership.tsv")),
    file.path(FAMILY_TREE_DIR, "AA95_membership", paste0(fam, ".Rep_conserved.AA95_membership.tsv")),
    file.path(FAMILY_TREE_DIR, paste0(fam, ".AA95_membership.tsv")),
    file.path(FAMILY_TREE_DIR, paste0(fam, ".nr95.membership.tsv"))
  )
  tree_candidates <- c(
    file.path(FAMILY_TREE_DIR, paste0(fam, ".treefile")),
    file.path(FAMILY_TREE_DIR, "trees", fam, paste0(fam, ".treefile"))
  )
  memfile <- mem_candidates[file.exists(mem_candidates)][1]
  treefile <- tree_candidates[file.exists(tree_candidates)][1]
  if (is.na(memfile) || is.na(treefile)) {
    stop(
      "Missing 03E output for ", fam, ".\n",
      "  membership tried: ", paste(mem_candidates, collapse = "\n                    "), "\n",
      "  treefile tried:   ", paste(tree_candidates, collapse = "\n                    ")
    )
  }
  mem <- fread(memfile, sep="\t", showProgress=FALSE)
  mem[, virus_family := fam]
  all_mem[[fam]] <- mem
  tr <- read.tree(treefile)
  if (anyDuplicated(tr$tip.label)) stop("Duplicated tree tip labels in ",fam)
  tr <- ladderize(tr, right=FALSE)
  tips <- tr$tip.label
  mem_tree <- mem[representative_sequence_id %in% tips]
  if (uniqueN(mem_tree$representative_sequence_id) != length(tips)) stop("Some tree tips lack CD-HIT membership for ",fam)

  # Provenance long table: tree tip -> CD-HIT95 members -> all host associations of the cluster members.
  long <- merge(mem_tree, seq_meta, by.x="member_sequence_id", by.y="sequence_id", all.x=TRUE)
  long <- merge(long, assoc_use, by.x="member_sequence_id", by.y="sequence_id", all.x=TRUE, allow.cartesian=TRUE)
  long <- unique(long)
  # The released AA95 membership tables carry cluster id / representative / member
  # but not the percent identity to the representative, which the cluster-side
  # CD-HIT output had. It is provenance decoration only and no downstream step
  # reads it, so keep the column present and empty rather than fail the run.
  if (!"identity_to_representative_percent" %in% names(long)) {
    long[, identity_to_representative_percent := NA_real_]
  }
  setcolorder(long, c("virus_family","cdhit95_cluster_id","representative_sequence_id","member_sequence_id","is_representative","identity_to_representative_percent",
                     "accession","virus_taxid","virus_species","virus_genus","source_label","VC80_id",host_cols))
  all_long[[fam]] <- long

  # Per-tip visualization master. All host/geographic traits are aggregated from every cluster member.
  ann <- long[, .(
    cluster_member_n=uniqueN(member_sequence_id),
    cluster_accession_n=uniqueN(accession[present(accession)]),
    cluster_virus_taxid_n=uniqueN(virus_taxid[present(virus_taxid)]),
    cluster_genus_n=uniqueN(virus_genus[present(virus_genus)]),
    member_sequence_ids=collapse_unique(member_sequence_id),
    member_accessions=collapse_unique(accession),
    virus_taxids=collapse_unique(virus_taxid),
    virus_species_set=collapse_unique(virus_species),
    virus_genus_set=collapse_unique(virus_genus),
    source_label_set=collapse_unique(source_label),
    VC80_n=uniqueN(VC80_id[present(VC80_id)]),
    host_taxid_n=uniqueN(host_taxid[present(host_taxid)]),
    host_species_n=uniqueN(host_species[present(host_species)]),
    host_family_n=uniqueN(host_family[present(host_family)]),
    host_order_n=uniqueN(host_order[present(host_order)]),
    host_families=collapse_unique(host_family),
    host_orders=collapse_unique(host_order),
    dominant_host_family=mode_value(host_family),
    dominant_host_order=mode_value(host_order),
    country_n=uniqueN(country[present(country)]), continent_n=uniqueN(continent[present(continent)]),
    countries=collapse_unique(country), continents=collapse_unique(continent)
  ), by=.(virus_family,cdhit95_cluster_id,representative_sequence_id)]
  repmeta <- seq_meta[sequence_id %in% tips, .(representative_sequence_id=sequence_id, rep_accession=accession,
                                               rep_virus_taxid=virus_taxid, rep_virus_species=virus_species,
                                               rep_virus_genus=virus_genus, rep_source_label=source_label)]
  ann <- merge(ann, repmeta, by="representative_sequence_id", all.x=TRUE)
  hosttax <- long[complete_host_taxonomy(long), .(max_host_distance=max_host_distance_one_group(.SD)),
                  by=representative_sequence_id,
                  .SDcols=c("host_taxid","host_species","host_genus","host_family","host_order","host_class","host_phylum","host_kingdom")]
  ann <- merge(ann,hosttax,by="representative_sequence_id",all.x=TRUE); ann[is.na(max_host_distance),max_host_distance:=0]
  ann[, plot_order := match(representative_sequence_id, tips)]
  ann[, host_annotated := host_family_n > 0]
  setorder(ann,plot_order)
  all_ann[[fam]] <- ann

  # Primary 03F test: patristic distance vs host-family Jaccard distance.
  testtips <- ann[host_family_n>0, representative_sequence_id]
  if (length(testtips) < 4) stop("Too few host-annotated tree tips for 03F in ",fam)
  tr2 <- keep.tip(tr, testtips)
  tip_order <- tr2$tip.label
  hs <- lapply(tip_order, function(x) {
    z <- unique(std_text(long[representative_sequence_id==x & present(host_family),host_family])); z[z!=""]
  }); names(hs) <- tip_order
  Dp <- cophenetic.phylo(tr2)
  Dh <- jaccard_matrix(hs)
  breadth <- vapply(hs,length,integer(1))
  pairs <- sample_pairs(length(tip_order), max_pair_test, seed0+fi*100)
  pairs[, tip1 := tip_order[i]]; pairs[, tip2 := tip_order[j]]
  pairs[, patristic_distance := Dp[cbind(i,j)]]
  pairs[, host_jaccard_distance := Dh[cbind(i,j)]]
  pairs[, host_breadth_1 := breadth[i]]; pairs[, host_breadth_2 := breadth[j]]
  pairs[, host_breadth_difference := abs(host_breadth_1-host_breadth_2)]
  keep <- is.finite(pairs$patristic_distance) & is.finite(pairs$host_jaccard_distance)
  pairs <- pairs[keep]
  obs_rho <- suppressWarnings(cor(pairs$patristic_distance,pairs$host_jaccard_distance,method="spearman"))
  obs_breadth <- suppressWarnings(cor(pairs$patristic_distance,pairs$host_breadth_difference,method="spearman"))
  if (!is.finite(obs_rho)) stop("Non-finite observed host Jaccard correlation for ",fam)

  set.seed(seed0+fi*1000)
  pr <- numeric(perm_iter); pbr <- numeric(perm_iter)
  ii <- pairs$i; jj <- pairs$j
  for (b in seq_len(perm_iter)) {
    p <- sample.int(length(tip_order))
    yy <- Dh[cbind(p[ii],p[jj])]
    bd <- abs(breadth[p[ii]]-breadth[p[jj]])
    pr[b] <- suppressWarnings(cor(pairs$patristic_distance,yy,method="spearman"))
    pbr[b] <- suppressWarnings(cor(pairs$patristic_distance,bd,method="spearman"))
  }
  perm <- data.table(virus_family=fam,iteration=seq_len(perm_iter),perm_host_jaccard_rho=pr,perm_host_breadth_difference_rho=pbr)
  all_perm[[fam]] <- perm
  p_greater <- (sum(pr >= obs_rho,na.rm=TRUE)+1)/(sum(is.finite(pr))+1)
  p_two <- (sum(abs(pr) >= abs(obs_rho),na.rm=TRUE)+1)/(sum(is.finite(pr))+1)
  pb_greater <- (sum(pbr >= obs_breadth,na.rm=TRUE)+1)/(sum(is.finite(pbr))+1)
  summary <- data.table(
    virus_family=fam, tree_tip_n=length(tr$tip.label), host_annotated_tip_n=length(testtips),
    total_possible_pair_n=length(testtips)*(length(testtips)-1)/2, tested_pair_n=nrow(pairs),
    pair_sampling_mode=unique(pairs$sampling_mode)[1],
    observed_spearman_rho_host_jaccard=obs_rho,
    null_median_rho=median(pr,na.rm=TRUE), null_ci_low=quantile(pr,.025,na.rm=TRUE), null_ci_high=quantile(pr,.975,na.rm=TRUE),
    empirical_p_greater=p_greater, empirical_p_two_sided_abs=p_two,
    observed_spearman_rho_breadth_difference=obs_breadth,
    breadth_null_median_rho=median(pbr,na.rm=TRUE), breadth_empirical_p_greater=pb_greater,
    permutations=perm_iter
  )
  all_sum[[fam]] <- summary
  pairs[, virus_family := fam]
  pairs[, `:=`(i=NULL,j=NULL)]
  setcolorder(pairs,c("virus_family","tip1","tip2","patristic_distance","host_jaccard_distance","host_breadth_1","host_breadth_2","host_breadth_difference","sampling_mode","total_pair_n"))
  all_pair[[fam]] <- pairs
  qc[[fam]] <- data.table(virus_family=fam,tree_tip_n=length(tr$tip.label),membership_cluster_n=uniqueN(mem_tree$cdhit95_cluster_id),
                          membership_member_n=nrow(mem_tree),host_annotated_tip_n=length(testtips),tested_pair_n=nrow(pairs),permutations=perm_iter)
}

write_tsv(rbindlist(all_mem,use.names=TRUE,fill=TRUE), file.path(base,"tables","03E_cluster95_membership_all.tsv"))
write_tsv(rbindlist(all_ann,use.names=TRUE,fill=TRUE), file.path(base,"tables","03E_tree_tip_annotation_master.tsv"))
write_tsv(rbindlist(all_long,use.names=TRUE,fill=TRUE), file.path(base,"tables","22_03E_cluster_member_host_long.tsv"))
write_tsv(rbindlist(all_pair,use.names=TRUE,fill=TRUE), file.path(base,"tables","23_03F_phylogeny_host_pairwise_tested.tsv"))
write_tsv(rbindlist(all_perm,use.names=TRUE,fill=TRUE), file.path(base,"tables","24_03F_phylogeny_host_permutation_long.tsv"))
write_tsv(rbindlist(all_sum,use.names=TRUE,fill=TRUE), file.path(base,"tables","25_03F_phylogeny_host_distance_summary.tsv"))
write_tsv(rbindlist(qc,use.names=TRUE,fill=TRUE), file.path(base,"qc","03_Fig3_EF_QC.tsv"))
cat("Fig3 03F phylogeny-host spectrum structure PASS\n")
print(rbindlist(all_sum,use.names=TRUE,fill=TRUE))
