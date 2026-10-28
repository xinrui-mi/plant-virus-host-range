suppressPackageStartupMessages({library(data.table);library(Matrix);library(parallel)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"))
source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"02_host_relatedness_and_phylogeny");make_dirs(base)
o<-readRDS(file.path(base,"work","03_main_objects.rds"))

randomize_edges <- function(edge_dt,n_swaps,seed) {
  set.seed(seed); h<-edge_dt$host_index; v<-edge_dt$virus_index; ne<-length(h)
  env<-new.env(hash=TRUE,parent=emptyenv()); key<-function(a,b) paste0(a,"/",b)
  for (k in seq_len(ne)) env[[key(h[k],v[k])]]<-TRUE
  done<-0L; att<-0L; maxatt<-max(1000L,n_swaps*30L)
  while(done<n_swaps && att<maxatt){
    att<-att+1L; z<-sample.int(ne,2L);e1<-z[1];e2<-z[2];h1<-h[e1];h2<-h[e2];v1<-v[e1];v2<-v[e2]
    if(h1==h2||v1==v2)next
    k1<-key(h1,v2);k2<-key(h2,v1)
    if(exists(k1,env,inherits=FALSE)||exists(k2,env,inherits=FALSE))next
    rm(list=c(key(h1,v1),key(h2,v2)),envir=env);env[[k1]]<-TRUE;env[[k2]]<-TRUE;v[e1]<-v2;v[e2]<-v1;done<-done+1L
  }
  list(h=h,v=v,swaps_done=done,attempts=att)
}

run_null <- function(edges,nvirus,pair_mask,iter,label,seed0){
  pi<-o$pair_i[pair_mask];pj<-o$pair_j[pair_mask];dd<-o$dist[pair_mask];nh<-nrow(o$host);nswap<-min(FIG3_SWAP_MULTIPLIER*nrow(edges),FIG3_MAX_SWAPS)
  worker<-function(it){
    rr<-randomize_edges(edges,nswap,seed0+it)
    X<-sparseMatrix(i=rr$h,j=rr$v,x=1,dims=c(nh,nvirus));S<-tcrossprod(X);sh<-as.integer(as.numeric(S[cbind(pi,pj)])>0)
    sm<-data.table(taxonomy_distance=dd,share=sh)[,.(total_pairs=.N,positive_pairs=sum(share),sharing_probability=mean(share)),by=taxonomy_distance][order(taxonomy_distance)]
    fit<-fit_binomial_aggregated(sm,paste0(label,"_null_",it));fit[,`:=`(iteration=it,swaps_requested=nswap,swaps_done=rr$swaps_done,attempts=rr$attempts)]
    sm[,`:=`(iteration=it,swaps_requested=nswap,swaps_done=rr$swaps_done)]
    list(fit=fit,dist=sm)
  }
  # Iterations run in small parallel chunks, and every chunk is validated before it
  # is accepted. `mclapply` is not trusted on its own: in constrained containers a
  # fork worker can be reaped without returning a value, and mclapply then silently
  # substitutes NULL, so a run can report fewer iterations than it was asked for
  # while still exiting successfully. Missing iterations are recomputed serially, and
  # the run aborts rather than emitting a short null distribution.
  ncore <- max(1L, min(NCORES, iter))
  block_retry <- 0L
  chunk <- function(ids) {
    res <- mclapply(ids, function(it) try(worker(it), silent = TRUE),
                    mc.cores = min(ncore, length(ids)), mc.preschedule = TRUE)
    ok <- vapply(res, function(r) !inherits(r, "try-error") &&
                   is.data.frame(r$fit) && nrow(r$fit) == 1L, logical(1))
    if (any(!ok)) {
      message("  ", label, ": ", sum(!ok), " iteration(s) not delivered by fork; recomputing serially")
      block_retry <<- block_retry + sum(!ok)
      for (j in which(!ok)) {
        for (try_k in seq_len(3L)) {
          r <- try(worker(ids[j]), silent = TRUE)
          if (!inherits(r, "try-error") && is.data.frame(r$fit) && nrow(r$fit) == 1L) { res[[j]] <- r; break }
        }
        if (inherits(res[[j]], "try-error") || !is.data.frame(res[[j]]$fit))
          stop("iteration ", ids[j], " failed 3 serial retries; refusing to report a short null")
      }
    }
    res
  }
  # Memory: the `dist` table is ~7 rows per iteration; keeping all 999 of them plus the
  # fit rows is a ~50 MB object that fragments badly on a constrained heap. The per
  # iteration `dist` rows are therefore concatenated straight into the block accumulator
  # and the fit rows are the only per-iteration structure retained.
  blocks <- split(seq_len(iter), ceiling(seq_len(iter) / max(ncore * 10L, 20L)))
  fit_parts <- vector("list", length(blocks))
  dist_parts <- vector("list", length(blocks))
  n_retry <- 0L
  for (b in seq_along(blocks)) {
    part <- chunk(blocks[[b]])
    fit_parts[[b]] <- rbindlist(lapply(part, `[[`, "fit"), fill = TRUE)
    dist_parts[[b]] <- rbindlist(lapply(part, `[[`, "dist"), fill = TRUE)
    part <- NULL
    n_retry <- n_retry + block_retry
    cat("  ", label, "iterations", sum(vapply(fit_parts[seq_len(b)], nrow, integer(1))), "/", iter,
        if (n_retry > 0L) paste0(" (", n_retry, " recovered serially after fork loss)") else "",
        "\n"); flush.console()
  }
  fit_dt <- rbindlist(fit_parts, fill = TRUE); fit_parts <- NULL
  dist_dt <- rbindlist(dist_parts, fill = TRUE); dist_parts <- NULL
  if (uniqueN(fit_dt$iteration) != iter)
    stop("null for ", label, " produced ", uniqueN(fit_dt$iteration), " distinct iterations; expected ", iter)
  invisible(gc())
  list(fit = fit_dt, dist = dist_dt, n_retry = n_retry)
}

summarize_null <- function(z,observed_beta,prefix){
  f<-z$fit; d<-z$dist
  p<-(sum(f$beta<=observed_beta,na.rm=TRUE)+1)/(sum(is.finite(f$beta))+1)
  bsum<-data.table(observed_beta=observed_beta,null_iter=sum(is.finite(f$beta)),null_beta_median=median(f$beta,na.rm=TRUE),null_beta_q025=quantile(f$beta,.025,na.rm=TRUE),null_beta_q975=quantile(f$beta,.975,na.rm=TRUE),empirical_p_more_negative=p,median_swap_completion=median(f$swaps_done/pmax(f$swaps_requested,1),na.rm=TRUE))
  dsum<-d[,.(null_iteration_n=uniqueN(iteration),null_probability_mean=mean(sharing_probability),null_probability_median=median(sharing_probability),null_probability_q025=quantile(sharing_probability,.025),null_probability_q975=quantile(sharing_probability,.975)),by=taxonomy_distance][order(taxonomy_distance)]
  write_tsv(f,file.path(base,"tables",paste0(prefix,"_null_beta_iterations.tsv")))
  write_tsv(d,file.path(base,"tables",paste0(prefix,"_null_distance_iterations.tsv")))
  write_tsv(bsum,file.path(base,"tables",paste0(prefix,"_null_beta_summary.tsv")))
  write_tsv(dsum,file.path(base,"tables",paste0(prefix,"_null_distance_summary.tsv")))
  bsum
}

formal_mask<-rep(TRUE,length(o$pair_i))
formal<-run_null(o$formal_edges,o$formal_nvirus,formal_mask,FIG3_NULL_ITER,"formal",SEED+100000L)
fs<-summarize_null(formal,o$formal_observed$beta[1],"06_formal_taxid_degree_preserving")
rs<-NULL
if (FIG3_RDRP_NULL_ITER>0 && nrow(o$rdrp_edges)>1) {
  re<-copy(o$rdrp_edges);setnames(re,"rdrp_index","virus_index")
  rz<-run_null(re,o$rdrp_nvirus,o$rdrp_pair_mask,FIG3_RDRP_NULL_ITER,"RdRP65",SEED+200000L)
  rs<-summarize_null(rz,o$rdrp_observed$beta[1],"10_RdRP65_degree_preserving")
}
qc<-rbind(data.table(network="formal_taxid",iterations_requested=FIG3_NULL_ITER,iterations_completed=fs$null_iter,median_swap_completion=fs$median_swap_completion),if(!is.null(rs)) data.table(network="RdRP65",iterations_requested=FIG3_RDRP_NULL_ITER,iterations_completed=rs$null_iter,median_swap_completion=rs$median_swap_completion) else NULL)
write_tsv(qc,file.path(base,"qc","02_Fig3_null_QC.tsv"));cat("Fig3 degree-preserving null PASS\n");print(qc)
flush.console();invisible(gc())
# Long jobs are launched with Rscript. An explicit non-interactive exit avoids
# LSF retaining an already-finished R parent after fork-worker cleanup.
if (!interactive()) quit(save="no",status=0,runLast=FALSE)
