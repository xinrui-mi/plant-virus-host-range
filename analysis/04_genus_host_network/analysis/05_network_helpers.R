suppressPackageStartupMessages({library(data.table);library(igraph)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"))
source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))

network_metrics <- function(long_dt, virus_col="virus_node", seed=SEED, make_candidates=TRUE) {
  d <- copy(as.data.table(long_dt))
  if (virus_col != "virus_node") setnames(d, virus_col, "virus_node")
  d <- unique(d[present(virus_node)&present(host_family)&present(vc80_final_id),.(virus_node,host_family,vc80_final_id)])
  ew <- d[,.(weight=uniqueN(vc80_final_id)),by=.(virus_node,host_family)]
  if(nrow(ew)<2) stop("Network input has <2 edges")
  ed <- ew[,.(from=paste0("V:",virus_node),to=paste0("H:",host_family),weight=as.numeric(weight))]
  nodes <- data.table(name=unique(c(ed$from,ed$to))); nodes[,type:=startsWith(name,"H:")]
  g <- graph_from_data_frame(ed,directed=FALSE,vertices=nodes)
  set.seed(seed)
  clu <- cluster_louvain(g,weights=E(g)$weight)
  memb <- membership(clu); modmap <- data.table(name=names(memb),module_id=as.integer(memb))
  invw <- 1/pmax(E(g)$weight,1e-12)
  bet <- betweenness(g,directed=FALSE,weights=invw,normalized=TRUE)
  strn <- strength(g,weights=E(g)$weight)
  deg <- degree(g)
  virus_vnames <- paste0("V:",sort(unique(ew$virus_node)))
  em <- merge(ed,modmap,by.x="to",by.y="name",all.x=TRUE)
  part <- em[from%in%virus_vnames,.(w=sum(weight)),by=.(from,module_id)][,.(participation_coefficient=1-sum((w/sum(w))^2)),by=from]
  metric <- data.table(vname=virus_vnames,virus_node=sub("^V:","",virus_vnames),bipartite_betweenness=as.numeric(bet[virus_vnames]),bipartite_strength=as.numeric(strn[virus_vnames]),bipartite_degree=as.numeric(deg[virus_vnames]))
  metric <- merge(metric,part,by.x="vname",by.y="from",all.x=TRUE); metric[is.na(participation_coefficient),participation_coefficient:=0]
  breadth <- ew[,.(host_family_n=uniqueN(host_family)),by=virus_node]
  support <- d[,.(qc_pass_VC80_support=uniqueN(vc80_final_id)),by=virus_node]
  metric <- merge(metric,breadth,by="virus_node",all.x=TRUE); metric<-merge(metric,support,by="virus_node",all.x=TRUE)
  metric[,scaled_betweenness:=scale01(bipartite_betweenness)];metric[,scaled_participation:=scale01(participation_coefficient)];metric[,scaled_host_breadth:=scale01(host_family_n)];metric[,scaled_support:=scale01(qc_pass_VC80_support)]
  metric[,bridge_score:=0.6*scaled_betweenness+0.4*scaled_participation]
  metric[,candidate_score:=0.4*scaled_host_breadth+0.3*bridge_score+0.3*scaled_support]
  metric <- merge(metric,modmap[startsWith(name,"V:"),.(vname=name,module_id)],by="vname",all.x=TRUE)
  hb_raw <- safe_quantile(metric$host_family_n, FIG5_HOST_Q)
  hb_cut <- ceiling(hb_raw - 1e-12)
  b_nonzero <- metric$bridge_score[is.finite(metric$bridge_score) & metric$bridge_score > 0]
  b_cut <- if (length(b_nonzero)) safe_quantile(b_nonzero, FIG5_BRIDGE_Q) else 0
  thresholds <- data.table(
    host_breadth_quantile=FIG5_HOST_Q,
    host_breadth_cutoff_raw=hb_raw,
    host_breadth_cutoff=hb_cut,
    bridge_quantile=FIG5_BRIDGE_Q,
    bridge_cutoff=b_cut,
    support_cutoff=FIG5_SUPPORT_MIN
  )
  if(make_candidates){
    metric[,candidate := host_family_n>=thresholds$host_breadth_cutoff & bridge_score>=thresholds$bridge_cutoff & qc_pass_VC80_support>=thresholds$support_cutoff]
    setorder(metric,-candidate,-candidate_score,-host_family_n)
  } else metric[,candidate:=FALSE]
  list(metrics=metric[,vname:=NULL],edge_weights=ew,module_map=modmap,thresholds=thresholds,modularity=modularity(clu),graph=g,input_long=d)
}
