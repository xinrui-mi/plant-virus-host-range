suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"));source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"06_TSWV_evolution");make_dirs(base);assoc<-fread_any(ASSOC_FINAL,sep="\t",showProgress=FALSE)
manifest<-fread(file.path(base,"pipeline_results","03_alignment_and_phylogeny","07E_five_gene_tree_manifest.tsv"),sep="\t")
out<-list();k<-0L
for(g in c("RdRp","GnGc","NSm","N","NSs")){
 mp<-fread(file.path(base,"pipeline_results","02_exact_haplotypes",paste0(g,".haplotype_map.tsv")),sep="\t")
 reps<-mp[is_representative=="TRUE",representative_sequence_id]
 a<-assoc[sequence_id%in%reps]
 k<-k+1L;out[[k]]<-data.table(gene=g,tree_haplotype_n=length(reps),represented_accession_n=uniqueN(mp$member_sequence_id),host_family_n=uniqueN(a$host_family[present(a$host_family)]),host_order_n=uniqueN(a$host_order[present(a$host_order)]),country_n=uniqueN(a$country[present(a$country)]),continent_n=uniqueN(a$continent[present(a$continent)]))
}
z<-merge(rbindlist(out),manifest,by="gene",all.x=TRUE);write_tsv(z,file.path(base,"tables","06_TSWV_five_gene_evolutionary_input_summary.tsv"));cat("Fig7 five-gene metadata summary PASS\n");print(z)
