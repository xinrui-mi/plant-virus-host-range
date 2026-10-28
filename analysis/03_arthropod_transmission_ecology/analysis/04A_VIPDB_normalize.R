suppressPackageStartupMessages({library(data.table);library(readxl)})
source(file.path(Sys.getenv("ANALYSIS_ROOT"),"00_config","00_config.R"))
source(file.path(ANALYSIS_ROOT,"00_config","00_helpers.R"))
base<-file.path(ANALYSIS_ROOT,"03_arthropod_transmission_ecology");make_dirs(base)

# The VIPDB raw workbook is a third-party resource that is not redistributed with
# this release. When it is not available (VIPDB_RAW unset or file missing), the
# released normalized VIPDB outputs -- which are byte-identical to the frozen
# manuscript run -- are restored into <module>/tables so that every downstream
# module remains runnable. No normalization logic is altered.
if (!nzchar(VIPDB_RAW) || !file.exists(VIPDB_RAW)) {
  message("VIPDB raw workbook not available (", if (nzchar(VIPDB_RAW)) VIPDB_RAW else "VIPDB_RAW unset",
          "); restoring the released normalized VIPDB outputs.")
  restore <- function(src, dest) {
    if (!file.exists(src)) stop("Released VIPDB input missing: ", src)
    y <- fread_any(src, sep="\t", showProgress=FALSE)
    write_tsv(y, dest); y
  }
  norm <- restore(VIPDB_NORMALIZED, file.path(base,"tables","01_VIPDB_relation_normalized.tsv"))
  vs   <- restore(VIPDB_VIRUS_LEVEL, file.path(base,"tables","02_VIPDB_virus_level_summary.tsv"))
  audit <- norm[,.(relation_n=.N,virus_taxid_n=uniqueN(virus_taxid[present(virus_taxid)])),
                by=.(raw_transmission_mode,fine_mode,fine4_mode,retention_group)][order(-relation_n)]
  write_tsv(audit, file.path(base,"tables","03_VIPDB_transmission_mode_audit.tsv"))
  cat("Fig4 VIPDB normalize: released normalized outputs restored (n relations = ", nrow(norm), ")\n", sep="")
  quit(status=0)
}

x<-as.data.table(readxl::read_excel(VIPDB_RAW,sheet="Sheet1",col_types="text"))
need<-c("Virus Family","Virus Genus","Virus Name","Virus TaxID","Vector Order","Vector Family","Vector Genus","Vector","Vector TaxID","Virus Host","Host TaxID","DOI for Validation Virus-Insect Relationships","Virus Transmission Mode","DOI for Validation Virus Transmission Mode","Lineage")
if(length(setdiff(need,names(x))))stop("VIPDB missing columns: ",paste(setdiff(need,names(x)),collapse=", "))
normalize_mode_key<-function(z){z<-as.character(z);z[is.na(z)]<-"";z<-gsub("\u00A0"," ",z,fixed=TRUE);z<-gsub("[‐‑‒–—−]","-",z,perl=TRUE);z<-gsub("[[:punct:]]+"," ",z);z<-trimws(gsub("[[:space:]]+"," ",tolower(z)));z}
raw<-std_text(x[["Virus Transmission Mode"]]);key<-normalize_mode_key(raw)
allowed<-c("circulative persistent non propagative transmission","circulative persistent propagative transmission","non persistent transmission","non circulative semi persistent transmission","persistent circulative","na","")
unknown<-sort(unique(key[!key%in%allowed]));if(length(unknown)){write_tsv(unique(data.table(raw_value=raw[key%in%unknown],normalized_value=key[key%in%unknown])),file.path(base,"qc","VIPDB_UNMAPPED_TRANSMISSION_MODE.tsv"));stop("Unknown VIPDB transmission mode after normalization")}
fine<-fcase(key=="non persistent transmission","Non-persistent",key=="non circulative semi persistent transmission","Semi-persistent",key=="circulative persistent non propagative transmission","Persistent non-propagative",key=="circulative persistent propagative transmission","Persistent propagative",key=="persistent circulative","Persistent circulative (propagation unspecified)",default=NA_character_)
fine4<-fifelse(fine%in%c("Non-persistent","Semi-persistent","Persistent non-propagative","Persistent propagative"),fine,NA_character_)
ret<-fcase(fine%in%c("Non-persistent","Semi-persistent"),"short-retention / non-circulative",fine%in%c("Persistent non-propagative","Persistent propagative","Persistent circulative (propagation unspecified)"),"long-retention / circulative",default=NA_character_)
x[,virus_taxid:=std_text(`Virus TaxID`)];x[,vector_taxid:=std_text(`Vector TaxID`)]
norm<-x[,.(virus_family=std_text(`Virus Family`),virus_genus=std_text(`Virus Genus`),virus_name=std_text(`Virus Name`),virus_taxid,vector_order=std_text(`Vector Order`),vector_family=std_text(`Vector Family`),vector_genus=std_text(`Vector Genus`),vector_name=std_text(Vector),vector_taxid,virus_host=std_text(`Virus Host`),host_taxid=std_text(`Host TaxID`),relation_doi=std_text(`DOI for Validation Virus-Insect Relationships`),raw_transmission_mode=raw,transmission_mode_key=key,fine_mode=fine,fine4_mode=fine4,retention_group=ret,transmission_doi=std_text(`DOI for Validation Virus Transmission Mode`),lineage=std_text(Lineage))]
write_tsv(norm,file.path(base,"tables","01_VIPDB_relation_normalized.tsv"))
vs<-norm[present(virus_taxid),.(relation_row_n=.N,vector_taxid_n=uniqueN(vector_taxid[present(vector_taxid)]),vector_order_n=uniqueN(vector_order[present(vector_order)]),fine_mode_set=collapse_unique(fine_mode),fine4_mode_n=uniqueN(fine4_mode[present(fine4_mode)]),fine4_mode=if(uniqueN(fine4_mode[present(fine4_mode)])==1) mode_value(fine4_mode) else "",retention_group_n=uniqueN(retention_group[present(retention_group)]),retention_group=if(uniqueN(retention_group[present(retention_group)])==1) mode_value(retention_group) else ""),by=virus_taxid]
vs[,fine4_status:=fcase(fine4_mode_n==0,"no_four_class_mode",fine4_mode_n==1,"single_four_class_mode",default="multiple_four_class_modes")]
vs[,retention_status:=fcase(retention_group_n==0,"missing",retention_group_n==1,"single_group",default="mixed_groups")]
write_tsv(vs,file.path(base,"tables","02_VIPDB_virus_level_summary.tsv"))
audit<-norm[,.(relation_n=.N,virus_taxid_n=uniqueN(virus_taxid[present(virus_taxid)])),by=.(raw_transmission_mode,fine_mode,fine4_mode,retention_group)][order(-relation_n)]
write_tsv(audit,file.path(base,"tables","03_VIPDB_transmission_mode_audit.tsv"))
counts<-norm[,.(n=.N),by=fine_mode]
expected<-data.table(fine_mode=c("Persistent non-propagative","Persistent propagative","Non-persistent","Semi-persistent","Persistent circulative (propagation unspecified)",NA_character_),expected=c(600L,284L,227L,192L,21L,51L))
chk<-merge(expected,counts,by="fine_mode",all.x=TRUE);chk[is.na(n),n:=0L];chk[,pass:=n==expected];write_tsv(chk,file.path(base,"qc","01_VIPDB_mode_count_QC.tsv"))
stopifnot(nrow(norm)==1375L,all(chk$pass),sum(present(norm$fine4_mode))==1303L,sum(norm$retention_group=="short-retention / non-circulative",na.rm=TRUE)==419L,sum(norm$retention_group=="long-retention / circulative",na.rm=TRUE)==905L)
cat("Fig4 VIPDB normalize PASS\n");print(audit)
