#!/usr/bin/env Rscript
# HCC9 single-script pipeline: mtDNA markers + epiAneufinder + clone-level CNV summary.
# Merged, without changing statistical definitions, from user-provided 1.R and 2.R.
# Runs only HCC9_C1 / HCC9_C2 and the same 8 mtDNA-defined clones.
# Default: HCC9_MODE=full   => Stage 1 (mtDNA + epiAneufinder) and Stage 2 (clone report).
# Fast:    HCC9_MODE=report => Stage 2 only; needs existing per-cell 0/1/2 RDS and mtDNA marker CSV.
# Output:  pic/HCC9/HCC9_Markers_3Sheets_CloneRecalculated.xlsx (3 sheets).
# Dependencies: data.table, epiAneufinder, hg38 BSgenome; openxlsx OR writexl for Excel.
# NOTE: 0=Loss, 1=Normal, 2=Gain are inferred *relative* scATAC states,
#       not independently validated absolute DNA copy numbers.

options(stringsAsFactors = FALSE, timeout = 45)
RUN_MODE <- tolower(Sys.getenv("HCC9_MODE", "full"))
if (!(RUN_MODE %in% c("full", "report"))) stop("HCC9_MODE must be full or report")
R_LIB <- Sys.getenv("HCC9_R_LIBS_USER", "/mnt/data1/lin_env/R")
if (dir.exists(R_LIB)) .libPaths(unique(c(R_LIB, .libPaths())))
suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("openxlsx", quietly=TRUE) && !requireNamespace("writexl", quietly=TRUE)) {
  stop("A 3-sheet Excel file is required: install 'openxlsx' or 'writexl' before running")
}

# Server defaults; can be overridden via environment variables for reuse from GitHub.
ROOT <- Sys.getenv("HCC9_ROOT", "/mnt/data1/linyanlin/scATAC/CNV")
OUT <- file.path(ROOT, "pic", "HCC9")
TMP <- file.path(OUT, "tmp")
FRAG_ROOT <- Sys.getenv("HCC9_FRAG_ROOT", "/mnt/hdd/linyanlin/sc-ATAC/singlecellfilter/step1_filter")
META_PATH <- Sys.getenv("HCC9_META_PATH", file.path(OUT, "clone_chrCNV_similarity", "01_matched_cells_used.csv"))
VAF_PATH <- Sys.getenv("HCC9_VAF_PATH", file.path(ROOT, "突变亚克隆", "HCC9_results_8subclone", "heteroplasmy_matrix.rds"))
V2_PATH <- Sys.getenv("HCC9_V2_PATH", file.path(OUT, "HCC9_CNV_marker_v2.csv"))
WINDOW_SIZE <- 1000000L   # 1 Mb non-overlapping bins, for relatively sparse scATAC
MIN_FRAGS <- 5000L       # preliminary; retain exclusions in QC, don't claim all cells pass
NCORES <- as.integer(Sys.getenv("HCC9_NCORES", "2"))
MIN_CONSECUTIVE_BINS <- 3L
MIN_CLONE_EVENT_RATE <- 0.50
MIN_PER_CELL_BIN_FRACTION <- 0.70
MIN_SPECIFIC_RATE <- 0.60
MAX_OTHER_RATE <- 0.25
MIN_DELTA_RATE <- 0.35
MIN_AUROC <- 0.70
VAF_POSITIVE_PERCENT <- 5

if (!is.finite(NCORES) || NCORES < 1L) stop("Invalid HCC9_NCORES")
dir.create(TMP, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(TMPDIR=TMP, TMP=TMP, TEMP=TMP, SQLITE_TMPDIR=TMP,
           R_LIBS_USER=R_LIB)
setwd(ROOT)
logmsg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep="")
check_exists <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing input(s): ", paste(missing,collapse="; "))
}
auroc <- function(x,y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x,y), ties.method="average")
  as.numeric((sum(r[seq_along(x)])-length(x)*(length(x)+1)/2)/(length(x)*length(y)))
}

logmsg("HCC9 pipeline starting; mode=", RUN_MODE)

if (RUN_MODE == "full") {
  # ================== 1. Matched metadata (source of truth) ==================

  check_exists(c(META_PATH,VAF_PATH))
  m <- fread(META_PATH)
  req <- c("patient","Clone","sample","mtDNA_Cell","cnv_cell","match_status")
  if (!all(req %in% names(m))) stop("Mapping file lacks: ",paste(setdiff(req,names(m)),collapse=","))
  m <- m[,.(patient=as.character(patient),Clone=as.integer(Clone),
           sample=as.character(sample),mtDNA_Cell=as.character(mtDNA_Cell),
           cnv_cell=as.character(cnv_cell),match_status=as.character(match_status))]
  if (anyNA(m) || !all(m$patient=="HCC9") || !all(m$match_status=="strict_1to1")) stop("Mapping contains missing/non-strict HCC9 entries")
  if (nrow(m)!=3655L || !identical(sort(unique(m$sample)), c("HCC9_C1","HCC9_C2"))) stop("Unexpected HCC9 cell/sample counts")
  if (anyDuplicated(m$cnv_cell) || anyDuplicated(m$mtDNA_Cell)) stop("Duplicate mapping cell IDs")
  m[,barcode := sub("_(HCC9_C1|HCC9_C2)$", "", cnv_cell)]
  if (!all(paste0(m$barcode,"_",m$sample)==m$cnv_cell)) stop("Unexpected CNV barcode naming rule")
  if (anyDuplicated(m[,.(sample,barcode)])) stop("Duplicated sample/barcode")
  clones <- sort(unique(m$Clone))
  if (!identical(clones,1:8)) stop("Expected HCC9 Clone1..8")
  fwrite(m,file.path(OUT,"HCC9_v4_matched_cells.csv"))
  logmsg("Matched cells: ",nrow(m),"; Clone labels: ",paste(clones,collapse=","))

  # ================== 2. New mtDNA marker tables ==================
  logmsg("Computing mtDNA marker statistics from original VAF matrix...")
  v0 <- readRDS(VAF_PATH)
  if (is.null(dim(v0)) || is.null(rownames(v0)) || is.null(colnames(v0))) stop("VAF matrix lacks dimensions or names")
  if (!all(m$mtDNA_Cell %in% colnames(v0))) stop("Not all mapped mtDNA cells are present in VAF matrix")
  v <- as.matrix(v0[,m$mtDNA_Cell,drop=FALSE]); storage.mode(v) <- "double"
  rm(v0); gc(FALSE)
  if (max(v,na.rm=TRUE)<=1) v <- v*100
  if (any(v[is.finite(v)]<0 | v[is.finite(v)]>100)) stop("Invalid VAF values (expected percent 0..100)")
  if (anyDuplicated(rownames(v))) stop("Duplicate mtDNA mutation labels")
  positive <- is.finite(v) & v>=VAF_POSITIVE_PERCENT
  valid <- is.finite(v)
  mt_out <- vector("list",length(clones)*nrow(v)); ii <- 0L
  for (k in clones) {
    which_clone <- which(m$Clone==k)
    which_other <- which(m$Clone!=k)
    for (j in seq_len(nrow(v))) {
      ii <- ii+1L
      row <- v[j,]
      pos <- positive[j,]
      npos <- sum(pos[which_clone]); npos_other <- sum(pos[which_other])
      rate <- npos/length(which_clone); rate_other <- npos_other/length(which_other)
      per_s <- list()
      peer <- numeric()
      for (s in c("HCC9_C1","HCC9_C2")) {
        ix <- which(m$Clone==k & m$sample==s)
        other <- which(m$Clone!=k & m$sample==s)
        per_s[[s]] <- list(n=length(ix),pos=sum(pos[ix]),rate=mean(pos[ix]),
                            other=mean(pos[other]),auc=auroc(row[ix],row[other]))
        for (o in setdiff(clones,k)) {
          opp <- which(m$Clone==o & m$sample==s)
          peer <- c(peer,mean(pos[ix])-mean(pos[opp]))
        }
      }
      c1 <- per_s$HCC9_C1; c2 <- per_s$HCC9_C2
      pos_vaf <- row[which_clone][pos[which_clone]]
      mlabel <- paste0("m.",rownames(v)[j])
      site <- suppressWarnings(as.integer(sub("^([0-9]+).*$","\\1",rownames(v)[j])))
      caution <- if(!is.na(site) && site>=11180 && site<=11190) "REVIEW_ADJACENT_11180_11190" else "coverage_depth_not_available"
      priority <- (npos>=5L && rate>=0.10 && rate-rate_other>=0.10 &&
                     c1$rate-c1$other>=0.10 && c2$rate-c2$other>=0.10 &&
                     min(peer,na.rm=TRUE)>=0.05 && caution=="coverage_depth_not_available")
      mt_out[[ii]] <- data.table(Clone=paste0("Clone",k),Marker=mlabel,
        Positive_cells=npos,Clone_cells=length(which_clone),Other_positive=npos_other,
        Other_cells=length(which_other),Clone_positive_pct=round(100*rate,2),
        Other_positive_pct=round(100*rate_other,2),Enrichment_pp=round(100*(rate-rate_other),2),
        C1_positive=paste0(c1$pos,"/",c1$n),C2_positive=paste0(c2$pos,"/",c2$n),
        C1_positive_pct=round(100*c1$rate,2),C2_positive_pct=round(100*c2$rate,2),
        C1_other_positive_pct=round(100*c1$other,2),C2_other_positive_pct=round(100*c2$other,2),
        C1_VAF_AUROC=round(c1$auc,3),C2_VAF_AUROC=round(c2$auc,3),
        Min_peer_enrichment_pp=round(100*min(peer,na.rm=TRUE),2),
        Mean_positive_VAF_pct=if(length(pos_vaf)) round(mean(pos_vaf),2) else NA_real_,
        VAF_nonNA_cells=sum(valid[j,which_clone]),
        Coverage_verified=FALSE,
        Caution=caution,
        Priority=if(priority) "priority_candidate" else "exploratory_or_weak",
        Interpretation="VAF>=5% is an observed VAF threshold, NOT site-depth-verified mutation status")
    }
  }
  mtall <- rbindlist(mt_out)
  setorder(mtall,Clone,-Enrichment_pp)
  fwrite(mtall,file.path(OUT,"HCC9_mtDNA_marker_all_sites_v4.csv"))
  mt <- mtall[Positive_cells>=5 & Clone_positive_pct>=10 & Enrichment_pp>=10]
  mt[,rank_priority:=as.integer(Priority=="priority_candidate")]
  setorder(mt,Clone,-rank_priority,-Enrichment_pp,-Min_peer_enrichment_pp)
  mt[,Rank_in_clone:=seq_len(.N),by=Clone]
  mt[,rank_priority:=NULL]
  fwrite(mt,file.path(OUT,"HCC9_mtDNA_marker_new.csv"))
  logmsg("mtDNA candidate markers: ",nrow(mt),"; priority: ",sum(mt$Priority=="priority_candidate"))
  rm(v,positive,valid,mtall,mt_out); gc(FALSE)

  # ================== 3. Locate the hg38 ENCODE blacklist ==================
  # The input fragments were previously aligned on hg38 in this workflow; confirm upstream
  # reference/genome build before treating these CNV states as biologically interpretable.
  blacklist_dest <- file.path(OUT,"hg38-blacklist.v2.bed")
  blacklist_sources <- c(blacklist_dest,
     file.path("/mnt/data1/lin_env/R","epiAneufinder-main","sample_data","hg38-blacklist.v2.bed"),
     file.path("/mnt/data1/lin_env/R","hg38-blacklist.v2.bed"),
     system.file("sample_data","hg38-blacklist.v2.bed",package="epiAneufinder"))
  blacklist_sources <- unique(blacklist_sources[nzchar(blacklist_sources)])
  for (f in blacklist_sources) if (file.exists(f) && file.info(f)$size>5000) {
    if (normalizePath(f)!=normalizePath(blacklist_dest,mustWork=FALSE)) file.copy(f,blacklist_dest,overwrite=TRUE)
    break
  }
  if (!file.exists(blacklist_dest) || file.info(blacklist_dest)$size<5000) {
    archive <- "/mnt/data1/lin_env/R/epiAneufinder-main.zip"
    if(file.exists(archive)) {
      members <- utils::unzip(archive,list=TRUE)$Name
      x <- grep("/sample_data/hg38-blacklist\\.v2\\.bed$",members,value=TRUE)
      if(length(x)==1L) {
        utils::unzip(archive,files=x,exdir=TMP)
        f <- file.path(TMP,x)
        if(file.exists(f)) file.copy(f,blacklist_dest,overwrite=TRUE)
      }
    }
  }
  if (!file.exists(blacklist_dest) || file.info(blacklist_dest)$size<5000) {
    url <- "https://raw.githubusercontent.com/colomemaria/epiAneufinder/main/sample_data/hg38-blacklist.v2.bed"
    logmsg("Attempting to fetch official hg38 blacklist: ",url)
    try(suppressWarnings(utils::download.file(url,blacklist_dest,quiet=TRUE,mode="wb")),silent=TRUE)
  }
  if (!file.exists(blacklist_dest) || file.info(blacklist_dest)$size<5000)
    stop("hg38 blacklist unavailable. Download the official hg38-blacklist.v2.bed to ",blacklist_dest," and rerun; mtDNA output already saved.")
  bl <- fread(blacklist_dest,header=FALSE,nrows=10,fill=TRUE)
  if(ncol(bl)<3 || !all(grepl("^chr",bl[[1]]))) stop("Malformed hg38 blacklist: ",blacklist_dest)

  # ================== 4. Pre-filter original fragments BY SAMPLE ==================
  if (!requireNamespace("epiAneufinder",quietly=TRUE)) stop("epiAneufinder is not loadable")
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38",quietly=TRUE)) stop("hg38 genome package is not loadable")
  params <- list(windowSize=WINDOW_SIZE,minFrags=MIN_FRAGS,cores=NCORES,
                 k=4L,minsizeCNV=MIN_CONSECUTIVE_BINS,gc_correction=TRUE,
                 blacklist_md5=unname(tools::md5sum(blacklist_dest)))
  for (s in c("HCC9_C1","HCC9_C2")) {
    logmsg("Preparing filtered fragments and epiAneufinder for ",s)
    source_file <- file.path(FRAG_ROOT,paste0(s,"_fragments.tsv.gz"))
    check_exists(source_file)
    whitelist <- file.path(TMP,paste0(s,"_barcodes.txt"))
    writeLines(sort(unique(m[sample==s,barcode])),whitelist)
    filtered <- file.path(TMP,paste0(s,"_selected_fragments.tsv"))
    manifest_file <- file.path(TMP,paste0(s,"_filter_manifest.rds"))
    manifest <- list(whitelist_md5=unname(tools::md5sum(whitelist)),
                     input_size=as.numeric(file.info(source_file)$size),
                     input_mtime=as.character(file.info(source_file)$mtime))
    old <- if(file.exists(manifest_file)) tryCatch(readRDS(manifest_file),error=function(e) NULL) else NULL
    if(file.exists(filtered) && file.info(filtered)$size>0 && identical(old,manifest)) {
      logmsg("Reusing verified pre-filtered fragments: ",filtered)
    } else {
      epi_subdir <- file.path(OUT,paste0("epiAneufinder_",sub("HCC9_","",s)))
      if(file.exists(file.path(epi_subdir,"epiAneufinder_results","count_summary.rds")))
        stop("Cell whitelist/source changed but epiAneufinder cache exists in ",epi_subdir,
             ". Move/remove this specific cache before rerun (not done automatically).")
      target_tmp <- paste0(filtered,".incomplete")
      if(file.exists(target_tmp)) unlink(target_tmp)
      # Only matching barcodes, no header, fragment count column preserved.
      awk_code <- "NR==FNR{keep[$1]=1;next} ($4 in keep){print}"
      cmd <- paste("LC_ALL=C gzip -dc --",shQuote(source_file),"| LC_ALL=C awk -F '\\t'",
                   shQuote(awk_code),shQuote(whitelist),"- >",shQuote(target_tmp))
      status <- suppressWarnings(system2("bash",c("-o","pipefail","-c",shQuote(cmd)),stdout="",stderr=""))
      if(!identical(as.integer(status),0L) || !file.exists(target_tmp) || file.info(target_tmp)$size==0)
        stop("Fragment whitelist pre-filter failed for ",s,"; inspect file/disk space. status=",status)
      if(!file.rename(target_tmp,filtered)) stop("Cannot rename prefiltered fragment file: ",filtered)
      saveRDS(manifest,manifest_file)
      logmsg("Saved selected fragments, GiB: ",round(file.info(filtered)$size/1024^3,2))
    }
    epi_subdir <- file.path(OUT,paste0("epiAneufinder_",sub("HCC9_","",s)))
    dir.create(epi_subdir,recursive=TRUE,showWarnings=FALSE)
    params_path <- file.path(epi_subdir,"HCC9_epianeufinder_params.rds")
    p_this <- c(params,list(sample=s,whitelist_md5=manifest$whitelist_md5,
                            input_size=manifest$input_size,input_mtime=manifest$input_mtime))
    p_old <- if(file.exists(params_path)) readRDS(params_path) else NULL
    if (!is.null(p_old) && !identical(p_this,p_old)) stop("Existing epiAneufinder results have different parameters for ",s)
    saveRDS(p_this,params_path)
    result_path <- file.path(epi_subdir,"epiAneufinder_results","results_table.tsv")
    if (!file.exists(result_path) || file.info(result_path)$size==0) {
      logmsg("Running epiAneufinder on ",s," (GC correction + segmentation + Gain/Normal/Loss)...")
      epiAneufinder::epiAneufinder(input=filtered,outdir=epi_subdir,
         blacklist=blacklist_dest,windowSize=WINDOW_SIZE,
         genome="BSgenome.Hsapiens.UCSC.hg38",exclude=c("chrX","chrY","chrM"),
         reuse.existing=TRUE, selected_cells=whitelist, ncores=NCORES,
         minFrags=MIN_FRAGS,k=4L,minsize=1L,minsizeCNV=MIN_CONSECUTIVE_BINS,
         gc_correction=TRUE,plotKaryo=FALSE,save_removed_regions=TRUE,
         threshold_cells_nbins=0.05,threshold_blacklist_bins=0.85)
    }
    if(!file.exists(result_path) || file.info(result_path)$size==0) stop("epiAneufinder result missing for ",s)
    logmsg("Completed CNV inference: ",result_path)
    gc(FALSE)
  }

  # ================== 5. Parse official 0/1/2 per-cell outputs ==================
  read_epi_states <- function(s) {
    result_path <- file.path(OUT,paste0("epiAneufinder_",sub("HCC9_","",s)),
                             "epiAneufinder_results","results_table.tsv")
    d <- fread(result_path,check.names=FALSE,showProgress=FALSE)
    if(!all(c("seq","start","end") %in% names(d)))
      stop("Unexpected epiAneufinder results_table.tsv layout in ",result_path,
           "; columns: ",paste(head(names(d),12),collapse=","))
    cell_cols <- grep("^cell[-.]",names(d),value=TRUE)
    if(!length(cell_cols)) stop("No 'cell-' columns in ",result_path)
    bcs <- sub("^cell[-.]","",cell_cols)
    if (anyDuplicated(bcs)) stop("Duplicate barcodes in epiAneufinder output for ",s)
    info <- data.table(chr=as.character(d$seq),start=as.integer(d$start),end=as.integer(d$end))
    info[,key:=paste(chr,start,end,sep=":")]
    if(anyDuplicated(info$key)) stop("Duplicated genomic bins for ",s)
    if(anyNA(info$start) || anyNA(info$end)) stop("Missing coordinates in epiAneufinder result")
    mat <- as.matrix(d[,..cell_cols]); storage.mode(mat) <- "integer"
    if(any(!is.na(mat) & !(mat %in% 0:2))) stop("Unexpected states outside 0,1,2 for ",s)
    colnames(mat) <- bcs; rownames(mat) <- info$key
    expected <- m[sample==s,barcode]
    if(!all(colnames(mat)%in%expected)) stop("Unmatched epiAneufinder cells in ",s)
    qctab <- data.table(sample=s,Expected_cells=length(expected),Called_cells=ncol(mat),
               Missing_cells=length(setdiff(expected,colnames(mat))))
    fwrite(qctab,file.path(OUT,paste0(s,"_CNV_cell_QC.csv")))
    if(ncol(mat)<10) stop("Too few cells retained for sample ",s)
    list(info=info,mat=mat)
  }
  logmsg("Reading official epiAneufinder per-cell states...")
  a <- read_epi_states("HCC9_C1"); b <- read_epi_states("HCC9_C2")
  common <- intersect(a$info$key,b$info$key)
  if(length(common)<500L) stop("Too few shared valid CNV windows between C1/C2: ",length(common))
  info <- a$info[match(common,a$info$key)]
  chr_order <- paste0("chr",1:22)
  info[,chr_ord:=match(chr,chr_order)]
  setorder(info,chr_ord,start,end)
  info[,chr_ord:=NULL]
  common <- info$key
  sample_mats <- list(HCC9_C1=a$mat[common,,drop=FALSE],HCC9_C2=b$mat[common,,drop=FALSE])
  rm(a,b); gc(FALSE)
  for(s in names(sample_mats)) {
    kept <- colnames(sample_mats[[s]])
    logmsg(s," epiAneufinder retained ",length(kept)," / ",nrow(m[sample==s]),
           " mapped cells (",length(setdiff(m[sample==s,barcode],kept))," excluded)")
    for(k in clones) {
      x <- intersect(m[sample==s & Clone==k,barcode],kept)
      if(length(x)<10L) stop(s," Clone",k," has only ",length(x)," retained cells; cannot robustly compare. Inspect minFrags/coverage QC.")
    }
  }
  saveRDS(list(regions=info,samples=sample_mats,metadata=m,
               Note="epiAneufinder inferred relative states 0=Loss,1=Normal,2=Gain; no DNA validation"),
          file.path(OUT,"HCC9_epianeufinder_percell_states_v4.rds"),compress=TRUE)

  # ================== 6. CNV event markers matched to mtDNA-defined Clone ==================
  logmsg("Aggregating inferred CNV events within each mtDNA Clone...")
  # Find bins in which >=50% of the clone cells show a Gain or Loss in BOTH samples.
  # Group only consecutive 1Mb bins of the same state, require >=3 consecutive bins.
  all_rates <- list()
  for(s in names(sample_mats)) {
    mat <- sample_mats[[s]]
    rates <- list()
    for(k in clones) {
      cols <- intersect(m[sample==s & Clone==k,barcode],colnames(mat))
      rates[[paste0(k,"_0")]] <- rowMeans(mat[,cols,drop=FALSE]==0L,na.rm=TRUE)
      rates[[paste0(k,"_2")]] <- rowMeans(mat[,cols,drop=FALSE]==2L,na.rm=TRUE)
    }
    all_rates[[s]] <- rates
  }
  peer_ids <- function(k) setdiff(clones,k)
  rows <- list(); idx <- 0L
  v2 <- if(file.exists(V2_PATH)) fread(V2_PATH) else NULL
  for (k in clones) for (event in c(0L,2L)) {
    event_key <- paste0(k,"_",event)
    r1 <- all_rates$HCC9_C1[[event_key]]
    r2 <- all_rates$HCC9_C2[[event_key]]
    eligible <- which(is.finite(r1)&is.finite(r2)&r1>=MIN_CLONE_EVENT_RATE&r2>=MIN_CLONE_EVENT_RATE)
    if(!length(eligible)) next
    for(chrom in chr_order) {
      i_chr <- eligible[info$chr[eligible]==chrom]
      if(!length(i_chr)) next
      ord <- order(info$start[i_chr]); i_chr <- i_chr[ord]
      group <- cumsum(c(TRUE, diff(info$start[i_chr])!=WINDOW_SIZE))
      for (g in unique(group)) {
        bins <- i_chr[group==g]
        if(length(bins)<MIN_CONSECUTIVE_BINS) next
        spl <- list(); peer_rates <- list(); pairwise <- list()
        for(s in names(sample_mats)) {
          mm <- sample_mats[[s]]
          # A cell supports a region when >=70% of its valid bins have the same inferred state.
          cell_fraction <- colMeans(mm[bins,,drop=FALSE]==event,na.rm=TRUE)
          ix <- intersect(m[sample==s & Clone==k,barcode],names(cell_fraction))
          oo <- intersect(m[sample==s & Clone!=k,barcode],names(cell_fraction))
          x <- cell_fraction[ix]; y <- cell_fraction[oo]
          p <- vapply(peer_ids(k),function(o){
            ids <- intersect(m[sample==s & Clone==o,barcode],names(cell_fraction))
            mean(cell_fraction[ids]>=MIN_PER_CELL_BIN_FRACTION,na.rm=TRUE)
          },numeric(1))
          names(p) <- as.character(peer_ids(k))
          ax <- vapply(peer_ids(k),function(o){
            ids <- intersect(m[sample==s & Clone==o,barcode],names(cell_fraction))
            auroc(x,cell_fraction[ids])
          },numeric(1))
          spl[[s]] <- list(n=length(ix),support=mean(x>=MIN_PER_CELL_BIN_FRACTION,na.rm=TRUE),
             other=mean(y>=MIN_PER_CELL_BIN_FRACTION,na.rm=TRUE),auc=auroc(x,y))
          peer_rates[[s]] <- p; pairwise[[s]] <- ax
        }
        p1 <- spl$HCC9_C1; p2 <- spl$HCC9_C2
        peers <- peer_ids(k)
        shared <- peers[peer_rates$HCC9_C1>=0.50 & peer_rates$HCC9_C2>=0.50]
        min_peer_auc <- min(c(pairwise$HCC9_C1,pairwise$HCC9_C2),na.rm=TRUE)
        specific <- p1$support>=MIN_SPECIFIC_RATE && p2$support>=MIN_SPECIFIC_RATE &&
           p1$other<=MAX_OTHER_RATE && p2$other<=MAX_OTHER_RATE &&
           p1$support-p1$other>=MIN_DELTA_RATE && p2$support-p2$other>=MIN_DELTA_RATE &&
           min(p1$auc,p2$auc)>=MIN_AUROC && min_peer_auc>=MIN_AUROC &&
           length(shared)==0L
        classification <- if(specific) "Clone_specific_candidate" else if(length(shared)>=4) "Shared_CNV_candidate" else if(length(shared)>=1) "Subgroup_shared_CNV_candidate" else "Weak_or_unresolved_candidate"
        st <- min(info$start[bins]); en <- max(info$end[bins]); v2matches <- character()
        if(!is.null(v2) && all(c("Clone","Chromosome","Start","End")%in%names(v2))) {
          r <- v2[Clone==paste0("Clone",k)&Chromosome==chrom&as.numeric(Start)<=en&as.numeric(End)>=st]
          if(nrow(r)) v2matches <- paste0(r$Chromosome,":",format(r$Start,scientific=FALSE,trim=TRUE),"-",format(r$End,scientific=FALSE,trim=TRUE))
        }
        idx <- idx+1L
        rows[[idx]] <- data.table(Clone=paste0("Clone",k),Chromosome=chrom,
           Start=st,End=en,Region=paste0(chrom,":",st,"-",en),
           CNV_state=if(event==2L)"Putative_Gain" else "Putative_Loss",
           Bin_count=length(bins),
           C1_clone_cells=p1$n,C2_clone_cells=p2$n,
           C1_clone_support_pct=round(100*p1$support,1),C2_clone_support_pct=round(100*p2$support,1),
           C1_other_support_pct=round(100*p1$other,1),C2_other_support_pct=round(100*p2$other,1),
           C1_AUROC=round(p1$auc,3),C2_AUROC=round(p2$auc,3),
           Min_pairwise_AUROC=round(min_peer_auc,3),
           Shared_with_other_clones=paste0("Clone",shared,collapse=";"),
           N_peer_clones_same_state=length(shared),Marker_class=classification,
           Priority=if(specific)"priority_clone_specific" else "descriptive_or_weak",
           Overlap_with_v2=length(v2matches)>0,
           V2_overlapping_regions=paste(unique(v2matches),collapse=";"),
           Evidence="epiAneufinder per-cell inferred relative CNV (0/1/2); NOT validated absolute copy number",
           Rule=">=3 adjacent 1Mb bins; >=50% cells state-consistent in each sample; cell supports if >=70% bins match")
      }
    }
  }
  if (length(rows)) {
    cnv <- rbindlist(rows,fill=TRUE)
    cnv[,Clone_number:=as.integer(sub("Clone","",Clone,fixed=TRUE))]
    cnv[,priority_rank:=as.integer(Priority=="priority_clone_specific")]
    setorder(cnv,Clone_number,-priority_rank,-Min_pairwise_AUROC,Chromosome,Start)
    cnv[,Rank_in_clone:=seq_len(.N),by=Clone]
    cnv[,c("Clone_number","priority_rank"):=NULL]
  } else {
    cnv <- data.table(Clone=character(),Chromosome=character(),Start=integer(),End=integer(),
         Region=character(),CNV_state=character(),Marker_class=character(),Priority=character())
  }
  fwrite(cnv,file.path(OUT,"HCC9_CNV_marker_new.csv"))
  sm <- data.table(Clone=paste0("Clone",clones))
  sm[,Cells_C1:=vapply(clones,function(k)sum(m$sample=="HCC9_C1"&m$Clone==k),integer(1))]
  sm[,Cells_C2:=vapply(clones,function(k)sum(m$sample=="HCC9_C2"&m$Clone==k),integer(1))]
  sm[,Retained_C1:=vapply(clones,function(k)sum(colnames(sample_mats$HCC9_C1)%in%m[sample=="HCC9_C1"&Clone==k,barcode]),integer(1))]
  sm[,Retained_C2:=vapply(clones,function(k)sum(colnames(sample_mats$HCC9_C2)%in%m[sample=="HCC9_C2"&Clone==k,barcode]),integer(1))]
  sm[,CNV_events:=vapply(Clone,function(cl)sum(cnv$Clone==cl),integer(1))]
  sm[,Clone_specific_candidates:=vapply(Clone,function(cl)sum(cnv$Clone==cl&cnv$Priority=="priority_clone_specific"),integer(1))]
  sm[,mtDNA_candidates:=vapply(Clone,function(cl)sum(mt$Clone==cl),integer(1))]
  fwrite(sm,file.path(OUT,"HCC9_marker_summary_new.csv"))
  logmsg("HCC9 finished. CNV event rows: ",nrow(cnv),"; clone-specific candidates: ",sum(cnv$Priority=="priority_clone_specific"))
  print(sm)
  logmsg("Saved ",file.path(OUT,"HCC9_mtDNA_marker_new.csv"))
  logmsg("Saved ",file.path(OUT,"HCC9_CNV_marker_new.csv"))
  logmsg("Saved ",file.path(OUT,"HCC9_marker_summary_new.csv"))
  logmsg("All Gain/Loss states are epiAneufinder-inferred relative calls, not DNA sequencing validation.")


}

# ================== 7. Clone-level report (original 2.R) ==================
logmsg("Recalculating clone-level CNV burdens from saved 0/1/2 states...")
IN <- file.path(OUT, "HCC9_epianeufinder_percell_states_v4.rds")
stopifnot(file.exists(IN))
X <- readRDS(IN)
stopifnot(all(c("regions", "samples", "metadata") %in% names(X)))
meta <- as.data.table(X$metadata)
info <- as.data.table(X$regions)
need_meta <- c("sample", "Clone", "barcode")
need_info <- c("chr", "start", "end", "key")
stopifnot(all(need_meta %in% names(meta)), all(need_info %in% names(info)))
stopifnot(!anyDuplicated(info$key), !anyNA(info$key))

SS <- c("HCC9_C1", "HCC9_C2")
stopifnot(all(SS %in% names(X$samples)))
mats <- lapply(SS, function(s) {
  z <- X$samples[[s]]
  stopifnot(is.matrix(z), !is.null(rownames(z)), !is.null(colnames(z)),
            !anyDuplicated(rownames(z)), !anyDuplicated(colnames(z)),
            all(info$key %in% rownames(z)))
  z <- z[info$key, , drop=FALSE]
  if(any(!is.na(z) & !(z %in% 0:2))) stop("Unexpected CNV state outside 0/1/2 in ",s)
  z
})
names(mats) <- SS

CLONES <- 1:8
stopifnot(identical(sort(unique(as.integer(meta$Clone))),CLONES))
# Reuse the earlier threshold; a region requires >=3 adjacent 1Mb windows.
MIN_CLONE_RATE <- 0.50
MIN_CONSECUTIVE_BINS <- 3L

one_sample <- function(mat, ids) {
  ids <- intersect(ids, colnames(mat))
  if(length(ids)<10L) stop("Too few retained cells: ",length(ids))
  z <- mat[,ids,drop=FALSE]
  n <- rowSums(!is.na(z))
  getp <- function(state) {
    q <- rowSums(z==state,na.rm=TRUE)/n
    q[n==0] <- NA_real_
    q
  }
  # Mean over CELLS of per-cell fraction of analysable windows with this call.
  get_burden <- function(state) {
    vals <- colMeans(z==state,na.rm=TRUE)
    if(all(!is.finite(vals))) return(NA_real_)
    mean(vals[is.finite(vals)])
  }
  list(n=length(ids),loss=getp(0L),normal=getp(1L),gain=getp(2L),
       loss_burden=get_burden(0L),gain_burden=get_burden(2L))
}

window_tables <- vector("list",length(CLONES))
burden_tables <- vector("list",length(CLONES))
for(k in CLONES) {
  smp <- lapply(SS,function(s) {
    ids <- meta[sample==s & Clone==k,barcode]
    one_sample(mats[[s]],ids)
  })
  names(smp) <- SS
  c1 <- smp$HCC9_C1; c2 <- smp$HCC9_C2
  gg <- is.finite(c1$gain)&is.finite(c2$gain)&
        c1$gain>=MIN_CLONE_RATE&c2$gain>=MIN_CLONE_RATE
  ll <- is.finite(c1$loss)&is.finite(c2$loss)&
        c1$loss>=MIN_CLONE_RATE&c2$loss>=MIN_CLONE_RATE
  nn <- is.finite(c1$normal)&is.finite(c2$normal)&
        c1$normal>=MIN_CLONE_RATE&c2$normal>=MIN_CLONE_RATE
  # A tied Gain/Loss call is inconclusive, never force it into either state.
  state <- rep("Uncertain",nrow(info))
  state[gg & !ll] <- "Gain"
  state[ll & !gg] <- "Loss"
  state[nn & !gg & !ll] <- "Normal"
  window_tables[[k]] <- data.table(
    Clone=paste0("Clone",k),Chromosome=info$chr,Start=info$start,End=info$end,
    State=state,
    C1_Gain_pct=round(100*c1$gain,1),C2_Gain_pct=round(100*c2$gain,1),
    C1_Loss_pct=round(100*c1$loss,1),C2_Loss_pct=round(100*c2$loss,1),
    C1_Normal_pct=round(100*c1$normal,1),C2_Normal_pct=round(100*c2$normal,1))
  # Pool cells from C1 and C2 with equal weight PER CELL (not per sample).
  g <- weighted.mean(c(c1$gain_burden,c2$gain_burden),w=c(c1$n,c2$n),na.rm=TRUE)
  l <- weighted.mean(c(c1$loss_burden,c2$loss_burden),w=c(c1$n,c2$n),na.rm=TRUE)
  dominant <- if(!is.finite(g+l) || g+l==0) "No_Gain_or_Loss" else if(g>l) "Gain_dominant" else if(l>g) "Loss_dominant" else "Balanced"
  burden_tables[[k]] <- data.table(
    Clone=paste0("Clone",k),Cells_C1=c1$n,Cells_C2=c2$n,
    C1_Gain_burden_pct=round(100*c1$gain_burden,2),
    C2_Gain_burden_pct=round(100*c2$gain_burden,2),
    C1_Loss_burden_pct=round(100*c1$loss_burden,2),
    C2_Loss_burden_pct=round(100*c2$loss_burden,2),
    Gain_burden_pct=round(100*g,2),Loss_burden_pct=round(100*l,2),
    Gain_of_abnormal_pct=if(g+l>0) round(100*g/(g+l),2) else NA_real_,
    Overall_abnormal_tendency=dominant)
}
W <- rbindlist(window_tables)
B <- rbindlist(burden_tables)

# Merge adjacent same-state windows WITHIN each Clone, chromosome and state.
# A skipped/normal/uncertain window breaks the run; avoid counting disjoint bins
# as a single CNV event. Keep intervals for interpretation and gene annotation.
seg <- W[State %in% c("Gain","Loss"),.(Clone,Chromosome,Start,End,State)]
setorder(seg,Clone,Chromosome,State,Start)
seg[,PrevEnd:=shift(End),by=.(Clone,Chromosome,State)]
seg[,Run:=cumsum(is.na(PrevEnd) | (Start != PrevEnd + 1)),by=.(Clone,Chromosome,State)]
E <- seg[,.(Start=min(Start),End=max(End),Window_count=.N),
         by=.(Clone,Chromosome,State,Run)][Window_count>=MIN_CONSECUTIVE_BINS]
E[,Run:=NULL]
E[,`:=`(Length_Mb=round((End-Start+1)/1e6,2),
         Region=paste0(Chromosome,":",Start,"-",End))]
setcolorder(E,c("Clone","Chromosome","Start","End","Region","State","Window_count","Length_Mb"))

# For each candidate interval, compute the fraction of clone cells in C1 / C2
# that are in the selected state across >=70% of its non-missing windows.
if(nrow(E)) {
  rates <- lapply(seq_len(nrow(E)),function(i) {
    r <- E[i]
    idx <- which(info$chr==r$Chromosome & info$start>=r$Start & info$end<=r$End)
    code <- if(r$State=="Gain") 2L else 0L
    vapply(SS,function(s) {
      mm <- mats[[s]]
      k <- as.integer(sub("Clone","",r$Clone,fixed=TRUE))
      ids <- intersect(meta[sample==s & Clone==k,barcode],colnames(mm))
      frac <- colMeans(mm[idx,ids,drop=FALSE]==code,na.rm=TRUE)
      if(all(!is.finite(frac))) return(NA_real_)
      round(100*mean(frac[is.finite(frac)]>=0.70),1)
    },numeric(1))
  })
  pr <- do.call(rbind,rates)
  E[,`:=`(C1_cell_support_pct=pr[,1],C2_cell_support_pct=pr[,2])]
}

summary <- copy(B)
summary[,`:=`(Gain_regions=NA_character_,Loss_regions=NA_character_,
              Gain_region_count=0L,Loss_region_count=0L,
              Gain_region_Mb=0,Loss_region_Mb=0,
              Gain_consensus_window_count=0L,Loss_consensus_window_count=0L)]
for(i in seq_len(nrow(summary))) {
  cl <- summary$Clone[i]
  for(st in c("Gain","Loss")) {
    sub <- E[Clone==cl & State==st]
    prefix <- st
    summary[i,(paste0(prefix,"_region_count")):=nrow(sub)]
    summary[i,(paste0(prefix,"_region_Mb")):=round(sum(sub$Length_Mb),2)]
    summary[i,(paste0(prefix,"_regions")):=if(nrow(sub)) paste(sub$Region,collapse="; ") else "None"]
    summary[i,(paste0(prefix,"_consensus_window_count")):=W[Clone==cl & State==st,.N]]
  }
}
summary[,Note:="Gain/Loss dominant = greater CELL-weighted mean per-cell window fraction; intervals require >=50% cells in BOTH samples and >=3 adjacent 1Mb windows; inferred scATAC states only"]

f1 <- file.path(OUT,"HCC9_CNV_clone_level_summary.csv")
f2 <- file.path(OUT,"HCC9_CNV_clone_level_regions.csv")
f3 <- file.path(OUT,"HCC9_CNV_clone_level_windows.csv")
fwrite(summary,f1,bom=TRUE)
fwrite(E,f2,bom=TRUE)
fwrite(W,f3,bom=TRUE)
cat("\n===== Clone-level result from ORIGINAL per-cell 0/1/2 states =====\n")
print(summary[,.(Clone,Cells_C1,Cells_C2,Gain_burden_pct,Loss_burden_pct,
                Overall_abnormal_tendency,Gain_region_count,Loss_region_count)])
cat("\nWrote:\n",f1,"\n",f2,"\n",f3,"\n")

# Optional 3-sheet workbook if a spreadsheet writer is already installed.
mt_file <- file.path(OUT,"HCC9_mtDNA_marker_new.csv")
if (!file.exists(mt_file)) stop("Missing mtDNA marker CSV: ", mt_file, "; run with HCC9_MODE=full first")
{
  mt <- fread(mt_file)
  mtcols <- intersect(c("Clone","Marker","Positive_cells","Clone_cells",
    "Clone_positive_pct","Other_positive_pct","Enrichment_pp",
    "C1_positive_pct","C2_positive_pct","Mean_positive_VAF_pct",
    "C1_VAF_AUROC","C2_VAF_AUROC","Priority"), names(mt))
  mt <- mt[,..mtcols]
  # New overview uses recalculated CNV quantities, not old 66-marker counts.
  sm <- summary[,.(Clone,Cells_C1,Cells_C2,Total_cells=Cells_C1+Cells_C2,
                    Gain_region_count,Loss_region_count,
                    Gain_burden_pct,Loss_burden_pct,Overall_abnormal_tendency)]
  mcount <- mt[,.(mtDNA_candidates=.N),by=Clone]
  sm <- merge(sm,mcount,by="Clone",all.x=TRUE,sort=FALSE)
  sm[is.na(mtDNA_candidates),mtDNA_candidates:=0L]
  tabs <- list("01_总览"=sm,"02_CNV_按Clone汇总"=summary,"03_mtDNA_Marker"=mt)
  outfile <- file.path(OUT,"HCC9_Markers_3Sheets_CloneRecalculated.xlsx")
  if(requireNamespace("openxlsx",quietly=TRUE)) {
    wb <- openxlsx::createWorkbook()
    for (nm in names(tabs)) {
      openxlsx::addWorksheet(wb,nm)
      openxlsx::writeData(wb,nm,tabs[[nm]],withFilter=TRUE)
      openxlsx::freezePane(wb,nm,firstRow=TRUE)
      openxlsx::addStyle(wb,nm,openxlsx::createStyle(fgFill="#E7EEF6",textDecoration="bold"),
                         rows=1,cols=seq_len(ncol(tabs[[nm]])),gridExpand=TRUE)
      openxlsx::setColWidths(wb,nm,cols=seq_len(ncol(tabs[[nm]])),width=18)
      if(nm=="02_CNV_按Clone汇总") {
        openxlsx::setColWidths(wb,nm,cols=which(names(tabs[[nm]]) %in% c("Gain_regions","Loss_regions","Note")),width=45)
      }
    }
    openxlsx::saveWorkbook(wb,outfile,overwrite=TRUE)
  } else {
    writexl::write_xlsx(tabs,outfile)
  }
  cat("3-sheet Excel:\n",outfile,"\n")
}


logmsg("Final Excel: ", file.path(OUT, "HCC9_Markers_3Sheets_CloneRecalculated.xlsx"))
