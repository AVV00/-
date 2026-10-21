#!/usr/bin/env Rscript
# OXPHOS-only trajectory: all matched cells, nuclear genes only.
# Input is continuous gene activity, NOT raw counts; no additional log/size normalization.
# Gaussian CDS; no differentialGeneTest, cluster selection, or mitochondrial matrix.
# All cells are plotted. Root is automatic unless ROOT_STATE is explicitly supplied.
# Checkpoints are stage-level, NOT per-iteration restart states.
options(stringsAsFactors=FALSE, warn=1)
project_tmp <- "/mnt/data/linyanlin/Fig5/2/tmp"
dir.create(project_tmp, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(TMPDIR=project_tmp, TMP=project_tmp, TEMP=project_tmp,
           SQLITE_TMPDIR=project_tmp, OMP_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(monocle)
  library(igraph)
  library(ggplot2)
})
data.table::setDTthreads(8L)
set.seed(1234)
cat("Started:",format(Sys.time()),"\nR:",R.version.string,
    "\nMonocle:",as.character(packageVersion("monocle")),
    "\nDDRTree:",as.character(packageVersion("DDRTree")),
    "\nigraph:",as.character(packageVersion("igraph")),"\n")
ns <- asNamespace("monocle")
patch_one <- function(name, transform, helpers=list(), required=TRUE) {
  old <- get(name, envir=ns, inherits=FALSE)
  txt <- paste(deparse(old), collapse="\n")
  changed <- transform(txt)
  if (identical(txt, changed)) {
    if (required) stop("No expected source change: ", name)
    return(invisible(NULL))
  }
  env <- new.env(parent=environment(old))
  list2env(helpers, envir=env)
  fn <- eval(parse(text=changed)[[1]], envir=env)
  locked <- bindingIsLocked(name, ns)
  if (locked) unlockBinding(name, ns)
  tryCatch(assign(name, fn, envir=ns),
           finally={if (locked) lockBinding(name, ns)})
  cat("Patched:", name, "\n")
}
patch_one("project2MST", function(txt) {
  pattern <- paste0("V\\(dp_mst\\)\\[suppressWarnings\\(nei\\(",
                    "closest_vertex_names\\[i\\],[[:space:]]*",
                    "mode[[:space:]]*=[[:space:]]*\"all\"\\)\\)\\]")
  hits <- gregexpr(pattern, txt, perl=TRUE)[[1]]
  if (sum(hits > 0L) != 1L) stop("Unexpected project2MST source; no replacement made.")
  sub(pattern, 'igraph::neighbors(dp_mst, closest_vertex_names[i], mode="all")',
      txt, perl=TRUE)
})
patch_one("extract_ddrtree_ordering", function(txt) {
  txt <- gsub("graph\\.dfs[[:space:]]*\\(", "igraph::dfs(", txt, perl=TRUE)
  txt <- gsub("\\bneimode[[:space:]]*=", "mode =", txt, perl=TRUE)
  if ("parent" %in% names(formals(igraph::dfs))) {
    txt <- gsub("\\bfather[[:space:]]*=", "parent =", txt, perl=TRUE)
    txt <- gsub("$father", "$parent", txt, fixed=TRUE)
  }
  txt
})

# Direct plotting avoids the removed dplyr select_()/group_by_() APIs.
plot_cell_trajectory <- function(cds, color_by="Pseudotime") {
  md <- pData(cds)
  stopifnot(color_by %in% names(md))
  z <- reducedDimS(cds)
  y <- reducedDimK(cds)
  g <- minSpanningTree(cds)
  stopifnot(!is.null(colnames(z)), !is.null(colnames(y)),
            !is.null(igraph::V(g)$name))
  ci <- match(rownames(md), colnames(z))
  vi <- match(igraph::V(g)$name, colnames(y))
  stopifnot(!anyNA(ci), !anyNA(vi))
  d <- data.frame(DDRTree1=as.numeric(z[1,ci]), DDRTree2=as.numeric(z[2,ci]),
                  colour_value=md[[color_by]])
  e <- igraph::as_edgelist(g, names=FALSE)
  seg <- data.frame(x=y[1,vi[e[,1]]], y=y[2,vi[e[,1]]],
                    xend=y[1,vi[e[,2]]], yend=y[2,vi[e[,2]]])
  p <- ggplot(d, aes(DDRTree1, DDRTree2)) +
    geom_segment(data=seg, aes(x=x,y=y,xend=xend,yend=yend),
                 inherit.aes=FALSE, colour="grey35", linewidth=0.5) +
    geom_point(aes(colour=colour_value), size=0.5, alpha=0.75) +
    coord_equal() + theme_classic(base_size=12) + labs(colour=color_by)
  if (is.numeric(md[[color_by]])) p <- p + scale_colour_viridis_c()
  p
}


main <- function() {
  outdir <- "/mnt/data1/linyanlin/scATAC/picture/pic5/oxygen"
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  cache_dir <- file.path(outdir,"checkpoints_OXPHOS_nuclear_v1")
  dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
  inputs <- c(
    cells="/mnt/data1/linyanlin/scATAC/picture/pic1/heatmap/hall/1019_new/final2.csv",
    matrix="/mnt/data1/linyanlin/scATAC/all_data/step4_geneannotation/Malignant_Cells_gene_activity_gene_by_cell_sparse.rds",
    pathway="/mnt/data1/linyanlin/scATAC/all_data/step5_pathway/pathway.csv"
  )
  stopifnot(all(file.exists(inputs)))
  # Do not silently reuse caches after changing features, parameters or package versions.
  info <- file.info(inputs)
  signature <- list(pipeline="OXPHOS_nuclear_v1", inputs=normalizePath(inputs),
    sizes=info$size, mtimes=as.numeric(info$mtime),
    expected_cells=61145L, expected_genes=101L, norm_method="none",
    pseudo_expr=0, relative_expr=FALSE, scaling=TRUE,
    maxIter=20L, tol=0.01, seed=1234L,
    versions=vapply(c("monocle","DDRTree","igraph","Matrix","VGAM"),
                    function(x) as.character(packageVersion(x)),character(1)))
  root_env <- Sys.getenv("ROOT_STATE","")
  root_state <- if (nzchar(root_env)) suppressWarnings(as.integer(root_env)) else NULL
  if (!is.null(root_state) && (is.na(root_state)||root_state<1L))
    stop("ROOT_STATE must be a positive integer or left unset.")
  atomic_rds <- function(value, dest) {
    tmp <- paste0(dest,".writing.",Sys.getpid())
    on.exit(unlink(tmp),add=TRUE)
    saveRDS(value,tmp,compress=FALSE)
    if (!file.rename(tmp,dest)) stop("Could not finalize checkpoint: ",dest)
  }
  stage <- function(name, compute, key=signature) {
    dest <- file.path(cache_dir,paste0(name,".rds"))
    if (file.exists(dest)) {
      cat("Read checkpoint:",dest,"\n")
      x <- readRDS(dest)
      if (!identical(x$signature,key))
        stop("Checkpoint settings differ. Keep old results and use a new output/cache directory: ",dest)
      return(x$result)
    }
    t0 <- proc.time()[["elapsed"]]
    cat("Compute:",name,"at",format(Sys.time()),"\n")
    value <- compute()
    atomic_rds(list(signature=key,result=value),dest)
    cat("Saved checkpoint:",name,"; seconds:",proc.time()[["elapsed"]]-t0,"\n")
    value
  }
  fwrite(data.table(input=names(inputs),path=unname(inputs),bytes=info$size,
                    modified=as.character(info$mtime)),file.path(outdir,"input_files.csv"))
  writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
  fwrite(data.table(parameter=c("genes","cells","family","normalization",
      "scaling","maxIter","tol","seed","root","threads"),
    value=c("KEGG_OXIDATIVE_PHOSPHORYLATION; nuclear only","all matched; no sampling",
      "VGAM::uninormal","none; input values unchanged","per-gene z-score inside reduceDimension",
      "20","0.01","1234",if(is.null(root_state)) "automatic" else as.character(root_state),
      "data.table=8; BLAS/OpenMP=1; DDRTree is not an 8-worker differential test")),
    file.path(outdir,"run_parameters.csv"))

  cat("\nSTEP 1: Prepare OXPHOS matrix\n")
  cds <- stage("01_OXPHOS_prepared",function() {
    f <- fread(inputs[["cells"]])
    stopifnot(all(c("cell","patient","MSS")%in%names(f)))
    f[,cell:=as.character(cell)]
    f[,MSS:=suppressWarnings(as.numeric(MSS))]
    stopifnot(!anyNA(f$cell),!anyDuplicated(f$cell))
    p <- fread(inputs[["pathway"]])
    stopifnot(all(c("gmt_geneset_name","genes")%in%names(p)))
    s <- p[gmt_geneset_name=="KEGG_OXIDATIVE_PHOSPHORYLATION",genes]
    s <- s[!is.na(s)]
    if (!length(s)) stop("OXPHOS pathway missing.")
    genes <- unique(trimws(unlist(strsplit(s,";",fixed=TRUE))))
    genes <- genes[nzchar(genes)]
    mt_alias <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
                  "CO1","CO2","CO3","CYB","ATP6","ATP8")
    nuclear <- genes[!startsWith(genes,"MT-") & !genes%in%mt_alias]
    m <- readRDS(inputs[["matrix"]])
    stopifnot(inherits(m,"Matrix"),!is.null(rownames(m)),!is.null(colnames(m)),
              !anyDuplicated(rownames(m)),!anyDuplicated(colnames(m)))
    report <- data.table(cell=f$cell,in_nuclear_matrix=f$cell%in%colnames(m),
                         finite_MSS=is.finite(f$MSS))
    report[,included:=in_nuclear_matrix & finite_MSS]
    fwrite(report,file.path(outdir,"cell_matching.csv"))
    selected <- f[report$included]
    if (nrow(selected)!=signature$expected_cells)
      stop("Expected 61145 cells; found ",nrow(selected),". Inspect cell_matching.csv.")
    matched <- intersect(nuclear,rownames(m))
    if (length(matched)!=signature$expected_genes)
      stop("Expected 101 matched nuclear genes; found ",length(matched))
    x <- as.matrix(m[matched,match(selected$cell,colnames(m)),drop=FALSE])
    rm(m); invisible(gc())
    stopifnot(identical(colnames(x),selected$cell),all(is.finite(x)),all(x>=0))
    st <- data.table(gene=matched,mean_activity=rowMeans(x),
                     variance_activity=apply(x,1,var),detected_cells=rowSums(x>0))
    st[,used_for_ordering:=is.finite(variance_activity)&variance_activity>1e-12]
    fwrite(st,file.path(outdir,"OXPHOS_gene_statistics.csv"))
    gr <- data.table(gene=genes)
    gr[,nuclear:=gene%in%nuclear]
    gr[,in_matrix:=gene%in%matched]
    gr[,used_for_ordering:=gene%in%st[used_for_ordering==TRUE,gene]]
    fwrite(gr,file.path(outdir,"KEGG_OXPHOS_gene_matching.csv"))
    if (sum(st$used_for_ordering)<3L) stop("Too few variable OXPHOS genes.")
    # Raw arithmetic mean on the same 101 nuclear genes used in the rank check.
    meta <- data.frame(cell=selected$cell,patient=as.character(selected$patient),
      MSS=selected$MSS,OXPHOS_nuclear_mean=colMeans(x),row.names=selected$cell)
    k <- "KEGG_OXIDATIVE_PHOSPHORYLATION"
    if (k%in%names(selected)) meta$OXPHOS_existing_score <- selected[[k]]
    fwrite(as.data.table(meta),file.path(outdir,"input_cell_metadata.csv"))
    fwrite(st[used_for_ordering==TRUE,.(gene)],file.path(outdir,"ordering_genes.csv"))
    # Retain all 101 matched genes in CDS; ordering filter selects only variable genes.
    x <- methods::as(Matrix::Matrix(x,sparse=TRUE),"dgCMatrix")
    cds0 <- newCellDataSet(x,phenoData=new("AnnotatedDataFrame",data=meta),
      featureData=new("AnnotatedDataFrame",data=data.frame(gene_short_name=rownames(x),
        row.names=rownames(x))),expressionFamily=VGAM::uninormal(),lowerDetectionLimit=0)
    sizeFactors(cds0) <- rep(1,ncol(cds0))
    cds0 <- setOrderingFilter(cds0,st[used_for_ordering==TRUE,gene])
    cat("Cells:",ncol(cds0),"; matched nuclear genes:",nrow(cds0),
        "; variable ordering genes:",sum(st$used_for_ordering),"\n")
    cds0
  })

  cat("\nSTEP 2: DDRTree; maxIter=20; tol=0.01\n")
  cat("One dense cell-pair matrix requires about",
      round(ncol(cds)^2*8/1024^3,2),"GiB; orderCells can allocate several large objects.\n")
  cds <- stage("02_DDRTree_pre_orderCells",function() {
    set.seed(signature$seed)
    reduceDimension(cds,max_components=2,reduction_method="DDRTree",
      norm_method="none",pseudo_expr=0,relative_expr=FALSE,scaling=TRUE,
      maxIter=signature$maxIter,tol=signature$tol,verbose=TRUE)
  })

  objectives <- as.numeric(unlist(cds@auxOrderingData[["DDRTree"]]$objective_vals))
  if (length(objectives)) {
    delta <- c(NA_real_,abs(diff(objectives))/abs(head(objectives,-1L)))
    fwrite(data.table(record=seq_along(objectives),objective=objectives,
                      relative_change=delta,tol=signature$tol),
           file.path(outdir,"DDRTree_objective_history.csv"))
    cat("Final recorded relative objective change:",tail(delta,1),"; tol:",signature$tol,"\n")
  }

  cat("\nSTEP 3: Order cells\n")
  cds <- stage("03_ordered_auto",function() orderCells(cds))
  if (!is.null(root_state)) {
    if (!as.character(root_state)%in%as.character(unique(pData(cds)$State)))
      stop("Requested ROOT_STATE is not among the automatic State labels.")
    auto_cds <- cds
    cds <- stage(paste0("04_ordered_root_",root_state),
      function() orderCells(auto_cds,root_state=root_state),
      key=list(base=signature,root_state=root_state))
  }
  stopifnot(all(is.finite(pData(cds)$Pseudotime)))
  # Standalone final CDS, distinct from wrapper checkpoints.
  atomic_rds(cds,file.path(outdir,"OXPHOS_nuclear_ordered_cds.rds"))

  cat("\nSTEP 4: Export CSV\n")
  md <- pData(cds)
  z <- reducedDimS(cds)
  ci <- match(rownames(md),colnames(z))
  stopifnot(!anyNA(ci))
  tab <- as.data.table(md)
  tab[,DDRTree1:=as.numeric(z[1,ci])]
  tab[,DDRTree2:=as.numeric(z[2,ci])]
  fwrite(tab,file.path(outdir,"OXPHOS_nuclear_ordered_cells.csv"))
  fwrite(tab[,.(cells=.N),by=State],file.path(outdir,"state_counts.csv"))
  fwrite(tab[,.(cells=.N),by=.(patient,State)],file.path(outdir,"patient_state_counts.csv"))
  g <- minSpanningTree(cds); y <- reducedDimK(cds)
  vi <- match(igraph::V(g)$name,colnames(y))
  stopifnot(!anyNA(vi))
  fwrite(data.table(vertex=igraph::V(g)$name,DDRTree1=as.numeric(y[1,vi]),
                    DDRTree2=as.numeric(y[2,vi])),file.path(outdir,"tree_nodes.csv"))
  fwrite(as.data.table(igraph::as_data_frame(g,what="edges")),file.path(outdir,"tree_edges.csv"))

  cat("\nSTEP 5: Export PDFs; all cells plotted\n")
  for (colour in c("Pseudotime","MSS","OXPHOS_nuclear_mean","State","patient")) {
    p <- plot_cell_trajectory(cds,color_by=colour)
    pdf(file.path(outdir,paste0("OXPHOS_nuclear_trajectory_",colour,".pdf")),
        width=8,height=7,useDingbats=FALSE)
    tryCatch(print(p),finally=dev.off())
  }
  print(table(pData(cds)$State))
  print(summary(pData(cds)$Pseudotime))
  cat("\nCOMPLETE:",format(Sys.time()),"\nOutput:",outdir,
      "\nCells:",ncol(cds),"\nOrdering genes:",sum(fData(cds)$use_for_ordering),
      "\nRoot:",if(is.null(root_state)) "automatic; direction requires biological interpretation" else root_state,
      "\nOXPHOS changes along this OXPHOS-defined trajectory are not independent validation.\n")
}
main()

