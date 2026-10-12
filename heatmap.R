suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(UCell)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ragg)
})


##############################
# paths
##############################

outdir <- "/mnt/data1/linyanlin/scATAC/picture/pic1/heatmap/hall/group"

annotation_file <- "/mnt/data1/linyanlin/scATAC/all_data/new_celltype-261009/all_cells_new_celltype_annotation_metadata.csv"

gene_activity_file <- "/mnt/data1/linyanlin/scATAC/all_data/new_celltype-261009/all_cells_gene_activity_new_celltype_annotation.rds"

mss_file <- "/mnt/data/linyanlin/Fig5/2/03_malignant_mutation_OXPHOS_MSS_matched_all_mtDNA_genes.csv.gz"

pathway_file <- "/mnt/data1/linyanlin/scATAC/all_data/step5_pathway/pathway.csv"


dir.create(
  outdir,
  recursive=TRUE,
  showWarnings=FALSE
)



##############################
# 1 malignant
##############################

ann <- fread(
  annotation_file,
  select=c("cell","manual_celltype")
)


ann <- ann[
  manual_celltype %chin%
    c(
      "Malignant_Hepatocyte_like",
      "Malignant_Cholangiocyte_like"
    ),
  .(cell)
]


ann <- unique(
  ann,
  by="cell"
)


cat(
"Malignant:",
nrow(ann),
"\n"
)



##############################
# 2 MSS
##############################

mss <- fread(
  cmd=paste(
    "gzip -dc",
    shQuote(mss_file)
  ),
  select=c("cell","MSS")
)


mss <- mss[
  is.finite(MSS)
]


mss <- mss[
  ,
  .(MSS=first(MSS)),
  by=cell
]


meta <- merge(
  ann,
  mss,
  by="cell"
)


cat(
"MSS malignant:",
nrow(meta),
"\n"
)



##############################
# 3 gene activity
##############################

obj <- readRDS(
  gene_activity_file
)


ga <- GetAssayData(
  obj,
  assay="gene.activity"
)


cat(
"GA:",
paste(dim(ga),collapse=" x "),
"\n"
)


cells_use <- meta$cell[
  meta$cell %in% colnames(ga)
]


meta <- meta[
  match(cells_use,cell)
]


ga_sub <- ga[
  ,
  cells_use,
  drop=FALSE
]


cat(
"Matched:",
length(cells_use),
"\n"
)



##############################
# 4 KEGG sets
##############################

pw <- fread(
  pathway_file
)


pw <- pw[
  !is.na(genes)
]


sets <- lapply(
  pw$genes,
  function(x){

    g <- unique(
      trimws(
        unlist(
          strsplit(
            x,
            ";",
            fixed=TRUE
          )
        )
      )
    )

    intersect(
      g,
      rownames(ga_sub)
    )
  }
)


names(sets) <- pw$gmt_geneset_name


sets <- sets[
  lengths(sets)>=5
]


cat(
"KEGG used:",
length(sets),
"\n"
)



##############################
# 5 UCell
##############################

cat("UCell...\n")


scores <- ScoreSignatures_UCell(
  matrix=ga_sub,
  features=sets
)


scores <- as.data.table(scores)

scores[,cell:=cells_use]


score_cols <- grep(
  "_UCell$",
  names(scores),
  value=TRUE
)



res <- merge(
  meta,
  scores,
  by="cell"
)


setorder(
  res,
  MSS,
  cell
)


cat(
"Final cells:",
nrow(res),
"\n"
)



##############################
# pathway name
##############################

pretty_pathway <- function(x){

  x <- gsub(
    "_UCell$",
    "",
    x
  )

  x <- gsub(
    "^KEGG_",
    "",
    x
  )

  x <- gsub(
    "_",
    " ",
    x
  )

  x
}



##############################
# heatmap function
##############################

plot_heatmap <- function(n_group){


cat(
"\nGROUP:",
n_group,
"\n"
)


x <- copy(res)

setorder(
  x,
  MSS,
  cell
)


x[
  ,
  group_id :=
    floor(
      (seq_len(.N)-1)
      *
      n_group/.N
    )+1
]


long <- melt(
  x,
  id.vars=c(
    "cell",
    "MSS",
    "group_id"
  ),
  measure.vars=score_cols,
  variable.name="pathway",
  value.name="UCell_score"
)



grp <- long[
  ,
  .(
    n_cells=uniqueN(cell),
    MSS_mean=mean(MSS),
    MSS_min=min(MSS),
    MSS_max=max(MSS),
    UCell_mean=mean(UCell_score)
  ),
  by=.(group_id,pathway)
]



grp[
  ,
  UCell_z :=
    {
      s <- sd(UCell_mean)

      if(!is.finite(s)||s==0){
        0
      }else{
        (UCell_mean-mean(UCell_mean))/s
      }
    },
  by=pathway
]



cor_dt <- long[
  ,
  {
    ct <- suppressWarnings(
      cor.test(
        MSS,
        UCell_score,
        method="spearman",
        exact=FALSE
      )
    )

    .(
      rho=unname(ct$estimate),
      p_value=ct$p.value
    )
  },
  by=pathway
]


cor_dt[
  ,
  p_adj_BH:=p.adjust(
    p_value,
    method="BH"
  )
]


grp <- merge(
  grp,
  cor_dt,
  by="pathway"
)


grp[
  ,
  pathway_label:=pretty_pathway(pathway)
]


##############################
# CSV
##############################

csv <- file.path(
  outdir,
  paste0(
    "KEGG61_MSSgroup_",
    n_group,
    ".csv"
  )
)


fwrite(
  grp,
  csv
)



##############################
# matrix
##############################

order_path <- cor_dt[
  order(rho),
  pathway
]


wide <- dcast(
  grp,
  pathway~group_id,
  value.var="UCell_z"
)


wide <- wide[
  match(
    order_path,
    pathway
  )
]


mat <- as.matrix(
  wide[,-1]
)


rownames(mat) <-
  pretty_pathway(order_path)


mat[!is.finite(mat)] <- 0



mss_bar <- grp[
  ,
  .(MSS_mean=unique(MSS_mean)),
  by=group_id
][order(group_id)]$MSS_mean



rho <- cor_dt[
  match(order_path,pathway),
  rho
]


rho_mat <- matrix(
  rho,
  ncol=1
)


rownames(rho_mat)<-
  rownames(mat)



##############################
# colors
##############################

col_main <- colorRamp2(
  c(-2,0,2),
  c(
    "#6C7FB8",
    "#F3EEE8",
    "#F04A3A"
  )
)


col_mss <- colorRamp2(
  c(
    min(mss_bar),
    median(mss_bar),
    max(mss_bar)
  ),
  c(
    "#2C7FB8",
    "#D9D9D9",
    "#F03B20"
  )
)



col_rho <- colorRamp2(
  c(
    min(rho),
    0,
    max(rho)
  ),
  c(
    "#6C7FB8",
    "#F3EEE8",
    "#F04A3A"
  )
)



##############################
# heatmaps
##############################

ha <- HeatmapAnnotation(
  MSS=mss_bar,
  col=list(
    MSS=col_mss
  ),
  show_annotation_name=FALSE
)



ht1 <- Heatmap(
  mat,
  name="UCell\nZ-score",
  col=col_main,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_column_names=FALSE,
  show_row_names=TRUE,
  row_names_side="left",
  row_names_gp=gpar(fontsize=8),
  bottom_annotation=ha,
  border=TRUE,
  column_title=paste0(
    "KEGG pathway activity across MSS-ordered malignant cells (",
    n_group,
    " groups; n = ",
    nrow(res),
    ")"
  )
)



ht2 <- Heatmap(
  rho_mat,
  name="MSS\ncorrelation",
  col=col_rho,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  show_column_names=FALSE,
  width=unit(8,"mm")
)



##############################
# PDF output
##############################

pdf_file <- file.path(
  outdir,
  paste0(
    "KEGG61_MSSgroup_",
    n_group,
    ".pdf"
  )
)


pdf(
  pdf_file,
  width=15,
  height=11,
  useDingbats=FALSE
)


draw(
  ht1+ht2,
  heatmap_legend_side="right",
  annotation_legend_side="right"
)


dev.off()



cat(
"SAVED:",
pdf_file,
"\n"
)

}



##############################
# run
##############################

for(
  ng in c(1000,200,100)
){

  plot_heatmap(ng)

}


cat(
"\nDONE\n"
)
