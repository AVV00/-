#!/usr/bin/env Rscript
# 1009 MSS matched malignant cells: gene.activity data -> PCA -> Harmony -> UMAP
# output: PDF only (new filename, no overwrite of old file)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(ggplot2)
  library(data.table)
  library(grid)
})

input_rds <- "/mnt/data1/linyanlin/scATAC/all_data/new_celltype-261009/all_cells_gene_activity_new_celltype_annotation.rds"
mss_csv   <- "/mnt/data1/linyanlin/scATAC/picture/pic1/heatmap/hall/1019_new/final2.csv"
out_dir   <- "/mnt/data1/linyanlin/scATAC/picture/pic3/umap"
tmp_dir   <- "/mnt/data/linyanlin/Fig5/2/tmp"

# 鏂板悕瀛楋細涓嶈鐩栨棫鏂囦欢
out_pdf <- file.path(out_dir, "1009_MSS_matched_UMAP_blue_red_largerLegend.pdf")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir, SQLITE_TMPDIR = tmp_dir)

stopifnot(file.exists(input_rds), file.exists(mss_csv))
set.seed(1234)

cat("[1/6] Read MSS final.csv\n")
meta <- fread(mss_csv, select = c("cell", "patient", "MSS"))
meta <- meta[!is.na(cell) & nzchar(cell) & is.finite(MSS)]
stopifnot(!anyDuplicated(meta$cell))

cat("[2/6] Read Seurat object\n")
source_obj <- readRDS(input_rds)
stopifnot(inherits(source_obj, "Seurat"), "gene.activity" %in% Assays(source_obj))

matched <- meta$cell[meta$cell %in% colnames(source_obj)]
cat("Matched cells:", length(matched), "\n")
stopifnot(length(matched) >= 100L)

cat("[3/6] Extract gene.activity data layer only\n")
stopifnot("data" %in% SeuratObject::Layers(source_obj[["gene.activity"]]))
activity_all <- SeuratObject::LayerData(source_obj, assay = "gene.activity", layer = "data")
stopifnot(nrow(activity_all) > 100L, ncol(activity_all) == ncol(source_obj))

if (is.null(colnames(activity_all))) colnames(activity_all) <- colnames(source_obj)
if (is.null(rownames(activity_all))) rownames(activity_all) <- rownames(source_obj[["gene.activity"]])

stopifnot(
  identical(colnames(activity_all), colnames(source_obj)),
  identical(rownames(activity_all), rownames(source_obj[["gene.activity"]]))
)

activity <- activity_all[, matched, drop = FALSE]
stopifnot(ncol(activity) == length(matched))
if (!inherits(activity, "dgCMatrix")) activity <- as(activity, "dgCMatrix")
if (any(!is.finite(activity@x))) stop("Non-finite values found in gene.activity data layer")

rm(activity_all, source_obj); invisible(gc())
cat("Matrix:", nrow(activity), "genes x", ncol(activity), "cells\n")

cat("[4/6] Build temporary Seurat object\n")
obj <- CreateSeuratObject(counts = activity, assay = "gene.activity", min.cells = 0, min.features = 0)
obj <- SetAssayData(obj, assay = "gene.activity", slot = "data", new.data = activity)

idx <- match(colnames(obj), meta$cell)
stopifnot(!anyNA(idx))
obj$MSS <- meta$MSS[idx]
obj$patient <- as.character(meta$patient[idx])

DefaultAssay(obj) <- "gene.activity"

means <- Matrix::rowMeans(activity)
vars  <- Matrix::rowMeans(activity^2) - means^2
vars[!is.finite(vars)] <- 0
names(vars) <- rownames(activity)
features <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000L, sum(vars > 0)))]
stopifnot(length(features) >= 50L)
VariableFeatures(obj) <- features

rm(activity, means, vars); invisible(gc())

cat("[5/6] ScaleData / PCA / Harmony / UMAP\n")
obj <- ScaleData(obj, features = features, verbose = FALSE)
obj <- RunPCA(obj, features = features, npcs = 30, verbose = FALSE, seed.use = 1234)
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = "patient",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  project.dim = FALSE,
  verbose = TRUE
)

stopifnot("harmony" %in% Reductions(obj))

obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 1234,
  umap.method = "uwot",
  metric = "cosine",
  verbose = TRUE
)

cat("[6/6] Prepare plot and save PDF only\n")
coords <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "cell")
setnames(coords, names(coords)[2:3], c("UMAP_1", "UMAP_2"))
coords[, MSS := obj$MSS[match(cell, colnames(obj))]]
stopifnot(nrow(coords) == length(matched), all(is.finite(coords$MSS)))

# 浣?MSS 鍏堢敾锛岄珮 MSS 鍚庣敾锛岄伩鍏嶉珮鍊肩偣琚洊浣?setorder(coords, MSS)

mss_min <- min(coords$MSS, na.rm = TRUE)
mss_max <- max(coords$MSS, na.rm = TRUE)
hi_cut <- as.numeric(stats::quantile(coords$MSS, 0.90, na.rm = TRUE))
coords_hi <- coords[MSS >= hi_cut]

p <- ggplot(coords, aes(UMAP_1, UMAP_2, colour = MSS)) +
  geom_point(size = 0.40, alpha = 0.70, stroke = 0) +
  geom_point(data = coords_hi, aes(UMAP_1, UMAP_2, colour = MSS), size = 0.75, alpha = 1, stroke = 0) +
  scale_colour_gradient2(
    low = "#3B82FF",
    mid = "#FFFFFF",
    high = "#D62728",
    midpoint = (mss_min + mss_max) / 2,
    limits = c(mss_min, mss_max),
    breaks = c(mss_min, mss_max),
    labels = c("Low", "High"),
    name = "MSS"
  ) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 13) +
  theme(
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.key.height = unit(40, "mm"),
    legend.key.width  = unit(7, "mm"),
    legend.spacing.y  = unit(4, "mm"),
    plot.margin = margin(18, 28, 18, 22)
  ) +
  guides(
    colour = guide_colourbar(
      barwidth  = unit(7, "mm"),
      barheight = unit(42, "mm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  )

grDevices::pdf(out_pdf, width = 8.5, height = 7, useDingbats = FALSE, onefile = FALSE)
tryCatch(print(p), finally = grDevices::dev.off())

cat("Done.\n")
cat("Output PDF: ", out_pdf, "\n", sep = "")
cat("Cells plotted: ", nrow(coords), "\n", sep = "")
