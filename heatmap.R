#!/usr/bin/env Rscript

# KEGG 7 大类与 MSS 的 UCell 趋势热图
#
# 输入：
#   1. KEGG59_UCell_29974_per_cell_cache.rds
#      已计算的 59 条通路逐细胞 UCell 分数。
#   2. Figure3_malignant_MSS_1009_final.csv
#      细胞名称及 MSS。
#   3. Figure3_KEGG_7groups_59pathways.csv
#      59 条通路与 7 大类的对应关系。
#
# 汇总：
#   患者内按 MSS 排序 -> 50 个近似等细胞量箱
#   -> 每条通路的箱内均值 -> 每个大类的通路等权均值
#   -> 患者等权均值 -> 行 Z-score 与 50 箱 Spearman 相关性。
#
# 注意：
#   不同患者相同箱号代表相对 MSS 排名，不代表相同 MSS 数值区间。
#   rho_bin 是 50 个汇总点的描述性相关系数。
#   本脚本读取已有 UCell 缓存，不重新计算逐细胞 UCell 分数。
#
# 启动 R 前必须设置临时目录。Bash 运行示例：
#   mkdir -p /mnt/data/linyanlin/Fig5/2/tmp
#   TMPDIR=/mnt/data/linyanlin/Fig5/2/tmp \
#   TMP=/mnt/data/linyanlin/Fig5/2/tmp \
#   TEMP=/mnt/data/linyanlin/Fig5/2/tmp \
#   SQLITE_TMPDIR=/mnt/data/linyanlin/Fig5/2/tmp \
#   Rscript KEGG7_MSS_patientEqual_UCell_50bins.R

# 1. 加载依赖 ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# 2. 设置路径 ---------------------------------------------------------------

root <- "/mnt/data1/linyanlin/scATAC/picture/pic1/heatmap/hall/MSS"
setwd(root)

cache <- file.path(
  root,
  "KEGG59_UCell_29974_per_cell_cache.rds"
)

meta_file <- paste0(
  "/mnt/data1/linyanlin/scATAC/picture/pic1/heatmap/hall/group/",
  "Figure3_malignant_MSS_1009_final.csv"
)

mapping_file <- paste0(
  "/mnt/data1/linyanlin/scATAC/picture/pic1/heatmap/hall/group/",
  "Figure3_KEGG_7groups_59pathways.csv"
)

stopifnot(file.exists(cache))

# 3. 读取 MSS，提取患者编号并排序 -------------------------------------------

meta <- fread(
  meta_file,
  select = c("cell", "MSS")
)

meta <- meta[is.finite(MSS)]

stopifnot(
  nrow(meta) == 29974L,
  !anyDuplicated(meta$cell)
)

meta[, patient := sub(
  "^.*_(HCC[0-9]+)_.*$",
  "\\1",
  cell
)]

stopifnot(
  uniqueN(meta$patient) == 10L,
  all(grepl("^HCC[0-9]+$", meta$patient))
)

setorder(meta, patient, MSS, cell)

# 4. 读取 UCell 缓存，按细胞名称匹配 ----------------------------------------

sc <- as.data.table(readRDS(cache))

stopifnot(
  "cell" %in% names(sc),
  !anyDuplicated(sc$cell),
  setequal(sc$cell, meta$cell)
)

paths <- setdiff(names(sc), "cell")
stopifnot(length(paths) == 59L)

setkey(sc, cell)

dat <- cbind(
  meta,
  sc[J(meta$cell), ..paths]
)

# 5. 每个患者内按 MSS 分为 50 个近似等细胞量箱 -----------------------------

dat[
  ,
  group_id := floor((seq_len(.N) - 1L) * 50L / .N) + 1L,
  by = patient
]

# 每个患者、每个箱内，分别计算 59 条通路的 UCell 均值。
pb <- dat[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = .(patient, group_id),
  .SDcols = paths
]

# 6. 读取 59 条通路与 7 大类的对应关系 -------------------------------------

mp <- unique(
  fread(
    mapping_file,
    select = c("group", "gmt_geneset_name")
  ),
  by = "gmt_geneset_name"
)

setnames(mp, "gmt_geneset_name", "pathway")

u <- melt(
  pb,
  id.vars = c("patient", "group_id"),
  measure.vars = paths,
  variable.name = "pathway",
  value.name = "score",
  variable.factor = FALSE
)

u <- merge(
  u,
  mp,
  by = "pathway",
  sort = FALSE
)

# 7. 先对大类内通路等权，再对患者等权 -------------------------------------

# 每个患者、每个箱、每个大类：对该类内通路的均值再取均值。
pc <- u[
  ,
  .(patient_class_score = mean(score, na.rm = TRUE)),
  by = .(patient, group_id, group)
]

# 每个大类、每个箱：对患者等权取均值。
ca <- pc[
  ,
  .(
    UCell_mean = mean(patient_class_score),
    n_patients = .N
  ),
  by = .(group, group_id)
]

# 横轴 MSS 也先计算患者内箱均值，再按患者等权汇总。
mm <- dat[
  ,
  .(MSS_mean = mean(MSS)),
  by = .(patient, group_id)
][
  ,
  .(MSS_mean = mean(MSS_mean)),
  by = group_id
][order(group_id)]

ca <- merge(
  ca,
  mm,
  by = "group_id",
  sort = FALSE
)

# 8. 计算行 Z-score 与 50 箱 Spearman 相关系数 -----------------------------

# 每个大类在 50 个箱之间进行标准化；标准差为 0 时填 0。
ca[
  ,
  UCell_z := {
    s <- sd(UCell_mean)

    if (is.finite(s) && s > 0) {
      (UCell_mean - mean(UCell_mean)) / s
    } else {
      rep(0, .N)
    }
  },
  by = group
]

ca[
  ,
  rho_bin := if (uniqueN(UCell_mean) > 1L) {
    suppressWarnings(
      cor(MSS_mean, UCell_mean, method = "spearman")
    )
  } else {
    NA_real_
  },
  by = group
]

rk <- unique(
  ca[, .(group, rho_bin)],
  by = "group"
)

setorder(rk, -rho_bin)

rk[
  ,
  rank_positive := frank(
    -rho_bin,
    ties.method = "min",
    na.last = "keep"
  )
]

# 9. 导出汇总表与相关趋势排序表 -------------------------------------------

fwrite(
  ca[order(group, group_id)],
  file.path(root, "KEGG7_MSS_patientEqual_UCell_50bins.csv")
)

fwrite(
  rk,
  file.path(root, "KEGG7_MSS_patientEqual_UCell_50bins_rank.csv")
)

# 10. 构建热图矩阵 ---------------------------------------------------------

wide <- dcast(
  ca,
  group ~ group_id,
  value.var = "UCell_z"
)

wide <- wide[match(rk$group, wide$group)]

heat <- as.matrix(
  wide[, as.character(1:50), with = FALSE]
)

rownames(heat) <- as.character(wide$group)

rv <- rk[match(rownames(heat), group), rho_bin]

rh <- matrix(
  rv,
  ncol = 1,
  dimnames = list(rownames(heat), "rho")
)

# 11. 构建底部 MSS 注释 ----------------------------------------------------

x <- mm$MSS_mean
mid <- median(x)

if (!(min(x) < mid && mid < max(x))) {
  mid <- (min(x) + max(x)) / 2
}

ha <- HeatmapAnnotation(
  MSS = x,
  col = list(
    MSS = colorRamp2(
      c(min(x), mid, max(x)),
      c("#2C7FB8", "#D9D9D9", "#F03B20")
    )
  ),
  annotation_height = unit(4, "mm"),
  annotation_name_side = "left"
)

# 12. 绘制 UCell 行 Z-score 热图 ------------------------------------------

ht <- Heatmap(
  heat,
  name = "UCell row Z-score",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#6C7FB8", "#F3EEE8", "#F04A3A")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 13),
  row_names_max_width = unit(10, "cm"),
  width = unit(8.2, "in"),
  height = unit(1.15, "in"),
  rect_gp = gpar(col = "white", lwd = 0.35),
  bottom_annotation = ha,
  column_title = paste0(
    "KEGG 7 classes | Patient-equal weighted UCell | 50 MSS bins"
  ),
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  )
)

# 右侧单列展示各大类与 MSS 的 Spearman 相关系数。
ht2 <- Heatmap(
  rh,
  name = "rho",
  col = colorRamp2(
    c(-1, 0, 1),
    c("#6C7FB8", "#F3EEE8", "#F04A3A")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  width = unit(10, "mm")
)

# 13. 导出 PDF -------------------------------------------------------------

out <- file.path(
  root,
  "KEGG7_MSS_patientEqual_UCell_50bins_squareish.pdf"
)

pdf(
  out,
  width = 13.5,
  height = 4.8,
  useDingbats = FALSE
)

tryCatch(
  draw(
    ht + ht2,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(6, 8, 6, 10), "mm")
  ),
  finally = dev.off()
)

cat("\n===== 7大类相关趋势：从正到负 =====\n")
print(rk, nrows = Inf)
cat("\nPDF:", out, "\n完成\n")

