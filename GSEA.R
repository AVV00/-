```{r}
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(patchwork)

# ==== GSEA：使用全基因列表（ProxJej Sham vs SHAMWM）====
all_file <- "proxjej_SHAMWM_DE_results_all_annotated.csv"
res_all <- read.csv(all_file)

# ==== 构建 geneList（log2FC 排序 + jitter + 降序）====
gene_list <- res_all %>%
  filter(!is.na(log2FoldChange) & !is.na(symbol)) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  mutate(log2FoldChange = jitter(log2FoldChange, amount = 1e-6)) %>%
  select(symbol, log2FoldChange)

gene_list <- setNames(gene_list$log2FoldChange, gene_list$symbol)
gene_list <- sort(gene_list, decreasing = TRUE)

# ==== gseGO 分析 ====
gsea_res <- gseGO(geneList = gene_list,
                  OrgDb = org.Rn.eg.db,
                  ont = "BP",
                  keyType = "SYMBOL",
                  pvalueCutoff = 0.1,
                  verbose = FALSE)

# ==== 提取前 20 项（按 abs(NES) 降序），处理 GO 名 ====
gsea_df <- as.data.frame(gsea_res) %>%
  arrange(desc(abs(NES))) %>%
  slice(1:20)
gsea_df$Description <- gsub("GO_", "", gsea_df$Description)
gsea_df$Description <- gsub("_", " ", gsea_df$Description)

# ==== 绘图：条形图 + 左对齐标题和 y 轴文字 ====
p2_core <- ggplot(gsea_df, aes(x = reorder(Description, NES), y = NES)) +
  geom_col(fill = "dodgerblue", width = 0.75) +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    x = "GO Term",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    axis.title.x = element_text(face = "bold", hjust = 0),
    axis.text.y = element_text(size = 11, face = "bold", hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ==== 设置标题 + 保存 ====
p2_final <- p2_core +
  plot_annotation(
    title = "GSEA GO Pathways (ProxJej Sham vs SHAMWM)"
  )

# 保存图像
ggsave("Barplot_GSEA_NES_ProxJej_Sham_vs_SHAMWM.png",
       plot = p2_final, width = 10, height = 6, dpi = 300)

# 打印图像
print(p2_final)

```
