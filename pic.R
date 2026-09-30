##############################################################火山图格式################################################################################################################
```{r}
library(ggplot2)
library(ggrepel)
library(dplyr)

# 要处理的组织名称
tissues <- c("Colon", "duodenum", "distjej", "ileum", "proxjej")

for (tissue in tissues) {

  input_file <- paste0(tissue, "_DE_results_all_annotated.csv")
  output_file <- paste0(tissue, "_volcano_pvalue.png")
  plot_title <- paste0(tissue, ": SHAMChow vs Sham ")

  res <- read.csv(input_file, row.names = 1)
  
  # 去除没有 symbol 的行
  res <- res %>% filter(!is.na(symbol))
  
  # 检查关键列
  if (!"symbol" %in% colnames(res) || all(is.na(res$symbol))) {
    warning(paste("Skipping", tissue, "- no valid gene symbols"))
    next
  }
  if (!"pvalue" %in% colnames(res)) {
    warning(paste("Skipping", tissue, "- no pvalue column found"))
    next
  }

  # 设置分组
  res$group <- "NotSig"
  res$group[res$log2FoldChange > 1 & res$pvalue < 0.05] <- "Up"
  res$group[res$log2FoldChange < -1 & res$pvalue < 0.05] <- "Down"
  res$group <- factor(res$group, levels = c("Down", "NotSig", "Up"))

  # 只取在 -5 到 5 范围内的 Up/Down 基因作为候选标注
  res_in_range <- res %>%
    filter(log2FoldChange > -5 & log2FoldChange < 5)

  top_up <- res_in_range %>%
    filter(group == "Up") %>%
    arrange(pvalue) %>%
    head(10)

  top_down <- res_in_range %>%
    filter(group == "Down") %>%
    arrange(pvalue) %>%
    head(10)

  label_data <- bind_rows(top_up, top_down) %>% filter(!is.na(symbol))

  # 自定义颜色
  color_map <- c("Up" = "firebrick", "Down" = "dodgerblue", "NotSig" = "gray70")

  # 绘图
  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = group)) +
    geom_point(alpha = 0.8, size = 3.5) +
    scale_color_manual(values = color_map) +
    geom_vline(xintercept = c(-1, 1), linetype = "twodash", color = "gray60", linewidth = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "twodash", color = "gray60", linewidth = 0.8) +
    geom_text_repel(data = label_data,
                    aes(label = symbol),
                    size = 5.5, box.padding = 0.5, max.overlaps = 20,
                    color = "black") +  # 黑色字体
    labs(title = plot_title,
         x = expression(Log[2]~"Fold Change"),
         y = expression(-Log[10]~"P-value"),
         color = "Group") +
    xlim(-5, 5) +
    theme_minimal() +
    theme(
      text = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.6)  # ✅ 添加的部分
    )

  # 保存图像
  ggsave(output_file, plot = p, width = 6.5, height = 9, dpi = 300)
  message("✅ Volcano plot saved for: ", tissue, 
          " | labeled Up: ", nrow(top_up), 
          " | Down: ", nrow(top_down))
  print(p)
}

###########################################################################################################################################################################################
