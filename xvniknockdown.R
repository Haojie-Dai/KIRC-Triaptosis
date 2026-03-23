library(Seurat)
library(dplyr)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(remotes)
library(scTenifoldKnk)
library(svglite)
library(ggplot2)
# 假设 scobj 和 count 已经存在
scobj <- readRDS("scobj.rds")##scobj就是经理数据处理后的seurat文件
count <- GetAssayData(scobj, layer = "counts")  # 或 slot = "counts"

DefaultAssay(scobj) <- "RNA"

# 1) 找高变基因（你可以试 2000 / 3000 / 4000）
n_hvg <- 3000
scobj <- FindVariableFeatures(scobj, nfeatures = n_hvg)
hvgs  <- VariableFeatures(scobj)

# 2) 只保留这些高变基因
genes_use <- intersect(rownames(count), hvgs)
count_sub <- count[genes_use, , drop = FALSE]

# 3) 去掉几乎不表达的基因（至少在 min_cells 个细胞中表达）
min_cells <- 20
keep_genes <- rowSums(count_sub > 0) >= min_cells
count_sub  <- count_sub[keep_genes, , drop = FALSE]

cat("精简后基因数：", nrow(count_sub), "\n")  # 建议控制在 2000~4000 左右


set.seed(123)

max_cells <- 5000  # 你可以按自己机器内存试 4000 / 6000 / 8000

if (ncol(count_sub) > max_cells) {
  cells_use <- sample(colnames(count_sub), max_cells)
  count_sub <- count_sub[, cells_use, drop = FALSE]
}

cat("精简后细胞数：", ncol(count_sub), "\n")
# scobj 如果后面暂时不用，可以先删掉
rm(scobj, count)
gc()

library(scTenifoldKnk)

gKO <- "所敲除基因的基因symbol"

set.seed(133)

result <- scTenifoldKnk(
  countMatrix    = count_sub,
  gKO            = gKO,
  qc_mtThreshold = 0.1,
  qc_minLSize    = 1000,
  nc_nNet        = 5,    # 先从 5 个子网络开始，比 10 要省很多
  nc_nCells      = 500,  # 每个子网络采样 400 个细胞
  nc_nComp       = 2     # 主成分数量先设 2
)

saveRDS(result, "scTenifoldKnk_result_VCAN_shrink.rds")
result=readRDS("scTenifoldKnk_result_VCAN_shrink.rds")
# 提取差异表达基因
diff_reg <- result$diffRegulation

#按 FC 从大到小排序，取前 20 个基因
top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 27)



p <- top_genes %>%
  mutate(gene = fct_reorder(gene, FC)) %>%       # 按 FC 排序
  ggplot(aes(x = gene, y = FC, fill = FC)) +
  geom_col(width = 0.7, show.legend = FALSE) +   # 柱状图，柱宽0.7
  coord_flip() +
  geom_text(aes(label = round(FC, 1)),
            hjust = -0.15, size = 3.6, color = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
  labs(
    title = "Differentially Regulated Genes",
    subtitle = paste0("Ranked by fold change (n = ", nrow(top_genes), ")"),
    x = NULL, y = "Fold Change (FC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40")
  )
p
ggsave(
  filename = "top_genes_barplot.svg",
  plot = p,
  device = svglite,
  width = 7,
  height = 6,
  units = "in"
)


df <- result$diffRegulation
pretty_z_p_plot <- function(df,
                            z_col = "Z",
                            p_col = "p.adj",
                            gene_col = "gene",
                            z_cut = 2,
                            p_cut = 0.01,
                            p_line = 0.05,       # 参考线（可与 p_cut 不同）
                            max_label = 60) {
  
  # 1) 预处理：-log10(p) & 分组
  df <- df %>%
    mutate(
      log_pval = -log10(pmax(.data[[p_col]], 1e-300)), # 避免 Inf
      group = case_when(
        .data[[p_col]] < p_cut & abs(.data[[z_col]]) >= z_cut ~ "Hit (|Z|≥2 & q<0.01)",
        abs(.data[[z_col]]) >= z_cut ~ "|Z|≥2 only",
        .data[[p_col]] < p_cut ~ "q<0.01 only",
        TRUE ~ "Not significant"
      )
    )
  
  # 2) 需要标注的基因
  label_df <- df %>%
    filter(.data[[p_col]] < p_cut, abs(.data[[z_col]]) >= z_cut)
  
  # 3) 色板
  pal <- c("Not significant" = "#B8B8B8",
           "|Z|≥2 only"     = "#7AA6C2",
           "q<0.01 only"    = "#C49A6C",
           "Hit (|Z|≥2 & q<0.01)" = "#E15759")
  
  # 4) 阈值矩形范围（左右两侧 |Z|≥z_cut 且 q<p_cut 的象限）
  y_thr  <- -log10(p_cut)
  x_thr  <- z_cut
  
  p <- ggplot(df, aes(x = .data[[z_col]], y = log_pval, color = group)) +
    
    # 高亮显著区域（右侧）
    annotate("rect",
             xmin =  x_thr, xmax =  Inf,
             ymin =  y_thr, ymax =  Inf,
             alpha = 0.06, fill = "#E15759") +
    # 左侧
    annotate("rect",
             xmin = -Inf, xmax = -x_thr,
             ymin =  y_thr, ymax =  Inf,
             alpha = 0.06, fill = "#E15759") +
    
    # 散点
    geom_point(size = 1.9, alpha = 0.75, stroke = 0) +
    
    # 参考线：p 和 Z
    geom_hline(yintercept = -log10(p_line), linetype = "dashed", color = "#E15759", linewidth = 0.6) +
    geom_hline(yintercept =  y_thr,        linetype = "dotted", color = "#E15759", linewidth = 0.6) +
    geom_vline(xintercept =  c(-x_thr, x_thr), linetype = "dashed", color = "#2B6CB0", linewidth = 0.6) +
    
    # 标注
    ggrepel::geom_text_repel(
      data = head(label_df %>% arrange(desc(log_pval)), max_label),
      aes(label = .data[[gene_col]]),
      size = 3.3, seed = 123,
      box.padding = 0.35, point.padding = 0.2,
      max.overlaps = Inf, min.segment.length = 0,
      segment.color = alpha("grey20", 0.6), segment.size = 0.3
    ) +
    
    # 标题与轴
    labs(
      title = "Z vs -log10(q-value)",
      subtitle = glue::glue("Thresholds: |Z| ≥ {z_cut}, q < {p_cut}  ·  Hits = {nrow(label_df)}"),
      x = "Z-score", y = expression(-log[10](q))
    ) +
    
    # 颜色
    scale_color_manual(values = pal, name = NULL) +
    
    # 坐标边距更紧凑
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    
    # 主题
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(face = "bold")
    )
  
  return(p)
}

p1 <- pretty_z_p_plot(df,
                      z_col = "Z",
                      p_col = "p.adj",
                      gene_col = "gene",
                      z_cut = 2,
                      p_cut = 0.01,
                      p_line = 0.05,
                      max_label = 50)
p1
ggsave(
  filename = "top_genes_vocno.svg",
  plot = p1,
  device = svglite,
  width = 6,
  height = 6,
  units = "in"
)


# install.packages(c("patchwork", "svglite"))  # 若未安装
library(ggplot2)
library(patchwork)
library(svglite)

# 上下拼接：p 在上，p1 在下
#combined <- p / p1

# 如果想左右拼接：p | p1
combined <- p | p1

ggsave(
  filename = "combined.svg",
  plot = combined,
  device = svglite::svglite,
  width = 16, height = 8, units = "in"
)
write.csv(top_genes, "top_genes.csv", row.names = FALSE, fileEncoding = "UTF-8")

