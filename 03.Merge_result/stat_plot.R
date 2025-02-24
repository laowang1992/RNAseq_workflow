# 加载扩展???
library(tidyverse)
library(ggplot2)
library(ggsci)
library(cowplot)
# 读取数据并转换格???
gene_exp <- read.table(file = "genes.TMM.EXPR.matrix", 
                       header = T, row.names=1)
head(gene_exp)
d <- gene_exp %>% rownames_to_column("GeneID") %>% as_tibble() %>% gather(key = sample, value = tmm, -GeneID)

colourCount = ncol(gene_exp)
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))

p1 <- ggplot(data=d, aes(x=sample, y=log10(tmm+1), fill=sample)) +
  #geom_violin(size=.5) +
  geom_boxplot(width=.5, size=.5) +
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(x = NULL, y = "log10(TMM+1)", fill = NULL) +
  #scale_fill_aaas() +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p1, filename="genes.TMM.EXPR.boxplot.pdf", height = 5, width = ncol(gene_exp)*0.4+2)
ggsave(p1, filename="genes.TMM.EXPR.boxplot.png", height = 5, width = ncol(gene_exp)*0.4+2, dpi = 500)

p2 <-ggplot(data=d, aes(x=log10(tmm+1), color=sample)) +
  geom_density(linewidth=1) +
  scale_color_manual(values = getPalette(colourCount)) +
  labs(x = "log10(TMM+1)", y = "Density", color = NULL) +
  theme_half_open()
ggsave(p2, filename="genes.TMM.EXPR.density.pdf", height = 4, width = 6)
ggsave(p2, filename="genes.TMM.EXPR.density.png", height = 4, width = 6, dpi = 500)


metadata <- read.table(file = "../00.data/samples.txt")
rownames(metadata) <- metadata$V2
metadata <- metadata %>% select(Group = V1)
metadata$Group <- factor(metadata$Group)

# PCA
library(PCAtools)
p <- pca(log10(gene_exp+1), metadata = metadata)
pdf(file = "screeplot.pdf", width = 6, height = 5)
screeplot(p)
dev.off()
png(file = "screeplot.png", width = 6, height = 5, units = "in", res = 500)
screeplot(p)
dev.off()
#biplot(p, 
#       x = 'PC1',                 # x ???
#       y = 'PC2',                 # y ???
#       colby = 'Strain',          # 颜色映射
#       shape = 'Stage',
#       legendPosition = 'right',  # 图例位置
#       lab = rownames(metadata)                    # 样本名称显示
#)

pca_rotated_plus <- rownames_to_column(p$rotated, 
                                       var = 'sample_name') %>%
  left_join(rownames_to_column(metadata, var = 'sample_name'), 
            by = 'sample_name')

p_pca <- ggplot(data = pca_rotated_plus, aes(x = PC1, y = PC2)) +
  geom_point(size = 6, aes(fill = Group), shape = 21) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(7, "Set2"))(length(unique(metadata$Group)))) +
  scale_color_npg() +
  labs(x = paste("PC1 (", round(p[["variance"]][["PC1"]], 2), "% variance explained)", sep = ""),
       y = paste("PC2 (", round(p[["variance"]][["PC2"]], 2), "% variance explained)", sep = "")) +
  #guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme_half_open() + 
  theme(#legend.position = c(0.2, 0.9), 
    legend.background = element_rect(fill = NA)
  )
ggsave(p_pca, filename = "PCAplot_PC1&PC2.pdf", height = 4.5, width = 5.6)
ggsave(p_pca, filename = "PCAplot_PC1&PC2.png", height = 4.5, width = 5.6, dpi = 500)

# 相关???
sample_cor <- cor(log10(gene_exp+1), 
                  method = "pearson") # 算法??? pearson | kendall | spearman
write.table(x = as.table(sample_cor), file = "./correlation.txt")
write.csv(x = as.table(sample_cor), file = "./correlation.csv")

color <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(unique(metadata$Group)))
names(color) <- unique(metadata$Group)
ann_colors = list(color)
names(ann_colors) <- "Group"

pdf("correlation_with_Cluster.pdf", width = nrow(sample_cor)*0.24+3, height = nrow(sample_cor)*0.24+1.5)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = T,          # 是否对行（基因）聚类
                   cluster_cols = T,          # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = metadata, # 添加列注释
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   annotation_names_col = F,
                   cellwidth = 12,
                   cellheight = 12,
                   fontsize = 6,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
)
dev.off()
png("correlation_with_Cluster.png", width = nrow(sample_cor)*0.24+3, height = nrow(sample_cor)*0.24+1.5, units = "in", res = 500)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = T,       # 是否对行（基因）聚类
                   cluster_cols = T,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = metadata, # 添加列注释
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   annotation_names_col = F,
                   cellwidth = 12,
                   cellheight = 12,
                   fontsize = 6,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
)
dev.off()


pdf("correlation_without_Cluster.pdf", width = nrow(sample_cor)*0.24+2, height = nrow(sample_cor)*0.24+0.5)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = F,       # 是否对行（基因）聚类
                   cluster_cols = F,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = metadata, # 添加列注释
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   annotation_names_col = F,
                   cellwidth = 12,
                   cellheight = 12,
                   fontsize = 6,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
)
dev.off()
png("correlation_without_Cluster.png", width = nrow(sample_cor)*0.24+2, height = nrow(sample_cor)*0.24+0.5, units = "in", res = 500)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = F,       # 是否对行（基因）聚类
                   cluster_cols = F,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = metadata, # 添加列注释
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   annotation_names_col = F,
                   cellwidth = 12,
                   cellheight = 12,
                   fontsize = 6,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
)
dev.off()

