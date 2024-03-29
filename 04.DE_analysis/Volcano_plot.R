library(argparser, quietly=TRUE)

p <- arg_parser("create volcano plot")

# Add command line arguments
p <- add_argument(p, "--de_result", help = "DE_result obtained form DEseq2 or edgeR", type = "character")
p <- add_argument(p, "--de_method", help = "DEG method, DEseq2 or edgeR", type = "character")
p <- add_argument(p, "--select_genes", help = "a file contain GeneID to label on the plot, one gene one row", type = "character")
p <- add_argument(p, "--padj_cutoff", help = "significance level", type = "numeric", default = 0.05)
p <- add_argument(p, "--log2FC_cutoff", help = "if absolute value of log2FC is larger than log2FC_cutoff, It's a DEG, default: 1", type = "numeric", default = 1)
p <- add_argument(p, "--height", help="volcano plot height, units: 'in'", type="numeric", default = 5)
p <- add_argument(p, "--width", help="volcano plot width, units: 'in'", type="numeric", default = 6)
p <- add_argument(p, "--heatmap", help = "whether to draw a heatmap", type = "logical", default = FALSE)
# Parse the command line arguments
argv <- parse_args(p)

filename <- argv$de_result
de_method <- argv$de_method
select_geneID <- argv$select_genes
padj_cutoff <- argv$padj_cutoff
log2FC_cutoff <- argv$log2FC_cutoff
width <- argv$width
height <- argv$height

test <- FALSE
if(test){
  filename <- "./WT_vs_MUT.DESeq2.DE_results.txt"
  de_method <- "DEseq2"
  #select_geneID <- "./select.txt"
  padj_cutoff <- 0.05
  log2FC_cutoff <- 1
  width <- 6
  height <- 5
}


library(tidyverse)
library(ggsci)
library(cowplot)

if (de_method == "DEseq2") {
  logFC <- "log2FoldChange"
  padj <- "padj"
  if (require(ggtext)) {
    ylab <- "-log10 adjusted *p*-value"
  } else {
    ylab <- "-log10 adjusted p-value"
  }
} else if (de_method == "edgeR") {
  logFC <- "logFC"
  padj <- "FDR"
  if (require(ggtext)) {
    ylab <- "*FDR*"
  } else {
    ylab <- "FDR"
  }
}
# import data
de_result <- read.table(file = filename) %>% rownames_to_column("GeneID") %>%
  select(GeneID, sampleA, sampleB, log2FoldChange = all_of(logFC), padj = all_of(padj)) %>%
  mutate(direction = if_else(padj > padj_cutoff | abs(log2FoldChange) < log2FC_cutoff, "NS",
                             if_else(log2FoldChange >= log2FC_cutoff, "UP", "DOWN")))
stat <- de_result %>% group_by(direction) %>% count()
up <- stat %>% filter(direction == "UP") %>% pull(n)
down <- stat %>% filter(direction == "DOWN") %>% pull(n)
# plot
Pvolcano <- ggplot(data = de_result, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 2, aes(color = direction), show.legend = F) +
  geom_hline(yintercept = -log10(padj_cutoff), 
             linetype = 'dotdash', color = 'grey30') +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), 
             linetype = 'dotdash', color = 'grey30') + 
  labs(y = ylab) +
  annotate("text", x = max(de_result$log2FoldChange), y = max(-log10(de_result$padj)), 
           label = paste("Up ", up, "\n", "Down ", down, sep = ""), 
           vjust=1, hjust=1, colour="black", size=4) +
  scale_color_manual(values = c('#1613BF', '#727272', '#D5161A'), breaks = c("DOWN", "NS", "UP")) +
  theme_half_open()
# add selected gene
if (!is.na(select_geneID)) {
  library(ggrepel)
  selected_genes <- read_tsv(file = select_geneID, col_names = c("GeneID", "GeneName"))
  res_selected <- selected_genes %>% left_join(de_result, by = "GeneID")
  Pvolcano <- Pvolcano + 
    geom_point(data = res_selected, size = 2, shape = 21, stroke = 1.5) +
    geom_text_repel(data = res_selected, aes(label = GeneName), fontface='italic', 
                    nudge_x = ifelse(res_selected$log2FoldChange > 0, 1.5, -1.5), 
                    nudge_y = ifelse(res_selected$log2FoldChange > 0, 10, 10), force = 10)
}

if(require(ggtext)){
  p <- Pvolcano + theme(axis.title.y = ggtext::element_markdown())
}else{
  p <- Pvolcano
}

ggsave(p, filename = str_replace(filename, "DE_results.txt", "volcano.pdf"), width = width, height = height)
ggsave(p, filename = str_replace(filename, "DE_results.txt", "volcano.png"), width = width, height = height, dpi = 500)

if (argv$heatmap) {
  expr <- read.table(file = "../03.Merge_result/genes.TMM.EXPR.matrix")
  sampleInfo <- read_tsv("../00.data/samples.txt", col_names = c("Group", "Sample", "R1Path", "R2Path"))
  
  de_GeneID <- de_result %>% filter(direction %in% c("UP", "DOWN")) %>% pull(GeneID)
  sample <- tibble(Group = c(unique(de_result$sampleA), unique(de_result$sampleB))) %>% 
    left_join(sampleInfo, by = "Group") %>% pull(Sample)
  p_HM <- pheatmap::pheatmap(log10(expr[de_GeneID,sample]+1),
                cluster_cols = F,
                cluster_rows = T,
                show_rownames = F,
                scale = "row",
                fontsize_row = 6,
                color = colorRampPalette(c("#758EB7", "#FFFFFF", "#FF7B89"))(100),
                border_color = "white",
                angle_col = 45,
  )
  ggsave(p_HM, filename = str_replace(filename, "DE_results.txt", "heatmap.pdf"), width = 3.5, height = 6)
  ggsave(p_HM, filename = str_replace(filename, "DE_results.txt", "heatmap.png"), width = 3.5, height = 6, dpi = 500)
}
