#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("do GO/KEGG enrichment, need prepare kegg_info.RData(containing ko2gene pathway2gene and pathway2name) and org.*.eg.db_1.0.tar.gz in $work_dir/db")

# Add command line arguments
p <- add_argument(p, "--de_result", help="input de_result file, from run_DE_analysis.pl", type="character")
p <- add_argument(p, "--de_log2FoldChange", help="log2FoldChange cutoff", type="numeric", default = 1)
p <- add_argument(p, "--de_padj", help="adjust pvalue cutoff", type="numeric", default = 0.05)
p <- add_argument(p, "--enrich_pvalue", help="pvalue cutoff for enrichment", type="numeric", default = 0.05)
p <- add_argument(p, "--enrich_qvalue", help="qvalue cutoff for enrichment", type="numeric", default = 0.05)
p <- add_argument(p, "--orgdb", help="orgdb package name, must have been installed in ./R_library", type="character")
p <- add_argument(p, "--species", help="character, either the kegg code, scientific name or the common name of the target species.", type="character", default = "ko")


# Parse the command line arguments
argv <- parse_args(p)

########################################################################
## test
#argv <- list()
#argv$de_result <- ""
#argv$de_log2FoldChange <- 1
#argv$de_padj <- 0.05
#argv$enrich_pvalue <- 0.05
#argv$enrich_qvalue <- 0.05
#argv$orgdb <- ""
#argv$species <- "ko"
test <- FALSE
if (test) {
  argv <- list()
  argv$de_result <- "../04.DE_analysis/genes.counts.matrix.WT_vs_WTE.DESeq2.DE_results"
  argv$de_log2FoldChange <- 1
  argv$de_padj <- 0.05
  argv$enrich_pvalue <- 0.05
  argv$enrich_qvalue <- 0.05
  argv$orgdb <- "org.Bnapus.eg.db"
  argv$species <- "bna"
}
########################################################################

# load library ------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(enrichplot)

out_prefix <- str_replace(string = basename(argv$de_result), pattern = "\\.csv$|\\.txt$|\\.tsv$", replacement = "")
orgdb <- argv$orgdb

if (!dir.exists("../db/R_Library/")) {
  dir.create('../db/R_Library', recursive = T)
}
if (!requireNamespace(orgdb, lib.loc = "../db/R_Library/", quietly = TRUE)) {
  orgpkg <- list.files(path = "../db/", pattern = "org.*eg.db_.*tar\\.gz")[[1]]
  install.packages(paste("../db/", orgpkg, sep = ""), 
                   repos = NULL, #从本地安装
                   lib = '../db/R_Library') # 安装文件夹
}
library(orgdb, lib.loc = "../db/R_Library", character.only = TRUE)

#
source("./enrich_plot.R")

de_result <- read.table(file = argv$de_result, header = TRUE, row.names = 1)
gene <- filter(de_result,
               abs(log2FoldChange) > argv$de_log2FoldChange & padj < argv$de_padj) %>%
  rownames(id)

geneList <- de_result$log2FoldChange
names(geneList) <- rownames(de_result)
geneList <- sort(geneList, decreasing = TRUE)

de_ego <- enrichGO(gene = gene,
                   OrgDb = orgdb,
                   keyType = 'GID',
                   ont = 'ALL',
                   qvalueCutoff = argv$enrich_qvalue,
                   pvalueCutoff = argv$enrich_pvalue)
try(expr = {de_ego <- clusterProfiler::simplify(de_ego, cutoff=0.7, by="p.adjust", select_fun=min)})
de_ego_df <- as.data.frame(de_ego)
head(de_ego_df)
write_csv(de_ego_df, paste(out_prefix, "GO_result", "csv", sep = "."))
write_tsv(de_ego_df, paste(out_prefix, "GO_result", "txt", sep = "."))

## GO绘图
# 
if (nrow(de_ego_df)>0) {
  #GO_emap <- emapplot(pairwise_termsim(de_ego), cex_label_category=.8, cex_line=.5) + 
  #  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
  #                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  #GO_emap <- emapplot(pairwise_termsim(de_ego), cex.params = list(category_label = .8, line = .5)) + 
  #  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
  #                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  # Bioconductor 升级到 version 3.20 以后，参数名称修改，
  # 为啥要修改参数名啊，功能难道有不同吗？用了好多年老改参数很不方便啊
  GO_emap <- emapplot(pairwise_termsim(de_ego), size_category=.8, size_edge=.5) + 
    scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                          guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  ggsave(GO_emap, filename = paste(out_prefix, "GO_emapplot.pdf", sep = "."), width = 10, height = 9)
  ggsave(GO_emap, filename = paste(out_prefix, "GO_emapplot.png", sep = "."), width = 10, height = 9, dpi = 500)
  
  plotGO(de_ego_df, out_prefix, width = 7.5)
}

# GO DAG plot
for (ont in c("BP", "MF", "CC")) {
  tryCatch(
    {
      #ego <- enrichGO(gene = gene,
      #                OrgDb = orgdb,
      #                keyType = 'GID',
      #                ont = ont,
      #                qvalueCutoff = argv$enrich_qvalue,
      #                pvalueCutoff = argv$enrich_pvalue)
      #ego <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      ego <- de_ego %>% clusterProfiler::filter(ONTOLOGY == ont)
      ego@ontology <- ont
      pdf(file = paste(out_prefix, paste("GO", ont, "DAG.pdf", sep = "_"), sep = "."), width = 10, height = 10)
      plotGOgraph(ego)
      dev.off()
      png(filename = paste(out_prefix, paste("GO", ont, "DAG.png", sep = "_"), sep = "."), width = 10, height = 10, units = "in", res = 500)
      plotGOgraph(ego)
      dev.off()
    },
    error = function(e){print("may no term\n")}
  )
}

## KEGG
load(paste("../db/", argv$species, ".kegg_info.RData", sep = ""))
de_ekp <- enricher(gene,
                   TERM2GENE = pathway2gene,
                   TERM2NAME = pathway2name,
                   pvalueCutoff = argv$enrich_pvalue,
                   qvalueCutoff = argv$enrich_qvalue)
de_ekp_df <- as.data.frame(de_ekp)
head(de_ekp_df)
write_csv(de_ekp_df, paste(out_prefix, "KEGG_result", "csv", sep = "."))
write_tsv(de_ekp_df, paste(out_prefix, "KEGG_result", "txt", sep = "."))

# KEGG绘图
if (nrow(de_ekp_df)>0) {
  KEGG_emap <- emapplot(pairwise_termsim(de_ekp), size_category=.8, size_edge=.5) + 
    scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                          guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  ggsave(KEGG_emap, filename = paste(out_prefix, "KEGG_emapplot.pdf", sep = "."), width = 10, height = 9)
  ggsave(KEGG_emap, filename = paste(out_prefix, "KEGG_emapplot.png", sep = "."), width = 10, height = 9, dpi = 500)
  
  plotKEGG(de_ekp_df, out_prefix)
}

if(FALSE){
  # pathway view ------------------------------------------------------------
  # 此处并不完美，并不是对应物种专用的pathway，有改善空间
  #id.map <- select(org.My.eg.db, keys = names(geneList), columns = "Ko")
  #gene.ko <- mol.sum(mol.data = geneList, id.map = ko2gene)
  gene.ko <- geneList
  id <- tibble(GID = names(gene.ko)) %>% left_join(ko2gene, by = "GID")
  names(gene.ko) <- id$KO
  gene.ko <- as.matrix(gene.ko[!is.na(names(gene.ko))])
  
  sig.pathway <- as.character(filter(de_ekp_df, p.adjust < argv$enrich_qvalue)$ID)
  
  work_dir <- getwd()
  pathview_dir <- paste(out_prefix, 'pathwiew', sep = "_")
  try(dir.create(pathview_dir, recursive=T))
  setwd(pathview_dir)
  
  result <- tryCatch(pathview(gene.data  = gene.ko,
           pathway.id = str_replace(sig.pathway, "bna", "ko"),
           species    = "ko"), error = function(e){"Can't download pathway map, maybe network is not available!\n"})
  
  setwd(work_dir)
}
