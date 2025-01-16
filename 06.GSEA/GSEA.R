#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("GSEA, need prepare kegg_info.RData(containing ko2gene pathway2gene and pathway2name) and org.*.eg.db_1.0.tar.gz in $work_dir/db,")

# Add command line arguments
p <- add_argument(p, "--de_result", help = "input de_result file, from run_DE_analysis.pl", type="character")
p <- add_argument(p, "--gsePadj", help = "pvalue cutoff for enrichment", type="numeric", default = 0.01)
p <- add_argument(p, "--orgdb", help = "orgdb package name, must have been installed in ./R_library", type="character")
p <- add_argument(p, "--species", help = "character, either the kegg code, scientific name or the common name of the target species.", type="character", default = "ko")
p <- add_argument(p, "--drawPdf", help = "draw pdf image or not", type = "logical", default = "FALSE")

# Parse the command line arguments
argv <- parse_args(p)

########################################################################
## test
test <- F
if (test) {
  argv <- list()
  argv$de_result <- "../04.DE_analysis/genes.counts.matrix.WT_vs_WTE.DESeq2.DE_results" 
  argv$gsePadj <- 0.01
  argv$orgdb <- "org.Bnapus.eg.db"
  argv$species <- "bna"
  argv$drawPdf <- FALSE
  
  filename <- basename(argv$de_result)
  orgdb <- argv$orgdb
  gsePadj <- argv$gsePadj
  drawPdf <- argv$drawPdf
  species <- argv$species
}

########################################################################
## load R packages
library(clusterProfiler)
library(enrichplot)
library(tidyverse)

filename <- str_replace(string = basename(argv$de_result), pattern = "\\.csv$|\\.txt$|\\.tsv$", replacement = "")
orgdb <- argv$orgdb
gsePadj <- argv$gsePadj
drawPdf <- argv$drawPdf
species <- argv$species

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
source("./enrich_plot.R")

# read DE result
de_result <- read.table(file = argv$de_result, header = TRUE, row.names = 1)
# 
genes <- de_result$log2FoldChange
names(genes) <- rownames(de_result)
genes <- sort(genes, decreasing = T)

## GO
# GSEA
gseGO_res <- gseGO(geneList = genes,
                   OrgDb = get(orgdb),
                   ont = "ALL",
                   keyType = "GID",
                   #minGSSize = 5,
                   #maxGSSize = 1000,
                   pvalueCutoff = gsePadj,
                   pAdjustMethod = "BH")
try(expr = {gseGO_res <- clusterProfiler::simplify(gseGO_res, cutoff=0.7, by="p.adjust", select_fun=min)})
# write GO result
write_tsv(gseGO_res@result, paste(filename, "gseGO.txt", sep = "."))
write_csv(gseGO_res@result, paste(filename, "gseGO.csv", sep = "."))

if (nrow(gseGO_res@result)>0) {
  #GO_emap <- emapplot(pairwise_termsim(gseGO_res), cex_label_category=.8, cex_line=.5) + 
  #  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
  #                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  #GO_emap <- emapplot(pairwise_termsim(gseGO_res), cex.params = list(category_label = .8, line = .5)) + 
  #  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
  #                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  # Bioconductor 升级到 version 3.20 以后，参数名称修改，
  # 为啥要修改参数名啊，功能难道有不同吗？用了好多年老改参数很不方便啊
  GO_emap <- emapplot(pairwise_termsim(gseGO_res), size_category=.8, size_edge=.5) + 
    scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                          guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  ggsave(GO_emap, filename = paste(filename, "gseGO_emapplot.pdf", sep = "."), width = 10, height = 9)
  ggsave(GO_emap, filename = paste(filename, "gseGO_emapplot.png", sep = "."), width = 10, height = 9, dpi = 500)
  
  plotGO(gseGO_res@result, filename, width = 8)
}

# create GSEA GO dictionary and plot
if (!dir.exists(paste(filename, "GO", sep = "/"))) {
  dir.create(paste(filename, "GO", sep = "/"), recursive = T)
} else {
  # 指定文件夹路径  
  folder_path <- paste(filename, "GO", sep = "/") 
  
  # 获取文件夹中的文件列表  
  file_list <- list.files(folder_path)  
  
  # 遍历文件列表，删除每个文件  
  for (file in file_list) {  
    file_path <- file.path(folder_path, file)  
    file.remove(file_path)  
  }
}
for (i in seq_along(gseGO_res@result$ID)) {
  cat("Plot ", gseGO_res@result$ID[i], " ", gseGO_res@result$Description[i], "\n", sep = "")
    prefix <- str_replace(gseGO_res@result$ID[i], ":", "")
  p <- gseaplot2(gseGO_res, geneSetID = gseGO_res@result$ID[i], title = gseGO_res@result$Description[i])
  if(drawPdf){
    pdf(file = paste(filename, "/GO/", prefix, ".pdf", sep = ""), width = 6, height = 4.5)
    print(p)
    dev.off()
  }
  png(file = paste(filename, "/GO/", prefix, ".png", sep = ""), width = 6, height = 4.5, units = "in", res = 500)
  print(p)
  dev.off()
}

#gseaplot2(gseGO_res, geneSetID = 4, title = gseGO_res@result$Description[4])

## KEGG
load(paste("../db/", argv$species, ".kegg_info.RData", sep = ""))
# GSEA
gseKEGG_res <- GSEA(geneList = genes,
                    TERM2GENE = pathway2gene,
                    TERM2NAME = pathway2name,
                    #minGSSize = 5,
                    #maxGSSize = 1000,
                    pvalueCutoff = gsePadj,
                    pAdjustMethod = "BH")
# write KEGG result
write_tsv(gseKEGG_res@result, paste(filename, "gseKEGG.txt", sep = "."))
write_csv(gseKEGG_res@result, paste(filename, "gseKEGG.csv", sep = "."))

if (nrow(gseKEGG_res@result)>0) {
  KEGG_emap <- emapplot(pairwise_termsim(gseKEGG_res), size_category=.8, size_edge=.5) + 
    scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                          guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  ggsave(KEGG_emap, filename = paste(filename, "gseKEGG_emapplot.pdf", sep = "."), width = 10, height = 9)
  ggsave(KEGG_emap, filename = paste(filename, "gseKEGG_emapplot.png", sep = "."), width = 10, height = 9, dpi = 500)
  
  plotKEGG(gseKEGG_res@result, filename, width = 8)
}

# create GSEA KEGG dictionary and plot
if (!dir.exists(paste(filename, "KEGG", sep = "/"))) {
  dir.create(paste(filename, "KEGG", sep = "/"), recursive = T)
} else {
  # 指定文件夹路径  
  folder_path <- paste(filename, "KEGG", sep = "/") 
  
  # 获取文件夹中的文件列表  
  file_list <- list.files(folder_path)  
  
  # 遍历文件列表，删除每个文件  
  for (file in file_list) {  
    file_path <- file.path(folder_path, file)  
    file.remove(file_path)  
  }
}
for (i in seq_along(gseKEGG_res@result$ID)) {
  cat("Plot ", gseKEGG_res@result$ID[i], " ", gseKEGG_res@result$Description[i], "\n", sep = "")
  p <- gseaplot2(gseKEGG_res, geneSetID = gseKEGG_res@result$ID[i], title = gseKEGG_res@result$Description[i])
  if(drawPdf){
    pdf(file = paste(filename, "/KEGG/", gseKEGG_res@result$ID[i], ".pdf", sep = ""), width = 6, height = 4.5)
    print(p)
    dev.off()
  }
  png(file = paste(filename, "/KEGG/", gseKEGG_res@result$ID[i], ".png", sep = ""), width = 6, height = 4.5, units = "in", res = 500)
  print(p)
  dev.off()
}

