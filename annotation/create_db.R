library(argparser, quietly=TRUE)

p <- arg_parser("Parse a kegg json file")

# Add command line arguments
p <- add_argument(p, "--emapper", help = "emapper.annotations file", type = "character")
p <- add_argument(p, "--json", help = "kegg json file", type = "character")
p <- add_argument(p, "--organism", help = "organism", type = "character")
p <- add_argument(p, "--fasta", help = "pep fasta file", type = "character")

# Parse the command line arguments
argv <- parse_args(p)

library(jsonlite)
library(tidyverse)
library(seqinr)
library(clusterProfiler)

emapper <- argv$emapper
json <- argv$json
organism <- argv$organism
fasta <- argv$fasta

if (FALSE) {
  emapper <- "zs11.emapper.annotations.gz"
  json <- "bna00001.json"
  organism <- "bna"
  fasta <- "./zs11.longest_pep.fa"
}

source("./support_script.R")

# 解析json生成kegg_pathway
kegg_pathway <- parse_kegg_json(json = json, organism = organism)

# 读取emapper注释信息
emapper <- read_tsv(emapper, comment = "##", 
                    col_names = TRUE, trim_ws = TRUE) %>% 
  dplyr::select(GID = "#query", Gene_Symbol = "Preferred_name", GO = "GOs", 
                Ko = "KEGG_ko", Pathway = "KEGG_Pathway", 
                COG = "COG_category", Gene_Name = "Description")

# 生成gene2pathway
gene2ko <- emapper %>% dplyr::select(GID, Ko) %>%
  separate_rows(Ko, sep = ",", convert = F) %>% 
  filter(Ko != "-") %>% mutate(Ko = str_remove(Ko, "ko:"))

gene2pathway <- gene2ko %>% left_join(kegg_pathway, by = "Ko") %>%
  dplyr::select(GID, Ko, Pathway, Pathway_Name, Pathway_Class, Pathway_Subclass) %>%
  distinct() %>% filter(!is.na(Pathway))
# 生成pathway2gene
pathway2gene <- gene2pathway %>% dplyr::select(Pathway, GID) %>% dplyr::distinct()
pathway2name <- kegg_pathway %>% dplyr::select(Pathway, Pathway_Name) %>% dplyr::distinct()
ko2gene <- gene2ko %>% dplyr::select(Ko, GID) %>% dplyr::distinct()

save(ko2gene, pathway2gene, pathway2name, file = paste(organism, "kegg_info.RData", sep = "."))

# 生成gene2cog
cog_info <- read_delim(file = "./cog_funclass.tab", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
gene2cog <- emapper %>% dplyr::select(GID, COG) %>%
  separate_rows(COG, sep = "", convert = F) %>%
  filter(COG != "" & COG != "-") %>% distinct()
gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")


# 生成GO org.db
gene_info <- emapper %>% dplyr::select(GID, Gene_Name) %>%
  dplyr::filter(Gene_Name != "-") %>% dplyr::distinct()

gene2go <- emapper %>% dplyr::select(GID, GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  dplyr::filter(GO != "-") %>%
  mutate(EVIDENCE = 'IEA') %>% 
  dplyr::distinct()

create_orgdb(gene_info = gene_info, gene2go = gene2go, genus = "Brassica", species = "napus")

if (FALSE) {
  dir.create("R_Library")
  install.packages("org.Bnapus.eg.db_1.0.tar.gz", repos = NULL, lib = "R_Library")
  
  library("org.Bnapus.eg.db", character.only = TRUE, lib.loc = "R_Library")
  
  # 统计
  all_gene <- getName.list(read.fasta(file = fasta, 
                                      seqtype = 'AA'))
  annoStat(all_gene = all_gene, gene_info = gene_info, gene2go = gene2go, gene2cog = gene2cog, gene2pathway = gene2pathway)
}
