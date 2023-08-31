# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("UpSet plot")

# Add command line arguments]
p <- add_argument(p, "--de_log2FoldChange", help="log2FoldChange cutoff", type="numeric", default = 1)
p <- add_argument(p, "--de_padj", help="adjust pvalue cutoff", type="numeric", default = 0.05)
#p <- add_argument(p, "--upsetH", help = "the height of upset plot", type = "numeric", default = 5)
#p <- add_argument(p, "--upsetW", help = "the width of upset plot", type = "numeric", default = 7)

# Parse the command line arguments
argv <- parse_args(p)

de_log2FC <- argv$de_log2FoldChange
de_padj <- argv$de_padj
#width <- argv$upsetW
#height <- argv$upsetH

test <- FALSE
if (test) {
  de_log2FC <- 1
  de_padj <- 0.05
#  width <- 7
#  height <- 5
}

library(UpSetR)
library(tidyverse)

#setwd("..")
#dir <- list.dirs(path = ".")
#setwd(dir[str_detect(dir, "DESeq2")][1])

files <- list.files(path = "./", pattern = "^genes.+DE_results$")
if (length(files)>=2) {
  item <- files %>% str_remove("genes.counts.matrix.") %>% str_remove(".DESeq2.DE_results")
  
  dataforUpset <- list()
  for (i in 1:length(item)) {
    dataforUpset[[item[i]]] <- read.table(files[i]) %>% 
      filter(abs(log2FoldChange) > de_log2FC & padj < de_padj) %>% 
      rownames(id)
  }
  #p <- upset(fromList(dataforUpset))
  
  pdf(file = "upsetPlot.pdf", height = length(item) * 0.5 + 3.5, width = length(item) * 0.5 + 6)
  #print(p)
  upset(fromList(dataforUpset), nsets = length(item), order.by = "freq")
  dev.off()
  
  png(filename = "upsetPlot.png", height = length(item) * 0.5 + 3.5, width = length(item) * 0.5 + 6, units = "in", res = 500)
  #print(p)
  upset(fromList(dataforUpset), nsets = length(item), order.by = "freq")
  dev.off()
}

