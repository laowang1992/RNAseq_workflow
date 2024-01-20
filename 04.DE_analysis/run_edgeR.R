#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Run differential expression analysis using DESeq2.")

# Add command line arguments
p <- add_argument(p, "--matrix", help="matrix of raw read counts (not normalized!)", type="character")
p <- add_argument(p, "--samples_file", help="tab-delimited text file indicating biological replicate relationships.", type="character")
p <- add_argument(p, "--min_reps", help="At least min count of replicates must have cpm values > min cpm value.", type="numeric", default = 1)
p <- add_argument(p, "--min_cpm", help="At least min count of replicates must have cpm values > min cpm value.", type="numeric", default = 1)
p <- add_argument(p, "--contrasts", help="file (tab-delimited) containing the pairs of sample comparisons to perform.", type="character")
p <- add_argument(p, "--dispersion", help = "edgeR dispersion value.", type = "numeric", default = 0.1)

# Parse the command line arguments
argv <- parse_args(p)

matrix <- argv$matrix
samples_file <- argv$samples_file
min_reps <- argv$min_reps
min_cpm <- argv$min_cpm
contrasts <- argv$contrasts
dispersion <- argv$dispersion
# test
if (FALSE) {
  matrix <- "../03.Merge_result/genes.counts.matrix"
  samples_file <- "../00.data/samples.txt"
  min_reps <- 1
  min_cpm <- 1
  contrasts <- "./contrasts.txt"
  shrinkage <- TRUE
}

library(edgeR)

data <- read.table(file = matrix, header = TRUE, row.names = 1, com = '')
samples <- read.table(file = samples_file, header = FALSE, sep = "\t")
colnames(samples)[1:2] <- c("group", "sample")
contrasts <- read.table(file = contrasts, header = FALSE, sep = "\t")
colnames(contrasts)[1:2] <- c("treatment", "control")

for (i in 1:nrow(contrasts)) {
  treatment <- contrasts[i, 1][[1]]
  control <- contrasts[i, 2][[1]]
  col_ordering <- c(which(samples$group %in% treatment), which(samples$group %in% control))
  rnaseqMatrix <- data[,col_ordering]
  rnaseqMatrix <- round(rnaseqMatrix)
  rnaseqMatrix <- rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > min_cpm) >= min_reps, ]
  conditions <- factor(c(rep(treatment, 1), rep(control, 1)))
  
  exp_study <- DGEList(counts = rnaseqMatrix, group = conditions)
  exp_study <- calcNormFactors(exp_study)
  et <- exactTest(exp_study, pair = c(treatment, control), dispersion = 0.1)
  tTags <- topTags(et, n = NULL)
  result_table <- tTags$table
  result_table <- data.frame(sampleA = treatment, sampleB = control, result_table)
  result_table$logFC <- -1 * result_table$logFC
  
  write.table(result_table, file=paste0("genes.counts.matrix.", treatment, "_vs_", control, ".edgeR.DE_results.txt"), sep = "\t", quote = F, row.names = T)
  write.table(result_table, file=paste0("genes.counts.matrix.", treatment, "_vs_", control, ".edgeR.count_matrix.txt"), sep = "\t", quote = F, row.names = T)
}
