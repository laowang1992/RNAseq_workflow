#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Run differential expression analysis using DESeq2.")

# Add command line arguments
p <- add_argument(p, "--matrix", help="matrix of raw read counts (not normalized!)", type="character")
p <- add_argument(p, "--samples_file", help="tab-delimited text file indicating biological replicate relationships.", type="character")
p <- add_argument(p, "--min_reps", help="At least min count of replicates must have cpm values > min cpm value.", type="numeric", default = 2)
p <- add_argument(p, "--min_cpm", help="At least min count of replicates must have cpm values > min cpm value.", type="numeric", default = 1)
p <- add_argument(p, "--contrasts", help="file (tab-delimited) containing the pairs of sample comparisons to perform.", type="character")
p <- add_argument(p, "--shrinkage", help = "whether use shrunken LFC or not (we use \"ashr\" method here).", type = "logical", default = TRUE)

# Parse the command line arguments
argv <- parse_args(p)

matrix <- argv$matrix
samples_file <- argv$samples_file
min_reps <- argv$min_reps
min_cpm <- argv$min_cpm
contrasts <- argv$contrasts
shrinkage <- argv$shrinkage
# test
if (FALSE) {
  matrix <- "../03.Merge_result/genes.counts.matrix"
  samples_file <- "../00.data/samples.txt"
  min_reps <- 2
  min_cpm <- 1
  contrasts <- "./contrasts.txt"
  shrinkage <- TRUE
}

library(DESeq2)
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
  conditions <- data.frame(conditions = factor(c(rep(treatment, length(which(samples$group %in% treatment))), 
                                                 rep(control,   length(which(samples$group %in% control))))))
  rownames(conditions) <- colnames(rnaseqMatrix)
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
  dds <- DESeq(ddsFullCountTable)
  contrast <- c("conditions", treatment, control)
  if (shrinkage) {
    res <- lfcShrink(dds = dds, contrast = contrast, type = "ashr")
  } else {
    res <- results(dds, contrast)
  }
  baseMeanA <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == treatment])
  baseMeanB <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == control])
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
  res <- cbind(sampleA = treatment, sampleB = control, as.data.frame(res))
  res$padj[is.na(res$padj)] <- 1
  res <- as.data.frame(res[order(res$pvalue), ])
  write.table(x = res, file = paste(paste(treatment, "vs", control, sep = "_"), "DESeq2.DE_results.txt", sep = "."), sep = '\t', quote = FALSE)
  write.table(x = rnaseqMatrix, file = paste(paste(treatment, "vs", control, sep = "_"), "DESeq2.read_count.txt", sep = "."), sep = '\t', quote = FALSE)
}
