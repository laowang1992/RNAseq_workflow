#!/usr/bin/env Rscript

# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("run featureCounts and calculate FPKM/TPM")

# Add command line arguments
p <- add_argument(p, "--bam", help="input: bam file", type="character")
p <- add_argument(p, "--gtf", help="input: gtf file", type="character")
p <- add_argument(p, "--output", help="output prefix", type="character")
p <- add_argument(p, "--nthread", help="thread number", type="integer", default="1")
p <- add_argument(p, "--attrType", help="calculate expression for gene_id or transcript_id, default: gene_id", type="character", default="gene_id")
p <- add_argument(p, "--strandSpecific", help="if strand-specific read counting should be performed, it should be one of the following three values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded), default: 0", type="integer", default="0")
p <- add_argument(p, "--featureType", help="the feature type used to select rows in the GTF annotation which will be used for read summarization", type="character", default="exon")


# Parse the command line arguments
argv <- parse_args(p)

library(Rsubread)
library(limma)
library(edgeR)

bamFile <- argv$bam
gtfFile <- argv$gtf
nthreads <- argv$nthread
outFilePref <- argv$output
attrType <- argv$attrType
strandSpecific <- argv$strandSpecific
featureType <- argv$featureType

outStatsFilePath  <- paste(outFilePref, '.log',  sep = ''); 
outCountsFilePath <- paste(outFilePref, '.count', sep = ''); 

fCountsList = featureCounts(bamFile, annot.ext=gtfFile, isGTFAnnotationFile=TRUE, GTF.featureType=featureType, nthreads=nthreads, isPairedEnd=TRUE, GTF.attrType=attrType, strandSpecific=strandSpecific)
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
fpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts, fpkm, tpm)
colnames(featureCounts) = c('gene_id', 'counts', 'fpkm','tpm')
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
