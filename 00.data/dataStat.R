#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program to summarize fastp json")
## Add command line arguments
p <- add_argument(parser = p, arg = "--dir", short = "-d", help = "path to folder containing fastp json result", type = "character")
p <- add_argument(parser = p, arg = "--out", short = "-o", help = "output file prefix", type = "character")
# Parse the command line arguments
argv <- parse_args(p)

dir <- argv$dir
out <- argv$out

library(jsonlite)
library(tidyverse)

currwd <- getwd()
setwd(dir = dir)

all <- tibble(sample=character(0), 
              total_reads.raw=integer(0),   total_bases.raw=integer(0),   q20_bases.raw=integer(0),   q30_bases.raw=integer(0),   q20_rate.raw=numeric(0),   q30_rate.raw=numeric(0),   read1_mean_length.raw=numeric(0),   read2_mean_length.raw=numeric(0),   gc_content.raw=numeric(0), 
              total_reads.clean=integer(0), total_bases.clean=integer(0), q20_bases.clean=integer(0), q30_bases.clean=integer(0), q20_rate.clean=numeric(0), q30_rate.clean=numeric(0), read1_mean_length.clean=numeric(0), read2_mean_length.clean=numeric(0), gc_content.clean=numeric(0), 
              effective.rate=numeric(0))
#sample <- "DiParent"
for (sample in list.files(pattern = "json") %>% basename() %>% str_remove(".json")) {
  cat("read ", sample, ".json ...\n", sep = "")
  df <- jsonlite::fromJSON(paste(sample, "json", sep = "."))
  raw <- as_tibble(df$summary$before_filtering) %>% mutate(sample = sample)
  clean <- as_tibble(df$summary$after_filtering) %>% mutate(sample = sample)
  
  stat <- tibble(sample = sample) %>% 
    left_join(raw, by = "sample") %>%
    left_join(clean, by = "sample", suffix = c(".raw", ".clean"))
  
  all <- plyr::rbind.fill(all, stat)
}

all <- all %>% mutate(effective.rate = total_bases.clean / total_bases.raw)
setwd(dir = currwd)
write_csv(x = all, file = paste(out, "csv", sep = "."))
write_tsv(x = all, file = paste(out, "txt", sep = "."))
