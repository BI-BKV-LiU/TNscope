#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  in_file <- args[1]
  cov_count <- args[2]
  baseq <- args[3]
  input <- read_tsv(in_file)
  p<-ggplot(data=input, aes(x=coverage_or_base_quality, y=high_quality_coverage_count)) +
    geom_bar(stat="identity")
  
  ggsave(cov_count, plot = p)
  
  p<-ggplot(data=input, aes(x=coverage_or_base_quality, y=unfiltered_baseq_count)) +
    geom_bar(stat="identity")
  
  ggsave(baseq, plot = p)
}

main()

cat(args, sep = "\n")





