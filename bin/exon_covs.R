#!/usr/bin/env Rscript
'This script creates multiple barplots (one for each exon) of a given coverage sample file
The coverage sample file is created using similar to following command:
bedtools coverage -hist -a genes_of_interest.bed -b PVAL_68_S2.tumor_deduped.bam > covs_PVAL_65_S1.tsv

Usage:
  exon_covs.R ship new <name>...
  exon_covs.R ship <name> move <x> <y> [--speed=<kn>]
  exon_covs.R ship shoot <x> <y>
  exon_covs.R mine (set|remove) <x> <y> [--moored | --drifting]
  exon_covs.R (-h | --help)
  exon_covs.R --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  -s --sample=<kn>  Sample file name, e.g. covs_PVAL_65_S1.tsv.
  --moored      Moored (anchored) mine.
  --drifting    Drifting mine.

' -> doc

library(docopt)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)

arguments <- docopt(doc, version = 'Exon coverages plot 1.0')
print(arguments)

# Input arguments:
# 1. Sample file name, e.g. covs_PVAL_65_S1.tsv
# 2. NCBI accession ID of the transcript of interest for which bedtools got coverage values, e.g. NM_000110.4
# 3. Common name/symbol of the gene of interest, e.g. DPYD
# 4. Name of the sample, e.g. PVAL_65_S1
# 5. Output directory, e.g. temp/
#
#
# One example run could like the following:
# Rscript exp/create_barplot/exon_covs.R exp/create_barplot/cov/results_PVAL_65_S1.tsv NM_000110.4 DPYD PVAL_65_S1 temp/



#### Functions

plot_bar <- function(df, exon_no, length_of_exon) {
  # 
  plot_title <- paste("Exon number:", exon_no, "- length:", length_of_exon)
  ggplot(data=df, aes(x=depth, y=num_bases_at_depth)) +
    # ggplot(data=df, aes(x=sum_of_num_bases_at_depth, y=depth)) +
    geom_bar(stat="identity") + 
    xlab("Depth") + 
    ggtitle(plot_title) +
    ylab("Number of bases")
}

#### 

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  in_file <- args[1] # E.g. # results_PVAL_65_S1.tsv
  NCBI_name <- args[2] # E.g. "NM_000110.4"
  common_name <- args[3] # E.g. "DPYD"
  sample_name <- args[4] # E.g. "PVAL_65_S1"
  out_dir <- args[5] # E.g "temp/"

  input <- read_tsv(in_file, 
                    col_names = c("chrom",
                                  "start",
                                  "end",
                                  "name",
                                  "score",
                                  "strand",
                                  "depth",
                                  "num_bases_at_depth",
                                  "size_of_feature",
                                  "pros_of_feature_at_depth"))
  
  # Get rid off summary coverage values for all transcripts 
  input %<>% filter(chrom != "all")  %>% 
    # Parse out "_" concatenated transcript info to separate columns
    separate(col = "name",
             into = c("ID","rest"),
             sep = "_cds_") %>% 
    separate(col = "rest",
             into = c("exon_number",
                      "unknown",
                      "exon_chrom",
                      "exon_start_pos",
                      "exon_strand"),
             sep = "_")
  # Extract 1 user specified transcript into a list and choose only col:s of interest
  exons <- input %>% filter(ID == NCBI_name) %>% 
    select(one_of(c("exon_number",
                    "depth",
                    "num_bases_at_depth",
                    "size_of_feature",
                    "pros_of_feature_at_depth"))) %>% 
    group_by(exon_number) %>% 
    group_split()

  # Create a list of coverage plots for each exon in the transcript  
  all_plots <- vector(mode = "list", length = 0)
  i <- 0
  for(df in exons){
    i <- i + 1
    length_of_exon <- df[1,4][[1]]
    exon_no <- df[1,1][[1]]
    p <- df %>% plot_bar(exon_no, length_of_exon)
    all_plots[[i]] <- p
  }

  transcript_name <- paste(common_name, NCBI_name)

  p <- wrap_plots(all_plots) + 
    plot_annotation(title = paste0("Gene: ", common_name, ", sample: ", sample_name))
  
  out_file_base_name <- paste0(out_dir, common_name,"_",NCBI_name)

  ggsave(filename = paste0(out_file_base_name,".svg"), 
         plot = p, 
         device = "svg",
         width = 25,
         height = 25,
         units = "in"
  )

  # Write session info for reproducibility purposes
  writeLines(capture.output(sessionInfo()), paste0(out_file_base_name,"_","sessionInfo.txt"))
}

main()


