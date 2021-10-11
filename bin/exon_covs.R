#!/usr/bin/env Rscript

'This script creates multiple barplots (one for each exon) of a given coverage sample file
It can be created using command similar to the following:
bedtools coverage -hist -a genes_of_interest.bed -b PVAL_65_S1.tumor_deduped.bam > PVAL_65_S1.exon_cov.tsv

Usage:
exon_covs.R  (-n <id>|--ncbi <id>) (-g <gname>|--gene_name <gname>) (-s <sname>|--sample_name <sname>) [-o <outdir>|--outdir <outdir>] <sfile>
exon_covs.R (-h | --help)
exon_covs.R (-v | --version)

Arguments:
sfile    Coverage sample tsv file, e.g. PVAL_65_S1.exon_cov.tsv

Options:
-h --help                             Show this screen.
-v --version                          Show version.
-n <id>, --ncbi <id>                  NCBI accession ID of the transcript of interest for which bedtools got coverage values, e.g. NM_000110.4.
-g <gname>, --gene_name <gname>       Common name/symbol of the gene of interest, e.g. DPYD.
-s <sname>, --sample_name <sname>     Name of the sample, e.g. PVAL_65_S1.
-o <outdir>, --outdir <outdir>        Output directory with respect to the current working directory [default: temp/].
' -> doc


# Load the libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(here)
  library(docopt)
})

args <- docopt(doc, version = 'exon_covs 1.0\n')
print(args)

# https://here.r-lib.org/
source(here("bin", "utils.R"))

main <- function(cov_sample_input, NCBI_ID, gene_symbol, sample_name, out_results_dir) {  
  # Get rid off summary coverage values for all transcripts 
  cov_sample_input %<>% filter(chrom != "all") %>% 
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
  # Extract one user specified transcript into a list and choose only col:s of interest
  exons <- cov_sample_input %>% 
    filter(ID == NCBI_ID) %>% 
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
  
  transcript_name <- paste(gene_symbol, NCBI_ID)
  
  p <- wrap_plots(all_plots) + 
    plot_annotation(title = paste0("Gene: ", gene_symbol, ", sample: ", sample_name))
  
  out_file_basename <- here(out_results_dir, paste0(gene_symbol,"_",NCBI_ID))
  ggsave(filename = paste0(out_file_basename,".svg"),
         plot = p, 
         device = "svg",
         width = 25,
         height = 25,
         units = "in"
  )
  
  # Write session info for reproducibility purposes
  writeLines(capture.output(sessionInfo()), paste0(out_file_basename,"_","sessionInfo.txt"))
}

cov_sample_input <- args$sfile # E.g. # results_PVAL_65_S1.tsv or stdin
NCBI_name <- args$n # E.g. "NM_000110.4"
common_name <- args$g # E.g. "DPYD"
sample_name <- args$s # E.g. "PVAL_65_S1"
out_dir <- test_output(here(getwd(),args$o)) # E.g "temp/"

input <- read_tsv(cov_sample_input, 
                  col_names = c("chrom",
                                "start",
                                "end",
                                "name",
                                "score",
                                "strand",
                                "depth",
                                "num_bases_at_depth",
                                "size_of_feature",
                                "pros_of_feature_at_depth"),
                  col_types = "ciiciciiid")

class(out_dir)
print(out_dir)
main(cov_sample_input = input,
     NCBI_ID = NCBI_name, 
     gene_symbol = common_name, 
     sample_name = sample_name, 
     out_results_dir = out_dir)
