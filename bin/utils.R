#!/usr/bin/env Rscript

# From https://stackoverflow.com/a/26153325
# or https://stackoverflow.com/a/15785789

open_read <- function(arg) {
  if (arg %in% c("-", "/dev/stdin")) {
    file("stdin", open = "r")
  } else if (grepl("^/dev/fd/", arg)) {
    fifo(arg, open = "r")
  } else {
    file(arg, open = "r")
  }
}

# This test function was originally licensed under MIT by pyrevo and
# obtained from https://github.com/pyrevo/degs2venn/blob/main/fun.R

test_output <- function(out) {
  if (dir.exists(out)) {
    write("Output directory already exists!", stdout())
    return(out)
  } else {
    dir.create(out, recursive=TRUE)
    return(out)
  }
}

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