SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: clean, all, gather

.DEFAULT_TARGET: all

# Dirs to search and not search PVAL_*_hs_metrics.txt files
RESULTS_DIR = /home/rada/Documents/TNscope
EXCLUDE_DIR1 = /home/rada/Documents/TNscope/temp
EXCLUDE_DIR2 = /home/rada/Documents/TNscope/exp

# Where do all the PVAL_*_hs_metrics.txt end up
GATHERED_DIR = /home/rada/Documents/TNscope/exp/gather_Hs_metrics/deduped_hs_metrics

# Find paths to PVAL_*_hs_metrics.txt files in the RESULTS_DIR
METRICS_FP = $(shell find $(RESULTS_DIR) -type d \( -path $(EXCLUDE_DIR1) -o -path $(EXCLUDE_DIR2) \) -prune -false -o -type f -name '*_hs_metrics.txt')
#METRICS = $(notdir $(METRICS_FP))

# The file where all the single PVAL_*_hs_metrics.txt files are gathered into
JOINED = hsmetrics_all_samples.tsv

JOIN_SCRIPT = /home/rada/Documents/TNscope/exp/gather_Hs_metrics/join_tables.sh

all: clean gather $(JOINED)

# Copy PVAL_*_hs_metrics.txt files to GATHERED_DIR
$(GATHERED_DIR):$(METRICS_FP)
	mkdir -p $@
	cp $^ $@

gather: $(GATHERED_DIR)

# Join all gathered PVAL_*_hs_metrics.txt files into one tsv
$(JOINED): $(JOIN_SCRIPT) $(GATHERED_DIR) 
	bash $^ $@

clean:
	rm -rf $(GATHERED_DIR)
