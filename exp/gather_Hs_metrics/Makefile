SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: clean, all, gather#, join_tables

.DEFAULT_TARGET: all

RESULTS_DIR = /home/rada/Documents/TNscope
EXCLUDE_DIR1 = /home/rada/Documents/TNscope/temp
EXCLUDE_DIR2 = /home/rada/Documents/TNscope/exp
GATHERED_DIR = /home/rada/Documents/TNscope/exp/gather_Hs_metrics/deduped_hs_metrix
#METRICS_FP=$(wildcard $(RESULTS_DIR)/*_hs_metrics.txt)
METRICS_FP = $(shell find $(RESULTS_DIR) -type d \( -path $(EXCLUDE_DIR1) -o -path $(EXCLUDE_DIR2) \) -prune -false -o -type f -name '*_hs_metrics.txt')
METRICS = $(notdir $(METRICS_FP))
JOINED = hsmetrics_all_samples.tsv

JOIN_SCRIPT = /home/rada/Documents/TNscope/exp/gather_Hs_metrics/join_tables.sh

all: clean gather $(JOINED)

$(GATHERED_DIR):$(METRICS_FP)
	mkdir -p $@
	cp $^ $@


gather: $(GATHERED_DIR)

$(JOINED): $(JOIN_SCRIPT) $(GATHERED_DIR) 
	bash $^ $@

clean:
	rm -rf $(GATHERED_DIR)
