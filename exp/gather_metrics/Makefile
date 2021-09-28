SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

JOIN_SCRIPT = join_tables.sh
#RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/MergedProbe_ROstergotland_Onco_v2_TE-94002956_hg19
RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/pool1_pool2_nochr_3c.sort.merged.hg19.210311.met.annotated
#RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/coverage_cap_200
OUT_DIR = metrics

.PHONY: clean, all, run

.DEFAULT_TARGET: all

all: clean run

run:
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*_hs_metrics.txt" hs_metrics_all_samples.tsv $(OUT_DIR) 1 BAIT_SET
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*.aln.summary.metrics.txt" aln_metrics_all_samples.tsv $(OUT_DIR) 3 CATEGORY
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*.gcbias.summary_metrics.txt" gc_bias_metrics_all_samples.tsv $(OUT_DIR) 2 READS_USED
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*.insert_size_metrics.txt" insert_size_metrics_all_samples.tsv $(OUT_DIR) 1 MEDIAN_INSERT_SIZE

clean:
	rm -rf $(OUT_DIR)/hs_metrics_all_samples.tsv ; \
	rm -rf $(OUT_DIR)/aln_metrics_all_samples.tsv ; \
	rm -rf $(OUT_DIR)/gc_bias_metrics_all_samples.tsv ; \
	rm -rf $(OUT_DIR)/insert_size_metrics_all_samples.tsv
