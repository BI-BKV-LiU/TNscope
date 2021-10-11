SHELL := bash
.ONESHELL:
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

GNUP_CONDA_ENV_NAME = gnup
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
GNUP_CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(GNUP_CONDA_ENV_NAME)

JOIN_SCRIPT = join_tables.sh
#RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/MergedProbe_ROstergotland_Onco_v2_TE-94002956_hg19
RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/MergedProbe_ROstergotland_Onco_v2_TE-94002956_hg19.padded20.bed_with_duplicates_included
#RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/pool1_pool2_nochr_3c.sort.merged.hg19.210311.met.annotated.bed_with_existant_ref_dict_used
#RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/coverage_cap_200
OUT_DIR = metrics

.PHONY: clean, all, run, fetch_sample_coverage_files

.DEFAULT_TARGET: all

all: clean run fetch_sample_coverage_files

run:
	mkdir -p $(OUT_DIR) ; \
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*_hs_metrics.txt" hs_metrics_all_samples.tsv $(OUT_DIR) 1 BAIT_SET
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*.aln.summary.metrics.txt" aln_metrics_all_samples.tsv $(OUT_DIR) 3 CATEGORY
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*.gcbias.summary_metrics.txt" gc_bias_metrics_all_samples.tsv $(OUT_DIR) 2 READS_USED
	bash $(JOIN_SCRIPT) $(RESULTS_ROOT_DIR) "*.insert_size_metrics.txt" insert_size_metrics_all_samples.tsv $(OUT_DIR) 1 MEDIAN_INSERT_SIZE

fetch_sample_coverage_files:
	$(GNUP_CONDA_ACTIVATE) ; \
	mkdir -p $(OUT_DIR) ; \
	find $(RESULTS_ROOT_DIR) -name "PVAL_*_per_target_cov.tsv" -type f | parallel cp {} $(OUT_DIR)

clean:
	rm -rf $(OUT_DIR)
