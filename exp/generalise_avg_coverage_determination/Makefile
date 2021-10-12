SHELL = /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

CURRENT_CONDA_ENV_NAME = gen-cov
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

ROI = merged_probe_file_ROstergotland_Onco_v2_IntronsLargerThan5kbp_CNV_0p048X_TE-94002956_hg19_merge0_210328173615.bed
BAMS = bams
RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/MergedProbe_ROstergotland_Onco_v2_TE-94002956_hg19.bed_with_duplicates_included

MARK_DUP_SUFFIX = tumor_deduped.bam

.PHONY: link_bams

link_bams:
	$(CONDA_ACTIVATE) ; \
	mkdir -p $(BAMS) ; \
	cd $(BAMS) ; \
	find $(RESULTS_ROOT_DIR) -name "PVAL_*.$(MARK_DUP_SUFFIX)" -type f | parallel ln -s {} $$PWD

