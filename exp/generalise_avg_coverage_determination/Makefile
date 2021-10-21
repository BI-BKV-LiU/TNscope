SHELL = /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

GEN_COV_CONDA_ENV_NAME = gen-cov
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
GEN_COV_CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(GEN_COV_CONDA_ENV_NAME)

JUPYTER_CONDA_ENV_NAME = jupyter-plotly
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
JUPYTER_CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(JUPYTER_CONDA_ENV_NAME)


ROI = merged_probe_file_ROstergotland_Onco_v2_IntronsLargerThan5kbp_CNV_0p048X_TE-94002956_hg19_merge0_210328173615.bed
BAMS = bams
RESULTS_ROOT_DIR = /home/rada/Documents/TNscope/MergedProbe_ROstergotland_Onco_v2_TE-94002956_hg19.bed_with_duplicates_included
SUMMARY_SAMPLE_COVS_DIR = all

MARK_DUP_SUFFIX = tumor_deduped.bam
RAW_COVS_SUFFIX = tumor_deduped.cov_chrs.tsv
IDED_COVS_SUFFIX = tumor_deduped.cov_chrs.ids.tsv


COVS_DIR = covs
CHRS_COVS = chrs
PLOTS_DIR = plots

NUM_JOBS = 8

.PHONY: link_bams, \
get_covs, \
help

# link_bams: create symbolic links to duplicate marked bam files
link_bams:
	$(GEN_COV_CONDA_ACTIVATE) ; \
	mkdir -p $(BAMS) ; \
	cd $(BAMS) ; \
	find $(RESULTS_ROOT_DIR) -name "PVAL_*.$(MARK_DUP_SUFFIX)" -type f | parallel ln -s {} $$PWD

## get_covs: produce coverage in regions of interest with bedtools coverage
get_covs:
	$(GEN_COV_CONDA_ACTIVATE) ; \
	mkdir -p $(COVS_DIR) ; \
	find -L $(BAMS)/ -name "PVAL_*.$(MARK_DUP_SUFFIX)" -type f | parallel -N1 --jobs $(NUM_JOBS) "bedtools coverage -hist -a $(ROI) -b {} > $(COVS_DIR)/{/.}.cov.tsv"

separate_all_from_chrs:
	$(GEN_COV_CONDA_ACTIVATE) ; \
	mkdir -p $(CHRS_COVS) ; \
	cd $(CHRS_COVS) ; \
	find ../$(COVS_DIR)/ -name "PVAL_*.tumor_deduped.cov.tsv" -type f | parallel -N1 --jobs $(NUM_JOBS) "awk -f ../separate_by_first_column.awk {}" ; \
	mkdir -p ../$(SUMMARY_SAMPLE_COVS_DIR) ; \
	mv all_* ../$(SUMMARY_SAMPLE_COVS_DIR)

add_ids_sample_names_and_header: clean separate_all_from_chrs
	$(GEN_COV_CONDA_ACTIVATE) ; \
	cd $(CHRS_COVS) ; \
	find . -name "PVAL_*.$(RAW_COVS_SUFFIX)" -type f | parallel -N1 --jobs $(NUM_JOBS) "awk -f ../add_ids.awk {} > {.}.ids.tsv" ; \
	rm -f *.$(RAW_COVS_SUFFIX)

create_cov_plots:
	$(GEN_COV_CONDA_ACTIVATE) ; \
	cd $(CHRS_COVS) ; \
	find . -name "PVAL_*.$(IDED_COVS_SUFFIX)" -type f | parallel -N1 --jobs $(NUM_JOBS) "python ../cov_vis.py {} {= s:.$(IDED_COVS_SUFFIX):: =}" ; \
	mkdir -p ../$(PLOTS_DIR) ; \
	mv PVAL_*.html ../$(PLOTS_DIR) ; \
	mv temp ..

clean:
	rm -rf $(CHRS_COVS) $(SUMMARY_SAMPLE_COVS_DIR) temp

start_jupyterlab:
	$(JUPYTER_CONDA_ACTIVATE) ; \
	jupyter lab

## help: show this message
help:
	@grep '^##' ./Makefile
