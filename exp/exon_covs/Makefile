SHELL := bash
.ONESHELL:
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

RESULTS_SPREAD_OUT_DIR = /home/rada/Documents/TNscope/MergedProbe_ROstergotland_Onco_v2_TE-94002956_hg19_with_existant_ref_dict_used

GNUP_CONDA_ENV_NAME = gnup
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
GNUP_CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(GNUP_CONDA_ENV_NAME)

EXON_COVS_CONDA_ENV_NAME = exon_covs
# # Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
EXON_COVS_CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(EXON_COVS_CONDA_ENV_NAME)

#R_CONDA_ENV_NAME = rstudio-renv
#Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
#R_CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(R_CONDA_ENV_NAME)

ALL_EXON_COV_FILES_DIR = cov
OUT_FILES_DIR = temp

ARGS = --ncbi NM_000110.4 \
--gene_name DPYD \
--sample_name PVAL_65_S1 \
--outdir $(OUT_FILES_DIR)/ \
$(ALL_EXON_COV_FILES_DIR)/PVAL_65_S1.exon_cov.tsv

.PHONY: clean, all, fetch_sample_coverage_files, create_individual_exon_cov_plot

#.DEFAULT_TARGET: all

all: clean fetch_sample_coverage_files create_individual_exon_cov_plot

fetch_sample_coverage_files:
	$(GNUP_CONDA_ACTIVATE) ; \
	mkdir -p $(ALL_EXON_COV_FILES_DIR) ; \
	find $(RESULTS_SPREAD_OUT_DIR) -name "PVAL_*.exon_cov.tsv" -type f | parallel cp {} $(ALL_EXON_COV_FILES_DIR)

create_individual_exon_cov_plot:
	$(EXON_COVS_CONDA_ACTIVATE) ; \
	Rscript exon_covs.R $(ARGS)

clean:
	rm -rf $(ALL_EXON_COV_FILES_DIR); \
	rm -rf $(OUT_FILES_DIR)

