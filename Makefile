SHELL = /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

CURRENT_CONDA_ENV_NAME = TNscope
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

# Sample name for transcript exon coverage barplots
IN_FILE = exp/create_barplot/cov/results_PVAL_65_S1.tsv
NCBI_NAME = NM_000110.4
COMMON_NAME = DPYD
SAMPLE_NAME = PVAL_65_S1
OUT_DIR = res/

CPUS = 90

.PHONY: all, merge_targets, multiqc

all:
	@($(CONDA_ACTIVATE) ; \
	snakemake --cores $(CPUS) --config cpus=$(CPUS) $(ARGS))

# On why LC_ALL=C and the following command are run: 
# https://www.biostars.org/p/177653/#177764
# https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
merge_targets:
	$(CONDA_ACTIVATE) ; \
	cd references/ROstergotland_Onco_v2_TE-94002956_hg19/temp ; \
	LC_ALL=C ; \
	sort -t$$'\t' -k1,1 -k2,2n UCSC_combined_ROstergotland_Onco_v2_TE-94002956_hg19.bed > UCSC_combined_ROstergotland_Onco_v2_TE-94002956_hg19.sorted.bed
	bedtools merge -i UCSC_combined_ROstergotland_Onco_v2_TE-94002956_hg19.sorted.bed > UCSC_combined_ROstergotland_Onco_v2_TE-94002956_hg19.sorted.merged.bed

multiqc:
	$(CONDA_ACTIVATE) ; \
	multiqc . -f --ignore temp --ignore fastq-temp --ignore exp --ignore .snakemake

report:
	$(CONDA_ACTIVATE) ; \
	snakemake --report report.html

## exon_covs: Make barplots of all exons in a sample
exon_covs:
	$(CONDA_ACTIVATE) ; \
	Rscript bin/exon_covs.R $(IN_FILE) $(NCBI_NAME) $(COMMON_NAME) $(SAMPLE_NAME) $(OUT_DIR)