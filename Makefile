SHELL = /bin/bash

CURRENT_CONDA_ENV_NAME = TNscope
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

CPUS = 80

all:
	@($(CONDA_ACTIVATE) ; \
	snakemake --cores $(CPUS) --config cpus=$(CPUS))
