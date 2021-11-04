SHELL = /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

CONDA_ENV_NAME = off-target-annotation
# Note that the extra "conda activate" is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CONDA_ENV_NAME)

INTERSECT_GTF = int_out.gtf
INTERSECT_BED = int_out.bed
NAME_NAME2 = names.uniq.tsv
ANNOT = annotated.tsv


DROP = annotated.dropped.tsv
UNIQ = unique.tsv

UNFLAGGED = unflagged.tsv
FLAGGED = flagged.bed
FEATURES = features.tsv
UNUNIQUE = ununique_removed.tsv
REMOVED = features_removed.tsv
DEDUPED = deduped.tsv
DEDUP_ANNOT = target_bases_not_covered_ROstergotland_Onco_v2_TE-94002956_hg19_150XDS_cutoff_99_0.annotated.tsv

TWIST_99_OFF_TARGET = target_bases_not_covered_ROstergotland_Onco_v2_TE-94002956_hg19_150XDS_cutoff_99_0.bed

.PHONY: all, \
join_bed, \
join_tsv, \
take_names, \
join_names, \
innerjoin_names, \
drop_cols, \
take_unique_positions, \
check_missing, \
add_flags, \
check_annot_genes, \
check_GoI_regions, \
check_non_GoI_regions, \
join_gtf, \
remove_uninteresting_data, \
remove_ununique_regions, \
append_missing, \
check_how_many_regions, \
remove_cds_start_stop_codons, \
change_dots_to_introns, \
remove_duplicate_rows, \
check_that_we_have_same_regions, \
add_actions, \
clean, \
help

all: clean \
join_bed \
take_names \
join_names \
drop_cols \
take_unique_positions \
add_flags \
join_gtf \
remove_uninteresting_data \
remove_ununique_regions \
append_missing \
remove_cds_start_stop_codons \
change_dots_to_introns \
remove_duplicate_rows \
add_actions

## join_bed: Join UCSC table browser bed-format annotations to off-target effect regions (otr)
join_bed:
	$(CONDA_ACTIVATE) ; \
	bedtools intersect -a $(TWIST_99_OFF_TARGET) \
	-b off_target_regions_99.annot.bed \
	-wa \
	-wb > $(INTERSECT_BED)

## join_tsv: Join UCSC table browser tsv-format annotations to otr
join_tsv:
	$(CONDA_ACTIVATE) ; \
	bedtools intersect \
	-a $(TWIST_99_OFF_TARGET) \
	-b off_target_regions_99.annot.tsv \
	-wa \
	-wb > int_out.tsv

## take_names: Extract NCBI_id and NCBI_gene_name/symbol information from UCSC table browser tsv format annotations 
take_names:
	awk 'BEGIN{OFS="\t"}{print $$2, $$13}' off_target_regions_99.annot.tsv | \
	sort -u | \
	tail --lines=+2 > $(NAME_NAME2)

## join_names: Join NCBI_gene_name/symbol information to with UCSC table browser bed-format annotatated regions
join_names:
	$(CONDA_ACTIVATE) ; \
	csvtk join --left-join -t -f "7;1" $(INTERSECT_BED) $(NAME_NAME2) > $(ANNOT)

## innerjoin_names: Innerjoin NCBI_gene_name/symbol information to with UCSC table browser bed-format annotatated regions
innerjoin_names:
	$(CONDA_ACTIVATE) ; \
	csvtk join -t -f "7;1" $(INTERSECT_BED) $(NAME_NAME2) > inner_annotated.tsv

## drop_cols: Drop uninteresting bed-format columns from annotated regions file
drop_cols:
	cut -d $$'\t' -f1-3,7,16 $(ANNOT) > $(DROP)

## take_unique_positions: Drop regions that are ununique
take_unique_positions:
	$(CONDA_ACTIVATE) ; \
	gawk -f summarise.awk $(DROP) > $(UNIQ) ; \
	sort $(UNIQ) > $(UNFLAGGED)

## check_missing: Check what regions weren't annotated with gene ID
check_missing:
	comm -23 <(sort $(TWIST_99_OFF_TARGET)) <(cut --complement -d $$'\t' -f4 $(UNIQ) | sort)

## add_flags: Add "GoI" flag to all genes that were in a genes of interest list
add_flags:
	$(CONDA_ACTIVATE) ; \
	gawk -f flag_GoIs.awk $(UNFLAGGED) > $(FLAGGED)

## check_annot_genes: Check which of the genes of interest were found in the annotations
check_annot_genes:
	grep "GoI" $(FLAGGED) | cut -d $$'\t' -f4 | sort -u | cat -n

## check_GoI_regions: Check the GoI regions
check_GoI_regions:
	grep "GoI" $(FLAGGED) | cat -n

## check_non_GoI_regions: Check the regions which weren't GoIs
check_non_GoI_regions:
	grep -v "GoI" $(FLAGGED) | cat -n

## join_gtf: Add gtf annotations data to regions with flagged gene names
join_gtf:
	$(CONDA_ACTIVATE) ; \
	bedtools intersect -a $(FLAGGED) \
	-b  off_target_regions_99.annot.gtf \
	-loj > $(INTERSECT_GTF)

## remove_uninteresting_data: Drop uninteresting gtf-format annotation columns
remove_uninteresting_data:
	cut -d $$'\t' -f1-4,7-9,13 $(INTERSECT_GTF) > $(FEATURES)

## remove_ununique_regions: Remove regions that differ only on NCBI gene name or transcript name
remove_ununique_regions:
	cut -d $$'\t' -f1-7 $(FEATURES) | sort -u > $(UNUNIQUE)

## append_missing: Append regions that weren't in any genes to the end
append_missing:
	comm -23 <(sort $(TWIST_99_OFF_TARGET)) <(cut --complement -d $$'\t' -f4 $(UNIQ) | sort) >> $(UNUNIQUE)

## check_how_many_regions: Do a check on how many regions is left
check_how_many_regions:
	cut -d $$'\t' -f1-3 $(UNUNIQUE) | sort -u | wc -l

remove_cds_start_stop_codons:
	grep -v CDS $(UNUNIQUE) | grep -v stop_codon | grep -v start_codon > $(REMOVED)

change_dots_to_introns:
	sed -i -e 's/\t\.\t/\tintron\t/' -e 's/\t-1\t-1/\tNA\tNA/g' $(REMOVED)

remove_duplicate_rows:
	$(CONDA_ACTIVATE) ; \
	gawk -f remove_duplicate_rows.awk $(REMOVED) > $(DEDUPED)

check_that_we_have_same_regions:
	comm -12 <(sort $(TWIST_99_OFF_TARGET)) <(cut -d $$'\t' -f1-3 $(DEDUPED) | sort) | cat -n

## help: Show this message
help:
	@grep '^##' ./Makefile
