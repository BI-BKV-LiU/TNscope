#!/usr/bin/awk -f
BEGIN{
    OFS="\t";
	# Number of header lines
	header = 1;
	split(ARGV[1],file_name_array,"/"); # ../covs/PVAL_65_S1.tumor_deduped.cov.tsv
	sub(/\.tsv/,"_chrs.tsv",file_name_array[3]); # PVAL_65_S1.tumor_deduped.cov.tsv
	filename = file_name_array[3]
}
{
    if (NR <= header) { 
		next;
	}
	if ($0 ~ /^all/) {
		print $0 >> $1"_"filename
	}else{
		print $0 >> filename
	}
}
