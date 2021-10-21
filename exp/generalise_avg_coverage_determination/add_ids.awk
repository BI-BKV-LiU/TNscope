#!/usr/bin/awk -f
BEGIN{
    OFS="\t";
	# Parse sample name out from file name
	split(ARGV[1],dir_file_array,".tumor_deduped");
	split(dir_file_array[1],raw_sample_array,"/");
	split(raw_sample_array[2],sample_array,"_");
	sample_name = sample_array[1]"_"sample_array[2]"_"sample_array[3]

	# Number of header lines
	header = 1;
	# What should be prepended on the header line
	header_row = "sample\tname\tchromosome\tstart\tend\tcoverage\tnum_bases_at_depth\tlength\tprocent_of_coverage_in_region"
}
{
	name = $2"-"$3
    if (NR <= header) { 
		printf "%s\n%s\t%s\t%s\n", header_row, sample_name, name, $0;
	}
	#region_number = NR - 1
	#name = region_number"_"$2"-"$3
	print sample_name, name, $0
}