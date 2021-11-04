#!/usr/bin/awk -f

BEGIN{
    OFS="\t";
}
{
	chr = $1
	start = $2
	end = $3
	id = $4
	feature = $5
	feature_start = $6
	feature_end = $7
	current_row = ""

	# Handle first row separately: print and save as previous_row
	if (NR <= 1){
		print $0
		previous_row = chr"_"start"-"end"_"id"_"feature
		next
	}
	current_row = chr"_"start"-"end"_"id"_"feature

	# Check if the previous row is the same as the current row
	if(previous_row == current_row){
		# Skip row
		previous_row = chr"_"start"-"end"_"id"_"feature
		next;
	}else{
		# Print row
		previous_row = chr"_"start"-"end"_"id"_"feature
		print $0
	}
}

