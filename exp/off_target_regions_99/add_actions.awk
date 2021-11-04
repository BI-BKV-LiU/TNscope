#!/usr/bin/awk -f

function append_action(action_name, action_id){
	print $0, action_name, action_id
}

BEGIN{
    OFS="\t";
	print "chr","start","end","id","feature","feature_start","feature_end","action","action_id"
}
{
	chr = $1
	start = $2
	end = $3
	id = $4
	feature = $5
	feature_start = $6
	feature_end = $7

	## Rules for actions ##
	
	# region					action_id		action
	# not in a gene:			1				remove, if okay for cnv backbone?
	# exon + BRCA: 				2				relocate
	# intron + BRCA: 			3				remove
	# exon + nonBRCA-GoI: 		4				keep
	# intron + nonBRCA-GoI: 	5				remove, if okay for cnv backbone?
	# exon + not GoI: 			6				remove, purpose of design?
	# intron + not GoI: 		7				remove, purpose of design?

	# not in a gene: 1 remove, if okay for cnv backbone?
	if(length(id) == 0 && length(feature) == 0){
		print $0,"NA","NA","NA","NA","Remove, if okay for cnv backbone?","1"
	
	# exon + BRCA: 2 relocate
	}else if(feature ~ /exon/ && id ~ /BRCA._GoI/ ){
		append_action("Relocate", 2)
	
	# intron + BRCA: 3 remove
	}else if(feature ~ /intron/ && id ~ /BRCA._GoI/){
		append_action("Remove", 3)
	
	}else if(!/BRCA./ && id ~ /_GoI/){
	# exon + nonBRCA-GoI: 4	keep
		if(feature ~ /exon/){
			append_action("Keep", 4)	
	# intron + nonBRCA-GoI: 5 remove, if okay for cnv backbone?
		}else if(feature ~ /intron/){
			append_action("Remove, if okay for cnv backbone?", 5)	
		# Warning if the script is not matching anything
		}else{
			print $0, "Not a match anywhere: GoI", -1
		}

	}else if(!/_GoI/){
		# exon + not GoI: 6	remove, purpose of design?
		if(feature ~ /exon/){
			append_action("Remove, purpose of design?", 6)
		# intron + not GoI: 7 remove, purpose of design?
		}else if(feature ~ /intron/){
			append_action("Remove, purpose of design?", 7)
		# Warning if the script is not matching anything
		}else{
			print $0, "Not a match anywhere: Not GoI", -2
		}
	
	# Warning if the script is not matching anything
	}else{
		print $0, "Not a match anywhere", -3
	}
}
