#!/bin/bash
set -e
set -u
set -o pipefail

SEARCH_DIR=$1
SEARCH_TERM=$2
# Output file in which everything will be joined
ALL_SAMPLES=$3

# Remove the output file if it existed since before
rm -f $ALL_SAMPLES

# Create array of file names
read -r -a dir <<< $( find $SEARCH_DIR -name "$SEARCH_TERM" -and -type f -print0 | xargs -0 echo )
# Remove empty lines and move the sequences to 1 row
for f in "${dir[@]}"; do
	sample_name=$(basename -s ".tsv" "$f" | cut -d _ -f2-4)
	echo "Processing sample: $sample_name"
	# Add sample name as first elements in the tsv files
   	sed -i "s/^/$sample_name\t/" $f
	cat $f >> $ALL_SAMPLES
done
