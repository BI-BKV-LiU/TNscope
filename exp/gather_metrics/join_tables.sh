#!/bin/bash

# http://linuxcommand.org/lc3_man_pages/seth.html
#set -x
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -e
set -uo pipefail

INPUT_FILE_DIR=$1
FIND_PATTERNS=$2
JOINED_TBL=$3
OUT_DIR=$4
NUM_LINES=$5
DATA_HEADER_PATTERN=$6



mkdir -p $OUT_DIR
OUT_FILE="$OUT_DIR""/""$JOINED_TBL"
# Make sure that there is no out file before running anything
rm -f "$OUT_FILE"


start_time=$(date +"%c")
echo "Start time: $start_time"
start=$(date +%s)

input=$( find $INPUT_FILE_DIR -name "$FIND_PATTERNS" -and -type f -print0 | xargs -0 echo )
read -r -a input_files_array <<< $input

# This suffix is needed for extracting the sample id/name
SUFFIX=$(echo $FIND_PATTERNS | cut -c2-) # Remove "*" from the beginning of the string

# Loop through all the samples
for i in "${!input_files_array[@]}"; do
	FULL_FILE_PATH=$(echo ${input_files_array[$i]})
	SAMPLE_NAME=$(basename -s "$SUFFIX" $(echo $FULL_FILE_PATH))
	# Grab the Metrics data
	ALL=$(grep -A"$NUM_LINES" "$DATA_HEADER_PATTERN" "$FULL_FILE_PATH")
	# Grab header and data as they are
	HEADER=$(echo "$ALL" | head -n1 | sed 's/[[:blank:]]/\\t/g')
	DATA=$(echo "$ALL" | tail -n+2)
	# Augment with sample info
	HEADER=$(echo -e "SAMPLE\t")$(echo $HEADER)
	DATA=$(echo "$ALL" | tail -n+2 | awk -v OFS="\t" -v sample_name=$SAMPLE_NAME '{print sample_name,$0}')
	# Append it to the data file
	echo "$DATA" >> "$OUT_FILE"
done

# Remove blank lines from the final output file
sed -i "1 i $HEADER" "$OUT_FILE"
end=$(date +%s)
runtime=$((end-start))
echo "The run took: $runtime s"