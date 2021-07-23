#!/bin/bash

# http://linuxcommand.org/lc3_man_pages/seth.html
#set -x
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -e
set -uo pipefail

RESULTS_DIR=$1
JOINED_TBL=$2

CWD=$(pwd)

# Clear log file for the total run
rm -f logs.log

# Grab tsv header from the first listed file in the RESULTS_DIR
HEADER=$(echo -e "SAMPLE\t")$(awk '/^BAIT_SET.+/{print $0}' "$RESULTS_DIR"/$(ls "$RESULTS_DIR" | head -n1))

rm -f $JOINED_TBL

input=$( find $RESULTS_DIR -name "*.txt" -and -type f -print0 | xargs -0 echo )
read -r -a input_files_array <<< $input
# Loop through all the samples
for i in "${!input_files_array[@]}"; do
  FULL_FILE_PATH=$(echo ${input_files_array[$i]})
  STEM=$(basename -s '.txt' $(echo $FULL_FILE_PATH))
  # Insert header on the first row of the file
  echo "$HEADER" >> "$JOINED_TBL"
  HEADER=''
  
  start_time=$(date +"%c")
  echo "Start time: $start_time - $STEM" >> logs.log
  start=$(date +%s)
  
  # Grab the HsMetrix data
  DATA=$(sed -n '8p' "$FULL_FILE_PATH")
  # Append it to the data file
  printf "$STEM\t$DATA\n" >> "$JOINED_TBL"
  end=$(date +%s)
  runtime=$((end-start))
  echo "The run took: $runtime s for $STEM" >> logs.log
  end_time=$(date +"%c")
  echo "End time: $end_time - $STEM" >> logs.log
  echo >> logs.log
done

# Remove blank lines from the final output file
sed -i '/^$/d' "$JOINED_TBL"