# TNscope
Sentieon TNscope basic pipeline

__Authors:__ Massimiliano Volpe, Jyotirmoy Das, Lauri Mesilaakso (Adjustments in "ro-implementation"-branch)\
__Email:__ _massimiliano.volpe@liu.se_, _jyotirmoy.das@liu.se_, _lauri.mesilaakso@regionostergotland.se_\
__Date:__ 28/01/2021

__Developed on behalf of the Bioinformatics Core Facility, Linköping University__

## Rules:
- Raw data from matched tumor/normal samples should be stored in different folders but named the same, e.g.:\
  /path_to_normal/256_S1_R1_001.fastq.gz\
  /path_to_tumor/256_S1_R1_001.fastq.gz
- Paths to raw data and reference must be set in the config.json file.
- Fastq filename suffix must be set in the config.json file, e.g.:\
  /path_to_normal/256_S1_R1_001.fastq.gz --> "_R1_001.fastq.gz"\
  /path_to_normal/256_S1_R2_001.fastq.gz --> "_R2_001.fastq.gz"


Some reference files can also be downloaded from [Google Cloud](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false 'Link to the resource').

The pipeline uses [Sentieon v. 202010.02](https://support.sentieon.com/versions/202010.02/manual/).