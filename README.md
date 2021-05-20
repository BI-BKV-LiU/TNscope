# TNscope
Sentieon TNscope basic pipeline

__Authors:__ Massimiliano Volpe, Jyotirmoy Das\
__Email:__ _massimiliano.volpe@liu.se_, _jyotirmoy.das@liu.se_\
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


  Download from [Google Cloud](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false 'Link to the resource') following reference files and extract them to `references/`:

  - hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf.gz
  - hg19_v0_1000G_omni2.5.b37.vcf.gz
  - hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf


