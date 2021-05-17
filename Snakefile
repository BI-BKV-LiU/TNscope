################################################################################
## TNscope Pipeline
## Somatic variant calling with Sentieon TNscope (from fastq files to vcf)
##
## Authors: Massimiliano Volpe, Jyotirmoy Das
## Email: massimiliano.volpe@liu.se, jyotirmoy.das@liu.se
## Date: 28/01/2021
## Developed on behalf of the Bioinformatics Core Facility, Linköping University
##
## Rules:
## - Raw data from matched tumor/normal samples should be stored in different
## folders but named the same, e.g.:
##      /path_to_normal/256_S1_R1_001.fastq.gz
##      /path_to_tumor/256_S1_R1_001.fastq.gz
## - Paths to raw data and reference must be set in the config.json file
## - Fastq filename suffix must be set in the config.json file, e.g.:
##      /path_to_normal/256_S1_R1_001.fastq.gz --> "_R1_001.fastq.gz"
##      /path_to_normal/256_S1_R2_001.fastq.gz --> "_R2_001.fastq.gz"
################################################################################


# Functions -------------------------------------------------------------------
def id_maker(path, sep, sample, suffix, d):
    f = "".join([path, sep, sample, suffix])
    # run:flowcell:lane
    # l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 2,3,4 | sed  s/:/./g | sed 's/@//'", shell=True).strip().decode()
    # flowcell:lane
    l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 3,4 | sed  s/:/./g | sed 's/@//'", shell=True).strip().decode()
    #d[sample] = '_'.join([sample, group, l])
    d[sample] = l
    return(d)


# Globals ---------------------------------------------------------------------

# Getting Sentieon license server running
# start licence server
import subprocess
# subprocess.call("/home/rada/miniconda3/pkgs/sentieon-201808.08-h14580f3_0/bin/sentieon licsrvr --start --log /home/rada/Documents/Sentieon/logs/log.txt /home/rada/Documents/Sentieon/licence/GMC_SOR_usb127.lic")
# Update with the location of the Sentieon software package and license file
# subprocess.call("export SENTIEON_INSTALL_DIR=/home/rada/miniconda3/pkgs/sentieon-201808.08-h14580f3_0")
# subprocess.call("export SENTIEON_LICENSE=Idril:8080")

configfile:
    "config.json"

#workdir:
#    config['workdir']


R1SUFFIX = config['R1_suffix']
R2SUFFIX = config['R2_suffix']

SAMPLES, = glob_wildcards(config['tumor'] + "/{sample}" + R1SUFFIX)
RESULTS = config['workdir'] + '/{sample}/'
BAMS = RESULTS + 'bams/'
LOGS = RESULTS + 'logs/'
METRICS = RESULTS + 'metrics/'
PLOTS = RESULTS + 'plots/'
MARKDUP = RESULTS + 'markdup/'
RECAL = RESULTS + 'baserecal/'

fasta = config['fasta']
dbsnp = config['dbsnp']
known_Mills_indels = config['known_Mills_indels']
known_1000G_indels = config['known_1000G_indels']


#for sample in SAMPLES:
# f = "".join([config['normal'], '/', sample, '_R1_001.fastq.gz'])
# l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 4 | sed  s/:/./ | sed 's/@//'", shell=True).strip().decode()
# d[sample] = '.'.join([sample, l])
n = {}
t = {}

for sample in SAMPLES:
    #id_maker(config['normal'], '/', sample, R1SUFFIX, n)
    id_maker(config['tumor'], '/', sample, R1SUFFIX, t)
#print(n)
#print(t)


# Rules -----------------------------------------------------------------------
rule all:
    input:
        #normal_bam = expand(BAMS + "{sample}.normal_sorted.bam", sample=SAMPLES),
        tumor_bam = expand(BAMS + "{sample}.tumor_sorted.bam", sample=SAMPLES),
        #normal_ml = expand(LOGS + '{sample}.normal_metrics.log', sample=SAMPLES),
        tumor_ml = expand(LOGS + '{sample}.tumor_metrics.log', sample=SAMPLES),
        #normal_deduped = expand(BAMS + '{sample}.normal_deduped.bam', sample=SAMPLES),
        tumor_deduped = expand(BAMS + '{sample}.tumor_deduped.bam', sample=SAMPLES),
        #normal_rdt = expand(RECAL + '{sample}.normal_recal_data.table', sample=SAMPLES),
        tumor_rdt = expand(RECAL + '{sample}.tumor_recal_data.table', sample=SAMPLES),
        vcf = expand(RESULTS + '{sample}.tnscope.vcf.gz', sample=SAMPLES)

'''
rule mapping_normal:
    input:
        R1 = config['normal'] + "/{sample}" + R1SUFFIX,
        R2 = config['normal'] + "/{sample}" + R2SUFFIX,
        fasta = config['fasta']
    output:
        sam = temp(BAMS + '{sample}.normal.sam'),
        bam = BAMS + '{sample}.normal_sorted.bam'
    log:
        bwa = LOGS + '{sample}.normal_bwa.log',
        sort = LOGS + '{sample}.normal_sort.log'
    params:
        #ID = subprocess.check_output("zcat {input.R1} | head -n 1 | cut -d ':' -f 1,4 | sed  s/:/./ | sed 's/@//'", shell=True).strip().decode(),
        #R = "@RG\\tID:" + config["normal_group"] + "\\tSM:{sample}_normal\\tPL:" + config["platform"]
        K = 10000000,
        ID = lambda wildcards: n[wildcards.sample],
        SM = "{sample}_normal",
        PL = config["platform"]
    threads:
        48 # set the maximum number of available cores
    shell:
        # sentieon bwa mem -M -R '{params.R}' -t {threads} -K {params.K} -o {output.sam} {input.fasta} {input.R1} {input.R2} >> {log.bwa} 2>&1
        """
        sentieon bwa mem -M -R '@RG\\tID:{params.ID}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam} {input.fasta} {input.R1} {input.R2} >> {log.bwa} 2>&1
        sentieon util sort -r {input.fasta} -i {output.sam} -o {output.bam} -t {threads} --sam2bam >> {log.sort} 2>&1
        """
'''

rule mapping_tumor:
    input:
        R1 = config['tumor'] + "/{sample}" + R1SUFFIX,
        R2 = config['tumor'] + "/{sample}" + R2SUFFIX,
        fasta = config['fasta']
    output:
        sam = temp(BAMS + '{sample}.tumor.sam'),
        bam = BAMS + '{sample}.tumor_sorted.bam'
    log:
        bwa = LOGS + '{sample}.tumor_bwa.log',
        sort = LOGS + '{sample}.tumor_sort.log'
    params:
        #R = "@RG\\tID:" + config["tumor_group"] + "\\tSM:{sample}_tumor\\tPL:" + config["platform"]
        K = 10000000,
        ID = lambda wildcards: t[wildcards.sample],
        SM = "{sample}_tumor",
        PL = config["platform"]
    threads:
        48 # set the maximum number of available cores
    shell:
        # sentieon bwa mem -M -R '{params.R}' -t {threads} -K {params.K} -o {output.sam} {input.fasta} {input.R1} {input.R2} >> {log.bwa} 2>&1
        """
        sentieon bwa mem -M -R '@RG\\tID:{params.ID}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam} {input.fasta} {input.R1} {input.R2} >> {log.bwa} 2>&1
        sentieon util sort -r {input.fasta} -i {output.sam} -o {output.bam} -t {threads} --sam2bam >> {log.sort} 2>&1
        """

'''
rule metrics_normal:
    input:
        bam = rules.mapping_normal.output.bam,
        fasta = config['fasta']
    output:
        mqm = METRICS + '{sample}.normal_mq_metrics.txt',
        qdm = METRICS + '{sample}.normal_qd_metrics.txt',
        gcs = METRICS + '{sample}.normal_gc_summary.txt',
        gcm = METRICS + '{sample}.normal_gc_metrics.txt',
        aln = METRICS + '{sample}.normal_aln_metrics.txt',
        ism = METRICS + '{sample}.normal_is_metrics.txt',
        gcp = PLOTS + '{sample}.normal_gc-report.pdf',
        qdp = PLOTS + '{sample}.normal_qd-report.pdf',
        mqp = PLOTS + '{sample}.normal_mq-report.pdf',
        isp = PLOTS + '{sample}.normal_is-report.pdf'
    log:
        LOGS + '{sample}.normal_metrics.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.fasta} -t {threads} -i {input.bam} --algo MeanQualityByCycle {output.mqm} --algo QualDistribution {output.qdm} --algo GCBias --summary {output.gcs} {output.gcm} --algo AlignmentStat --adapter_seq '' {output.aln} --algo InsertSizeMetricAlgo {output.ism} >> {log} 2>&1
        sentieon plot GCBias -o {output.gcp} {output.gcm}
        sentieon plot QualDistribution -o {output.qdp} {output.qdm}
        sentieon plot MeanQualityByCycle -o {output.mqp} {output.mqm}
        sentieon plot InsertSizeMetricAlgo -o {output.isp} {output.ism}
        """
'''

rule metrics_tumor:
    input:
        bam = rules.mapping_tumor.output.bam,
        fasta = config['fasta']
    output:
        mqm = METRICS + '{sample}.tumor_mq_metrics.txt',
        qdm = METRICS + '{sample}.tumor_qd_metrics.txt',
        gcs = METRICS + '{sample}.tumor_gc_summary.txt',
        gcm = METRICS + '{sample}.tumor_gc_metrics.txt',
        aln = METRICS + '{sample}.tumor_aln_metrics.txt',
        ism = METRICS + '{sample}.tumor_is_metrics.txt',
        gcp = PLOTS + '{sample}.tumor_gc-report.pdf',
        qdp = PLOTS + '{sample}.tumor_qd-report.pdf',
        mqp = PLOTS + '{sample}.tumor_mq-report.pdf',
        isp = PLOTS + '{sample}.tumor_is-report.pdf'
    log:
        LOGS + '{sample}.tumor_metrics.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.fasta} -t {threads} -i {input.bam} --algo MeanQualityByCycle {output.mqm} --algo QualDistribution {output.qdm} --algo GCBias --summary {output.gcs} {output.gcm} --algo AlignmentStat --adapter_seq '' {output.aln} --algo InsertSizeMetricAlgo {output.ism} >> {log} 2>&1
        sentieon plot GCBias -o {output.gcp} {output.gcm}
        sentieon plot QualDistribution -o {output.qdp} {output.qdm}
        sentieon plot MeanQualityByCycle -o {output.mqp} {output.mqm}
        sentieon plot InsertSizeMetricAlgo -o {output.isp} {output.ism}
        """

'''
rule markdup_normal:
    input:
        bam = rules.mapping_normal.output.bam,
        fasta = config['fasta']
    output:
        ns = MARKDUP + '{sample}.normal_score.txt',
        dm = MARKDUP + '{sample}.normal_dedup_metrics.txt',
        bam = BAMS + '{sample}.normal_deduped.bam',
        cm = MARKDUP + '{sample}.normal_coverage_metrics'
    log:
        LOGS + '{sample}.normal_dedup.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.ns} >> {log} 2>&1
        sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {output.ns} --metrics {output.dm} {output.bam} >> {log} 2>&1
        sentieon driver -r {input.fasta} -t {threads} -i {output.bam} --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """
'''

rule markdup_tumor:
    input:
        bam = rules.mapping_tumor.output.bam,
        fasta = config['fasta']
    output:
        ns = MARKDUP + '{sample}.tumor_score.txt',
        dm = MARKDUP + '{sample}.tumor_dedup_metrics.txt',
        bam = BAMS + '{sample}.tumor_deduped.bam',
        cm = MARKDUP + '{sample}.tumor_coverage_metrics'
    log:
        LOGS + '{sample}.tumor_dedup.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.ns} >> {log} 2>&1
        sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {output.ns} --metrics {output.dm} {output.bam} >> {log} 2>&1
        sentieon driver -r {input.fasta} -t {threads} -i {output.bam} --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """

'''
rule baserecal_normal:
    input:
        bam = rules.markdup_normal.output.bam,
        fasta = config['fasta']
    output:
        rdt = RECAL + '{sample}.normal_recal_data.table',
        post = RECAL + '{sample}.normal_recal_data.table.post',
        recal = RECAL + '{sample}.normal_recal.csv',
        rp = PLOTS + '{sample}.normal_recal_plots.pdf',

    log:
        LOGS + '{sample}.normal_recal.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.fasta} -t {threads} -i {input.bam} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.rdt} >> {log} 2>&1
        sentieon driver -r {input.fasta} -t {threads} -i {input.bam} -q {output.rdt} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.post} >> {log} 2>&1
        sentieon driver -t {threads} --algo QualCal --plot --before {output.rdt} --after {output.post} {output.recal} >> {log} 2>&1
        sentieon plot QualCal -o {output.rp} {output.recal}
        """
'''

rule baserecal_tumor:
    input:
        bam = rules.markdup_tumor.output.bam,
        fasta = config['fasta']
    output:
        rdt = RECAL + '{sample}.tumor_recal_data.table',
        post = RECAL + '{sample}.tumor_recal_data.table.post',
        recal = RECAL + '{sample}.tumor_recal.csv',
        rp = PLOTS + '{sample}.tumor_recal_plots.pdf',

    log:
        LOGS + '{sample}.tumor_recal.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.fasta} -t {threads} -i {input.bam} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.rdt} >> {log} 2>&1
        sentieon driver -r {input.fasta} -t {threads} -i {input.bam} -q {output.rdt} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.post} >> {log} 2>&1
        sentieon driver -t {threads} --algo QualCal --plot --before {output.rdt} --after {output.post} {output.recal} >> {log} 2>&1
        sentieon plot QualCal -o {output.rp} {output.recal}
        """


rule variant_calling:
    input:
        tumor_bam = rules.markdup_tumor.output.bam,
        #normal_bam = rules.markdup_normal.output.bam,
        tumor_rdt = rules.baserecal_tumor.output.rdt,
        #normal_rdt = rules.baserecal_normal.output.rdt
    output:
        vcf = RESULTS + '{sample}.tnscope.vcf.gz',
    params:
        tumor_sample = '{sample}_tumor',
        #normal_sample = '{sample}_normal'
    log:
        LOGS + '{sample}.tnscope.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {fasta} -t {threads} -i {input.tumor_bam} -q {input.tumor_rdt} --algo TNscope --tumor_sample {params.tumor_sample} {output.vcf} >> {log} 2>&1
        """
        #sentieon driver -r {fasta} -t {threads} -i {input.tumor_bam} -i {input.normal_bam} -q {input.tumor_rdt} -q {input.normal_rdt} --algo TNscope --tumor_sample {params.tumor_sample} --normal_sample {params.normal_sample} --dbsnp {dbsnp} {output.vcf} >> {log} 2>&1


#############################
# Not properly working code #
#############################

# Prints SAM in standard output and @RG is not properly written in sorted BAM
#'( sentieon bwa mem -M -R "{params.R}" '
#'-t {threads} -K {params.K} -o {output.sam} {input.fasta} {input.R1} {input.R2} || echo -n "error" ) | '
#'sentieon util sort -o {output.bam} -t {threads} --sam2bam -i - >> {log} 2>&1'

# Prints SAM into log file and @RG is still not properly written in sorted BAM
#( sentieon bwa mem -M -R '{params.R}' -t {threads} -K {params.K} -o {output.sam} {input.fasta} {input.R1} {input.R2} || echo -n "error" ) >> {log} 2>&1 | sentieon util sort -o {output.bam} -t {threads} --sam2bam -i - >> {log} 2>&1

# Modified by Massimiliano Volpe on 12/05/2021 to test the tumor-only pipeline.
