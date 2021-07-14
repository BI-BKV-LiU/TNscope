################################################################################
## TNscope Pipeline
## Somatic variant calling with Sentieon TNscope (from fastq files to vcf)
##
## Authors: Massimiliano Volpe, Jyotirmoy Das, Lauri Mesilaakso (Adjustments in "ro-implementation"-branch)
## Email: massimiliano.volpe@liu.se, jyotirmoy.das@liu.se, lauri.mesilaakso@regionostergotland.se
## Date: 14/07/2021
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

configfile:
    "config.yml"

# In case cpus is given by the user convert it to int
cpus = int(config['cpus'])

# Getting Sentieon license server running
# start licence server
# Update with the location of the Sentieon software package and license file
shell.prefix(f"export SENTIEON_INSTALL={config['params']['sentieon_install']}; export SENTIEON_LICENSE={config['params']['sentieon_license']}; export SENTIEON_TMPDIR={config['params']['tmp_dir']};")

import subprocess

subprocess.run(config['params']['sentieon_install'] + "/bin/sentieon licsrvr --start --log logs/logs.txt " + config['params']['sentieon_license'], shell=True)


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
FASTQC = RESULTS + 'fastqc/'

fasta = config['fasta']
dbsnp = config['dbsnp']
known_Mills_indels = config['known_Mills_indels']
known_1000G_indels = config['known_1000G_indels']

t = {}

for sample in SAMPLES:
    id_maker(config['tumor'], '/', sample, R1SUFFIX, t)


# Rules -----------------------------------------------------------------------
rule all:
    input:
        expand(BAMS + "{sample}.tumor_sorted.bam", sample=SAMPLES),
        expand(LOGS + '{sample}.tumor_metrics.log', sample=SAMPLES),
        expand(BAMS + '{sample}.tumor_deduped.bam', sample=SAMPLES),
        expand(RECAL + '{sample}.tumor_recal_data.table', sample=SAMPLES),
        expand(RESULTS + '{sample}.tnscope.vcf.gz', sample=SAMPLES),
        expand(LOGS + '{sample}.fastqc.log', sample=SAMPLES),
        expand(config['workdir'] + '/{sample}/multiqc_report.html', sample=SAMPLES)
       

rule fastqc:
    input:
        R1 = config['tumor'] + "/{sample}" + R1SUFFIX,
        R2 = config['tumor'] + "/{sample}" + R2SUFFIX
    output:
        FASTQC + '{sample}_R1_001_fastqc.zip',
        FASTQC + '{sample}_R2_001_fastqc.zip',
        FASTQC + '{sample}_R1_001_fastqc.html',
        FASTQC + '{sample}_R2_001_fastqc.html'
    log:
        LOGS + '{sample}.fastqc.log'
    params:
        outdir = FASTQC
    threads:
        cpus
    shell:
        """
        fastqc {input.R1} {input.R2} \
        --format fastq \
        --threads {threads} \
        --outdir {params.outdir} >> {log} 2>&1
        """

rule multiqc:
    input:
        FASTQC + '{sample}_R1_001_fastqc.zip'
    output:
        config['workdir'] + '/{sample}/' + 'multiqc_report.html'
    log:
        config['workdir'] + '/{sample}/' + 'logs/multiqc.log'
    params:
        config['workdir'] + '/{sample}/'
    shell:
        """
        multiqc \
        -f \
        --outdir {params} \
        . >> {log} 2>&1
        """

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
        cpus # set the maximum number of available cores
    shell:
        """
        $SENTIEON_INSTALL/bin/sentieon bwa mem -M \
        -R '@RG\\tID:{params.ID}\\tSM:{params.SM}\\tPL:{params.PL}' \
        -t {threads} \
        -K {params.K} \
        -o {output.sam} {input.fasta} {input.R1} {input.R2} >> {log.bwa} 2>&1
        $SENTIEON_INSTALL/bin/sentieon util sort -r {input.fasta} \
        -i {output.sam} \
        -o {output.bam} \
        -t {threads} \
        --sam2bam >> {log.sort} 2>&1
        """


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
        cpus # set the maximum number of available cores
    shell:
        """
        $SENTIEON_INSTALL/bin/sentieon driver \
        -r {input.fasta} \
        -t {threads} \
        -i {input.bam} \
        --algo MeanQualityByCycle {output.mqm} \
        --algo QualDistribution {output.qdm} \
        --algo GCBias \
        --summary {output.gcs} {output.gcm} \
        --algo AlignmentStat \
        --adapter_seq '' {output.aln} \
        --algo InsertSizeMetricAlgo {output.ism} >> {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon plot GCBias -o {output.gcp} {output.gcm}
        $SENTIEON_INSTALL/bin/sentieon plot QualDistribution -o {output.qdp} {output.qdm}
        $SENTIEON_INSTALL/bin/sentieon plot MeanQualityByCycle -o {output.mqp} {output.mqm}
        $SENTIEON_INSTALL/bin/sentieon plot InsertSizeMetricAlgo -o {output.isp} {output.ism}
        """


rule markdup_tumor:
    input:
        bam = rules.mapping_tumor.output.bam,
        fasta = config['fasta']
    output:
        ns = MARKDUP + '{sample}.tumor_score.txt',
        dm = MARKDUP + '{sample}.tumor_dedup_metrics.txt',
        bam = BAMS + '{sample}.tumor_deduped.bam',
        cm = METRICS + 'cov/{sample}'
    log:
        LOGS + '{sample}.tumor_dedup.log'
    params:
        targets = config["targets"],
        off_targets = config["off_targets"]
    threads:
        cpus # set the maximum number of available cores
    shell:
        """
        $SENTIEON_INSTALL/bin/sentieon driver \
        -t {threads} \
        -i {input.bam} \
        --algo LocusCollector \
        --fun score_info {output.ns} >> {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon driver \
        -t {threads} \
        -i {input.bam} \
        --algo Dedup \
        --rmdup \
        --score_info {output.ns} \
        --metrics {output.dm} {output.bam} >> {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon driver \
        -r {input.fasta} \
        -t {threads} \
        --interval {params.targets} \
        --interval {params.off_targets} \
        -i {output.bam} \
        --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """


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
        cpus # set the maximum number of available cores
    shell:
        """
        $SENTIEON_INSTALL/bin/sentieon driver -r {input.fasta} \
        -t {threads} \
        -i {input.bam} \
        --algo QualCal \
        -k {dbsnp} \
        -k {known_Mills_indels} \
        -k {known_1000G_indels} {output.rdt} >> {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon driver -r {input.fasta} \
        -t {threads} \
        -i {input.bam} \
        -q {output.rdt} \
        --algo QualCal \
        -k {dbsnp} \
        -k {known_Mills_indels} \
        -k {known_1000G_indels} {output.post} >> {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon driver -t {threads} \
        --algo QualCal \
        --plot \
        --before {output.rdt} \
        --after {output.post} {output.recal} >> {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon plot QualCal -o {output.rp} {output.recal}
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
        cpus # set the maximum number of available cores
    shell:
        """
        $SENTIEON_INSTALL/bin/sentieon driver -r {fasta} \
        -t {threads} \
        -i {input.tumor_bam} \
        -q {input.tumor_rdt} \
        --algo TNscope \
        --tumor_sample {params.tumor_sample} {output.vcf} >> {log} 2>&1
        """
