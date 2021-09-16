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

import subprocess
from pathlib import Path

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

def get_sample_name(f_path):
    '''Extract samplename from string: /home/rada/Documents/TNscope/PVAL_65_S1/exon_cov/PVAL_65_S1.exon_cov.tsv'''
    return str(Path(f_path).stem).split(".")[0]

def create_dir_if_not_exist(dir_name):
    '''Create a directory based on input string'''
    if not Path(dir_name).is_dir():
        Path(dir_name).mkdir(parents=True, exist_ok=True)


# Globals ---------------------------------------------------------------------

configfile:
    "config.yml"

# In case cpus is given by the user convert it to int
cpus = int(config['cpus'])

# Getting Sentieon license server running
# start licence server
# Update with the location of the Sentieon software package and license file
shell.prefix(f"export SENTIEON_INSTALL={config['params']['sentieon_install']}; export SENTIEON_LICENSE={config['params']['sentieon_license']}; export SENTIEON_TMPDIR={config['params']['tmp_dir']};")



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
ILS = RESULTS + 'intervalLists/'

ref_genome = config['ref_genome']
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
        expand(RESULTS + 'variant_calling/{sample}.tnscope.vcf.gz', sample=SAMPLES),
        expand(RESULTS + 'exon_cov/{sample}.exon_cov.tsv', sample=SAMPLES),
        config['workdir'] + '/exon_cov_analysis/bar_plot_all_samples.html',
        expand(LOGS + '{sample}.fastqc.log', sample=SAMPLES),
        expand(RESULTS + 'samtools_stats/{sample}.stats.tsv', sample=SAMPLES),
        expand(RESULTS + 'samtools_stats/{sample}.idxstats.tsv', sample=SAMPLES),
        expand(RESULTS + 'metrics/{sample}.dup.metrics.txt', sample=SAMPLES),
        expand(RESULTS + 'metrics/{sample}.aln.summary.metrics.txt', sample=SAMPLES),
        expand(LOGS + 'collectHsMetrics.log', sample=SAMPLES),
        expand(LOGS + '{sample}.insertSize.metrics.log', sample=SAMPLES),
        expand(LOGS + '{sample}.gcbias.metrics.log', sample=SAMPLES)

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

rule mapping_tumor:
    input:
        R1 = config['tumor'] + "/{sample}" + R1SUFFIX,
        R2 = config['tumor'] + "/{sample}" + R2SUFFIX,
        ref_genome = config['ref_genome']
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
        -o {output.sam} {input.ref_genome} {input.R1} {input.R2} >> {log.bwa} 2>&1
        $SENTIEON_INSTALL/bin/sentieon util sort -r {input.ref_genome} \
        -i {output.sam} \
        -o {output.bam} \
        -t {threads} \
        --sam2bam >> {log.sort} 2>&1
        """


rule metrics_tumor:
    input:
        bam = rules.mapping_tumor.output.bam,
        ref_genome = config['ref_genome']
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
        -r {input.ref_genome} \
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
        ref_genome = config['ref_genome']
    output:
        ns = MARKDUP + '{sample}.tumor_score.txt',
        dm = MARKDUP + '{sample}.tumor_dedup_metrics.txt',
        bam = BAMS + '{sample}.tumor_deduped.bam',
        cm = METRICS + 'cov/{sample}'
    log:
        LOGS + '{sample}.tumor_dedup.log'
    params:
        targets = config["targets"]
        #off_targets = config["off_targets"]
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
        -r {input.ref_genome} \
        -t {threads} \
        --interval {params.targets} \
        -i {output.bam} \
        --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """

rule CreateSequenceDictionary:
    input:
        config["ref_genome"]
    output:
        "references/ILS/reference.dict"
    log:
        'logs/CreateSequenceDictionary.log'
    params:
    shell:
        """
        picard CreateSequenceDictionary \ 
        --REFERENCE {input} \ 
        --OUTPUT {output} &> {log}
        """

rule bed2IntervalList:
    input:
        baits = config["baits"],
        targets = config["targets"],
        ref_dict = rules.CreateSequenceDictionary.output
    output:
        baits_IL = '/home/rada/Documents/TNscope/references/ILS/bait.interval_list',
        target_IL = '/home/rada/Documents/TNscope/references/ILS/target.interval_list'
    log:
        'logs/interval_list.log'
    shell:
        """
        gatk BedToIntervalList \
        I={input.baits} \
        O={output.baits_IL} \
        SD={input.ref_dict} &> {log}
        
        gatk BedToIntervalList \
        I={input.targets} \
        O={output.target_IL} \
        SD={input.ref_dict} &>> {log}
        """

rule samtools_stats:
    input:
        rules.markdup_tumor.output.bam
    output:
        RESULTS + 'samtools_stats/{sample}.stats.tsv'
    log:
        LOGS + '{sample}.samtools_stats.log'
    shell:
        """
        samtools stats {input} > {output}
        """

rule samtools_idxstats:
    input:
        rules.markdup_tumor.output.bam
    output:
        RESULTS + 'samtools_stats/{sample}.idxstats.tsv'
    log:
        LOGS + '{sample}.samtools_idxstats.log'
    shell:
        """
        samtools stats {input} > {output}
        """

rule collectHsMetrics:
    input:
        bam = rules.markdup_tumor.output.bam,
        baits_IL = rules.bed2IntervalList.output.baits_IL,
        target_IL = rules.bed2IntervalList.output.target_IL
    output:
        per_target_cov = METRICS + '{sample}_per_target_cov.tsv',
        per_base_cov = METRICS + '{sample}_per_base_cov.tsv',
        hs_metrics = METRICS + '{sample}_hs_metrics.txt'
    log:
        LOGS + 'collectHsMetrics.log'
    params:
        ref = config["ref_genome"],
        coverage_cap = 2000
    shell:
        """
        gatk CollectHsMetrics \
        I={input.bam} \
        O={output.hs_metrics} \
        PER_TARGET_COVERAGE={output.per_target_cov} \
        PER_BASE_COVERAGE={output.per_base_cov} \
        R={params.ref} \
        COVERAGE_CAP={params.coverage_cap} \
        BAIT_INTERVALS={input.baits_IL} \
        TARGET_INTERVALS={input.target_IL} &> {log}
        """   

rule collectDuplicateMetrics:
    input:
        rules.markdup_tumor.output.bam
    output:
        RESULTS + 'metrics/{sample}.dup.metrics.txt'
    log:
        LOGS + '{sample}.dup.metrics.log'
    shell:
        """
        picard CollectDuplicateMetrics \
        --INPUT {input} \
        --METRICS_FILE {output} \
        --REFERENCE_SEQUENCE {ref_genome} >> {log} 2>&1
        """

rule CollectAlignmentSummaryMetrics:
    input:
        rules.markdup_tumor.output.bam
    output:
        RESULTS + 'metrics/{sample}.aln.summary.metrics.txt'
    log:
        LOGS + '{sample}.aln_summary.metrics.log'
    shell:
        """
        picard CollectAlignmentSummaryMetrics \
        --INPUT {input} \
        --REFERENCE_SEQUENCE {ref_genome} \
        --OUTPUT {output} >> {log} 2>&1
        """

rule CollectGcBiasMetrics:
    input:
        rules.markdup_tumor.output.bam
    output:
        txt = RESULTS + 'metrics/{sample}.gcbias.metrics.txt',
        chart = RESULTS + 'metrics/{sample}.gcbias.metrics.pdf',
        summary = RESULTS + 'metrics/{sample}.gcbias.summary_metrics.txt'
    log:
        LOGS + '{sample}.gcbias.metrics.log'
    shell:
        """
        picard CollectGcBiasMetrics \
        --INPUT {input} \
        --REFERENCE_SEQUENCE {ref_genome} \
        --OUTPUT {output.txt} \
        --CHART_OUTPUT {output.chart} \
        --SUMMARY_OUTPUT {output.summary} \
        --ALSO_IGNORE_DUPLICATES true >> {log} 2>&1
        """

rule CollectInsertSizeMetrics:
    input:
        rules.markdup_tumor.output.bam
    output:
        ins = RESULTS + 'metrics/{sample}.insert_size_metrics.txt',
        hist = RESULTS + 'metrics/{sample}.insert_size_histogram.pdf'
    log:
        LOGS + '{sample}.insertSize.metrics.log'
    shell:
        """
        picard CollectInsertSizeMetrics \
        --INPUT {input} \
        --REFERENCE_SEQUENCE {ref_genome} \
        --OUTPUT {output.ins} \
        --Histogram_FILE {output.hist} >> {log} 2>&1
        """

rule baserecal_tumor:
    input:
        bam = rules.markdup_tumor.output.bam,
        ref_genome = config['ref_genome']
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
        $SENTIEON_INSTALL/bin/sentieon driver -r {input.ref_genome} \
        -t {threads} \
        -i {input.bam} \
        --algo QualCal \
        -k {dbsnp} \
        -k {known_Mills_indels} \
        -k {known_1000G_indels} {output.rdt} > {log} 2>&1
        $SENTIEON_INSTALL/bin/sentieon driver -r {input.ref_genome} \
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
        vcf = RESULTS + 'variant_calling/{sample}.tnscope.vcf.gz',
    params:
        tumor_sample = '{sample}_tumor',
        #normal_sample = '{sample}_normal'
    log:
        LOGS + '{sample}.tnscope.log'
    threads:
        cpus # set the maximum number of available cores
    shell:
        """
        $SENTIEON_INSTALL/bin/sentieon driver \
        -i {input.tumor_bam} \
        -r {ref_genome} \
        -t {threads} \
        -q {input.tumor_rdt} \
        --algo TNscope \
        --min_base_qual 20 \
        --tumor_sample {params.tumor_sample} {output.vcf} > {log} 2>&1
        """

rule exon_coverages:
    input:
        eoi = config['exons_of_interest'],
        tumor_bam = rules.markdup_tumor.output.bam
    output:
        RESULTS + 'exon_cov/{sample}.exon_cov.tsv'
    log:
        LOGS + '{sample}.exon_cov.log'
    params:
        "exon_cov_analysis/all_samples.tsv"
    shell:
        """
        SAMPLE=$(basename -s .tumor_deduped.bam {input.tumor_bam})
        bedtools coverage -hist \
        -a {input.eoi} \
        -b {input.tumor_bam} > {output}
        sed 's/^/{wildcards.sample}\t/' {output} >> {params}
        """

rule summarise_exon_coverages:
    input:
        samples = expand(RESULTS + 'exon_cov/{sample}.exon_cov.tsv', sample=SAMPLES),
        eoi = config['exon_names']
    output:
        bar_plot = config['workdir'] + '/exon_cov_analysis/bar_plot_all_samples.html',
    params:    
        metrics_tables = config['workdir'] + '/exon_cov_analysis/samples_with_transcript_inputs/'
    log:
        config['workdir'] + '/logs/summarise_exon_covs.log'
    run:
        import pandas as pd
        import plotly.express as px
        import plotly.io as io
        import os
 
        gene_names = pd.read_csv(input.eoi, 
                                sep="\t", 
                                names=["name", "symbol"],
                                header=0)

        header = ["chrom", "start", "end", "name", "score", "strand","depth","num_bases_at_depth","size_of_feature","pros_of_feature_at_depth"]

        fig_list = []
        
        for f_path in input.samples:
            sample = pd.read_csv(f_path, 
                                 sep="\t", 
                                 names=header)
            
            sample_name = get_sample_name(f_path)

            # Remove summary 'all' sections from the df since they have only 6 fields and
            # aren't really useful now otherwise either
            df = sample[sample.chrom != "all"].copy()
        
            # Parse and expand the transcript data column
            df[['ID', 'rest']] = df['name'].str.split('_cds_', -1, expand=True) # https://stackoverflow.com/a/39358924
            df[["exon_number", "unknown", "exon_chrom", "exon_start_pos", "exon_strand"]] = df['rest'].str.split('_', -1, expand=True)
            df[['base_ID', 'version']] = df['ID'].str.split('.', 2, expand=True)
            #sample_name = df.iloc[0][header[0]]
            df = df.drop(["name","rest"], axis = 1)
            # Assign correct datatypes to each column
            df = df.astype({
                'chrom':'str',
                'start':'int',
                'end':'int',
                'score':'float',
                'strand':'str',
                'depth':'int',
                'num_bases_at_depth':'int',
                'size_of_feature':'int',
                'pros_of_feature_at_depth':'float',
                'ID':'category',
                "exon_number":'category', 
                "unknown":'category', 
                "exon_chrom":"category", 
                "exon_start_pos":"int", 
                "exon_strand":"category"
                }
                )
            # Add gene names to the current df
            df = df.merge(right=gene_names,
                          how='left',
                          left_on='base_ID',
                          right_on='name')
            
            # Join symbols and IDs into one column so they can be visualised as one entity
            df['symbol_ID'] = df['symbol'] + " " + df['ID'].astype(str)
            # Sort the symbol_IDs so they aren't randomly presented in the figure
            df = df.sort_values(by=['symbol_ID'])

            # Duplicate rows with several bases, this enables calculating averages in a more easier way
            df = df.loc[df.index.repeat(df.num_bases_at_depth)] # https://stackoverflow.com/a/57009491
            
            # Create dir for the current sample, into which all the transcript coverage tsv:s end up
            create_dir_if_not_exist(params.metrics_tables + sample_name)

            # Get a list of symbol_IDs so they can be looped through later on
            transcripts = df['symbol_ID'].unique()
            # Split the df into "symbol_ID" subgroups
            grouped_symbol_IDs = df.groupby(df.symbol_ID)
            # This will hold a list of df:s with depth metrics data such as 
            transcripts_list = []
            
            for j in transcripts:
                c = grouped_symbol_IDs.get_group(j)
                # Scrape NCBI transcript ID
                NCBI_id = c.iloc[0]['symbol_ID']
                c.index.name = "row_no"
                c.to_csv(params.metrics_tables + sample_name + "/" + NCBI_id.replace(" ", "_") + "_cov.tsv")
                # Extract metrics values for the depth column
                c = c.describe()['depth'].to_frame(NCBI_id).T
                transcripts_list.append(c)

            # Join all metrics data into one df
            all_transcripts_metrics = (pd.concat(transcripts_list, axis=0)
                                       .rename(columns={'count': 'total_length_of_exons'})
                                       .astype({
                                             'total_length_of_exons':'int',
                                             'max':'int',
                                             'min':'int'}
                                              ))
            
            all_transcripts_metrics.index.name = "ID"
            # Create bar plots
            bar_fig = px.bar(all_transcripts_metrics.reset_index(), 
                             y='mean', 
                             x='ID', 
                             hover_data=["total_length_of_exons", "mean", 'std', 'min', '25%', '50%', '75%', 'max',"ID"],
                             title="Sample name: " + sample_name)
            
            fig_list.append(bar_fig)
            del df

        # Join all metrics tables into one
        #all_metrics = pd.concat(metrics_df_list, axis=0, ignore_index=True)
        # Give the index column own name so it looks nice in the csv
        #all_metrics.index.name = "ID"

        # https://stackoverflow.com/a/59869358
        # Write all bar plots into one html page
        with open(output.bar_plot, 'a') as f:
            for fig in fig_list:
                f.write(fig.to_html(full_html=False, 
                                    include_plotlyjs='cdn', 
                                    config= {'displaylogo': False}))
