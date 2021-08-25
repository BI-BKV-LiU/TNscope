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
ILS = RESULTS + 'intervalLists/'

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
        expand(config['workdir'] + '/{sample}/multiqc_report.html', sample=SAMPLES),
        # expand(LOGS + '{sample}.DeOCov.log', sample=SAMPLES),
        expand(METRICS + '{sample}_hs_metrics.txt', sample=SAMPLES)
       

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
        config['workdir'] + '/{sample}/' + 'logs/{sample}.multiqc.log'
    params:
        config['workdir'] + '/{sample}/'
    shell:
        """
        multiqc \
        -f \
        --outdir {params} \
        {params} >> {log} 2>&1
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
        -r {input.fasta} \
        -t {threads} \
        --interval {params.targets} \
        -i {output.bam} \
        --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """

rule bed2IntervalList:
    input:
        baits = config["baits"],
        # targets = config["UCSC"]
        targets = config["targets"]
    output:
        baits_IL = '/home/rada/Documents/TNscope/references/ILS/bait.interval_list',
        # target_IL = '/home/rada/Documents/TNscope/references/ILS/UCSC.interval_list'
        target_IL = '/home/rada/Documents/TNscope/references/ILS/target.interval_list'
    log:
        'logs/interval_list.log'
    params:
        ref_dict = config["fasta"]
    threads:
        cpus
    shell:
        """
        gatk BedToIntervalList \
        I={input.baits} \
        O={output.baits_IL} \
        SD={params.ref_dict}
        
        gatk BedToIntervalList \
        I={input.targets} \
        O={output.target_IL} \
        SD={params.ref_dict}
        """

rule collectHsMetrics:
    input:
        bam = rules.markdup_tumor.output.bam,
        baits_IL = rules.bed2IntervalList.output.baits_IL,
        target_IL = rules.bed2IntervalList.output.target_IL
    output:
        hs_metrics = METRICS + '{sample}_hs_metrics.txt'
    log:
        LOGS + 'collectHsMetrics.log'
    params:
        ref = config["fasta"]
    threads:
        cpus
    shell:
        """
        gatk CollectHsMetrics \
        I={input.bam} \
        O={output.hs_metrics} \
        R={params.ref} \
        BAIT_INTERVALS={input.baits_IL} \
        TARGET_INTERVALS={input.target_IL}
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

rule exon_coverages:
    input:
        eoi = config['exons_of_interest'],
        tumor_bam = rules.markdup_tumor.output.bam
    output:
        exon_covs = RESULTS + 'exon_cov/{sample}.exon_cov.tsv'
    log:
        LOGS + '{sample}.exon_cov.log'
    shell:
        """
        SAMPLE=$(basename -s .tumor_deduped.bam {input.tumor_bam})
        bedtools coverage -hist \
        -a {input.eoi} \
        -b {input.tumor_bam} > {output}/"$SAMPLE".exon_cov.tsv
        """

rule summarise_exon_coverages:
    input:
        # exon_cov = RESULTS + 'exon_cov/{sample}.exon_cov.tsv',
        exon_cov = rules.exon_coverages.output.exon_covs,
        eoi = config['exon_names']
    params:
        bar_plot = config['workdir'] + '/bar_plot_all_samples.html',
        metrics_tables = config['workdir'] + '/all_cov_metrics.csv'
    log:       
        LOGS + '{sample}.summarise_exon_covs.rdy.log'
    run:
        import plotly
        import pandas as pd
        from pathlib import Path
        import plotly.express as px
        import plotly.io as io
        import os

        #sample_covs = Path('cov/')
        DOUBLE_MATCH = glob_wildcards("/home/rada/Documents/TNscope/{sample1}/exon_cov/{sample2}.exon_cov.tsv")
        SAMPLES = expand("/home/rada/Documents/TNscope/{sample1}/exon_cov/{sample2}.exon_cov.tsv", zip, **DOUBLE_MATCH._asdict())
        
        gene_names = pd.read_csv(input.eoi, 
                                sep="\t", 
                                names=["name", "symbol"],
                                header=0
                                )

        fig_list = []
        metrics_df_list = []
        out_path = "sample_coverages/"
        #print(type(input.exon_cov))

        #for path in sorted(sample_covs.glob("*.hist")):
        for path in SAMPLES:
            path = Path(path)
            sample_name = str(path.stem).split(".")[0]
            df = pd.read_csv(path, 
                            sep="\t", 
                            names=["chrom", "start", "end", "name", "score", "strand","depth","num_bases_at_depth","size_of_feature","pros_of_feature_at_depth"])
            df = df[df.chrom != "all"].copy()
            df[['ID', 'rest']] = df['name'].str.split('_cds_', -1, expand=True) # https://stackoverflow.com/a/39358924
            df[["exon_number", "unknown", "exon_chrom", "exon_start_pos", "exon_strand"]] = df['rest'].str.split('_', -1, expand=True)
            df[['base_ID', 'version']] = df['ID'].str.split('.', 2, expand=True)
            df = df.drop(["name","rest"], axis = 1)
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
            df = df.merge(right=gene_names,
                        how='left',
                        left_on='base_ID',
                        right_on='name'
                        )
            df['symbol_ID'] = df['symbol'] + " " + df['ID'].astype(str)
            df = df.sort_values(by=['symbol_ID'])
            df = df.loc[df.index.repeat(df.num_bases_at_depth)] # https://stackoverflow.com/a/57009491
            
            # Extract a list of all unique NCBI ID:s
            transcripts = df['symbol_ID'].unique()
            
            grouped = df.groupby(df.symbol_ID)

            transcripts_list = []
            
            for i in transcripts:
                c = grouped.get_group(i)
                # Scrape NCBI transcript ID
                NCBI_id = c.iloc[0]['symbol_ID']
                # Duplicate rows with several bases 
                c = c.loc[c.index.repeat(c.num_bases_at_depth)] # https://stackoverflow.com/a/57009491
                # Extract key values for the depth column
                c = c.describe()['depth'].to_frame(NCBI_id).T
                transcripts_list.append(c)

            # Join all key data into one df
            all_transcripts = (pd.concat(transcripts_list, axis=0)
                            .rename(columns={'count': 'total_length_of_exons'})
                            .astype({
                                    'total_length_of_exons':'int',
                                    'max':'int',
                                    'min':'int'}
                                    )                  
                            )
            
            all_transcripts.index.name = "ID"
            # Create bar plots
            bar_fig = px.bar(all_transcripts.reset_index(), 
                            x='ID', 
                            y='mean', 
                            hover_data=["total_length_of_exons", "mean", 'std', 'min', '25%', '50%', '75%', 'max',"ID"],
                            title="Sample name: " + sample_name
                            )
            #bar_fig.write_html(out_path + "bar/" + sample_name + "_bar" + ".html", 
            #                   config= {'displaylogo': False})
            
            fig_list.append(bar_fig)
            
            # Create tables with metrics
            all_transcripts['sample_name'] = sample_name
            metrics_df_list.append(all_transcripts)


        all_metrics = pd.concat(metrics_df_list, axis=0, ignore_index=True)
        all_metrics.index.name = "ID"
        out_path = Path(out_path + "metrics/all_metrics.csv")
        all_metrics.to_csv(params.metrics_tables)


        print("Printed")

        all_bar_plots = params.bar_plot
        # # Remove the previous version of html file, very inelegant but works...
        if os.path.exists(all_bar_plots):
            os.remove(all_bar_plots)

        # https://stackoverflow.com/a/59869358
        # Write all bar plots into one html page
        with open(all_bar_plots, 'a') as f:
            for fig in fig_list:
                f.write(fig.to_html(full_html=False, 
                                    include_plotlyjs='cdn', 
                                    config= {'displaylogo': False}))
