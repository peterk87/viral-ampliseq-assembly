"""
Remove duplicated reads, get BAM alignment stats and output the reads in FASTQ format.
"""


rule mark_duplicates:
    input:
        get_bam_file
    output:
        bam="preprocess/dedup/{sample}.bam",
        metrics="preprocess/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.27.1/bio/picard/markduplicates"


rule samtools_index_bam_initial:
    input:
        "preprocess/dedup/{sample}.bam"
    output:
        'preprocess/samtools/index/{sample}.done'
    threads: config['samtools']['threads']
    benchmark:
        'benchmarks/samtools_index_bam/{sample}.tsv'
    shell:
        """
        samtools index -@ {threads} {input}
        touch {output}
        """


rule samtools_flagstat_bam_initial:
    input:
        "preprocess/dedup/{sample}.bam"
    output:
        'preprocess/samtools/flagstat/{sample}.flagstat'
    threads: config['samtools']['threads']
    shell:
        'samtools flagstat -@ {threads} {input} > {output}'


rule samtools_idxstats_bam_initial:
    """
    Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index. 
    The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
    """
    input:
        bam="preprocess/dedup/{sample}.bam",
        bai_done='preprocess/samtools/index/{sample}.done'
    output:
        'preprocess/samtools/idxstats/{sample}.tsv'
    shell:
        'samtools idxstats {input.bam} > {output}'


rule process_samtools_idxstats_bam_initial:
    input:
        'preprocess/samtools/idxstats/{sample}.tsv'
    output:
        sorted='preprocess/samtools/idxstats/{sample}-sorted.tsv',
        top_mapped='preprocess/samtools/idxstats/{sample}-top_mapped.txt'
    script:
        '../scripts/process_samtools_idxstats.py'


rule samtools_depth_bam_initial:
    """
    Calculate depth for all positions including those with zero coverage and no
    limit on coverage.
    """
    input:
        "preprocess/dedup/{sample}.bam"
    output:
        'preprocess/samtools/depth/{sample}.tsv'
    benchmark:
        'benchmarks/samtools_depth/{sample}.tsv'
    shell:
        'samtools depth -a -d 0 {input} > {output}'


rule process_samtools_depth_bam_initial:
    input:
        depth='preprocess/samtools/depth/{sample}.tsv',
        idxstats='preprocess/samtools/idxstats/{sample}-sorted.tsv'
    output:
        genome_extent='preprocess/samtools/depth/{sample}-genome_extent.tsv',
        extent='preprocess/samtools/depth/{sample}-extent.tsv'
    script:
        '../scripts/process_samtools_depth.py'


rule bam_to_fastq:
    input:
        "preprocess/dedup/{sample}.bam"
    output:
        'preprocess/fastqs/{sample}.fastq'
    log:
        'logs/samtools/bam2fq-{sample}.log'
    shell:
        'samtools bam2fq {input} > {output} 2> {log}'


rule trimmomatic_se:
    input:
        'preprocess/fastqs/{sample}.fastq'
    output:
        'preprocess/trimmed_fastqs/{sample}.fastq'
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"],
        # optional parameters
        extra=""
    wrapper:
        "0.27.1/bio/trimmomatic/se"
