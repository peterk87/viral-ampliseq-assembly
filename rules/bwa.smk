"""
Read mapping to top reference rules
"""

rule replace_ambiguous_bases_in_reference:
    '''
    Replace ambiguous bases with N so that bcftools consensus won't have issues
    with the REF allele and reference sequence not matching up.
    '''
    input:
        'references/{sample}/reference.fasta'
    output:
        'references/{sample}/reference-no_ambig.fasta'
    log:
        'logs/scripts/replace_ambiguous_bases_in_reference/{sample}.log'
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/replace_ambiguous_bases_in_reference.py'


rule bwa_index:
    input:
        'references/{sample}/reference-no_ambig.fasta'
    output:
        'references/{sample}/reference-no_ambig.fasta.bwt'
    log:
        'logs/bwa_index/{sample}.log'
    benchmark:
        'benchmarks/bwa_index/{sample}.tsv'
    conda:
        '../envs/bwa.yaml'
    shell:
        '''
        bwa index {input} > {log} 2>&1
        test -r {output}
        '''


rule samtools_faidx:
    input:
        'references/{sample}/reference-no_ambig.fasta'
    output:
        'references/{sample}/reference-no_ambig.fasta.fai'
    conda:
        '../envs/bwa.yaml'
    shell:
        '''
        samtools faidx {input}
        test -r {output}
        '''


rule bwa_mem:
    """
    Read mapping with BWA-MEM, filtering SAM file for soft and hard clipped 
    alignments with samclip, samtools sorting and filtering of duplicate reads.
    """
    input:
        ref='references/{sample}/reference-no_ambig.fasta',
        ref_bwa_index='references/{sample}/reference-no_ambig.fasta.bwt',
        ref_fai='references/{sample}/reference-no_ambig.fasta.fai',
        reads='preprocess/fastqs/{sample}.fastq'
    output:
        'mapping/{sample}/{sample}.bam'
    params:
        maxsoft=10 # max soft clipping to allow
    threads: config['bwa']['threads']
    log:
        'logs/bwa_mem/{sample}.log'
    benchmark:
        'benchmarks/bwa_mem/{sample}.tsv'
    conda:
        '../envs/bwa.yaml'
    shell:
        '''
        (bwa mem -t {threads} {input.ref} {input.reads} \
        | samclip --max {params.maxsoft} --ref {input.ref_fai} \
        | samtools sort --threads {threads} \
        | samtools markdup -r -s - - \
        > {output}) \
          2> {log}
        '''


rule samtools_index:
    input:
        'mapping/{sample}/{sample}.bam'
    output:
        'mapping/{sample}/{sample}.bam.bai'
    threads: config['samtools']['threads']
    conda:
        '../envs/bwa.yaml'
    shell:
        '''
        samtools index -@ {threads} {input}
        test -r {output}
        '''


rule samtools_flagstat:
    input:
        'mapping/{sample}/{sample}.bam'
    output:
        'mapping/{sample}/{sample}.flagstat'
    threads: config['samtools']['threads']
    conda:
        '../envs/bwa.yaml'
    shell:
        'samtools flagstat -@ {threads} {input} > {output}'


rule samtools_idxstats:
    """
    Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index. 
    The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
    """
    input:
        bam='mapping/{sample}/{sample}.bam',
        bai_done='mapping/{sample}/{sample}.bam.bai'
    output:
        'mapping/{sample}/{sample}-idxstats.tsv'
    conda:
        '../envs/bwa.yaml'
    shell:
        'samtools idxstats {input.bam} > {output}'


rule process_samtools_idxstats:
    input:
        'mapping/{sample}/{sample}-idxstats.tsv'
    output:
        sorted='mapping/{sample}/{sample}-idxstats-sorted.tsv',
        top_mapped='mapping/{sample}/{sample}-idxstats-top_mapped.txt'
    conda:
        '../envs/python_pandas.yaml'
    script:
        '../scripts/process_samtools_idxstats.py'


rule samtools_depth:
    """
    Calculate depth for all positions including those with zero coverage and no
    limit on coverage.
    """
    input:
        'mapping/{sample}/{sample}.bam'
    output:
        'mapping/{sample}/{sample}-depth.tsv'
    conda:
        '../envs/bwa.yaml'
    shell:
        'samtools depth -a -d 0 {input} > {output}'


rule process_samtools_depth:
    input:
        depth='mapping/{sample}/{sample}-depth.tsv',
        idxstats='mapping/{sample}/{sample}-idxstats-sorted.tsv'
    output:
        genome_extent=report('mapping/{sample}/{sample}-genome_extent.tsv',
                             category='Read mapping to reference genome'),
        extent='mapping/{sample}/{sample}-extent.tsv'
    conda:
        '../envs/python_pandas.yaml'
    script:
        '../scripts/process_samtools_depth.py'

