"""
This is the main Snakemake entry point for the viral-ampliseq-assembly workflow.

New rules must be included here and target files specified under rule all or another rule. 
"""

# use miniconda3 docker container if --use-singularity cmdline arg set
singularity: "docker://continuumio/miniconda3:4.4.10"

# include common functions to all workflows
include: "rules/common.smk"


rule all:
    input:
        'qc/multiqc.html',
        expand('preprocess/samtools/flagstat/{sample}.flagstat', sample=samples.index),
        expand('preprocess/samtools/idxstats/{sample}.tsv', sample=samples.index),
        expand('preprocess/samtools/depth/{sample}-extent.tsv', sample=samples.index),
        expand('preprocess/fastqs/{sample}.fastq', sample=samples.index),
        expand('preprocess/mash/{sample}-screen_references-sorted.tsv', sample=samples.index),
        expand('mapping/{sample}/{sample}-extent.tsv', sample=samples.index),
        expand('assembly/spades/{sample}/contigs.fasta', sample=samples.index),
        expand('variant_calling/{sample}-vcf.tsv', sample=samples.index),
        expand('consensus/{sample}.fasta', sample=samples.index)


# include rules for each step in workflow
include: "rules/preprocess.smk"
include: "rules/mash.smk"
include: "rules/download_references.smk"
include: "rules/bwa.smk"
include: "rules/spades.smk"
include: "rules/variant_calling.smk"
include: "rules/qc.smk"
include: "rules/consensus.smk"
