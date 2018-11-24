# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

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
        expand('variant_calling/{sample}-vcf.tsv', sample=samples.index)


# include rules for each step in workflow
include: "rules/preprocess.smk"
include: "rules/mash.smk"
include: "rules/download_references.smk"
include: "rules/bwa.smk"
include: "rules/spades.smk"
include: "rules/variant_calling.smk"
include: "rules/qc.smk"
