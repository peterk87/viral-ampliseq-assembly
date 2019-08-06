"""
This is the main Snakemake entry point for the viral-ampliseq-assembly workflow.

New rules must be included here and target files specified under rule all or another rule. 
"""

# --use-singularity cmdline arg set
singularity: "shub://peterk87/viral-ampliseq-assembly"
# --use-conda will only create separate envs for each tool wrapper and not 
# create an environment common to all rules
# In order to run workflow with workflow global Conda env: 
# conda env create -f environment.yml && \
#   conda activate viral-ampliseq-assembly-1.0.0 && \
#   conda install snakemake
# then run your workflow inside your activated env with no --use-conda flag
conda: "environment.yml"

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
        expand('assembly/spades-{sample}.fasta', sample=samples.index),
        expand('variant_calling/{sample}-vcf.tsv', sample=samples.index),
        expand('consensus/{sample}.fasta', sample=samples.index),
        'msa/alignment.fasta',
        'msa/alignment.fasta.treefile',
        'phylogeny/tree.html'


# include rules for each step in workflow
include: "rules/preprocess.smk"
include: "rules/mash.smk"
include: "rules/download_references.smk"
include: "rules/bwa.smk"
include: "rules/spades.smk"
include: "rules/variant_calling.smk"
include: "rules/qc.smk"
include: "rules/consensus.smk"
include: "rules/msa.smk"
include: "rules/phylogeny.smk"
