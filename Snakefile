# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

# include common functions to all workflows
include: "rules/common.smk"


rule all:
    input:
        expand('preprocess/samtools/flagstat/{sample}.flagstat', sample=samples.index),
        expand('preprocess/samtools/idxstats/{sample}.tsv', sample=samples.index)


# include rules for each step in workflow
include: "rules/preprocess.smk"
