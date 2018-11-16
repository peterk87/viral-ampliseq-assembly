
rule bwa_index:
    input:
        "references/{sample}/" + config['organism'] "-references.fasta"
    output:
        "mapped/{sample}/bwa_index.done"
    log:
        "logs/bwa_index/{sample}-" + config['organism'] + ".log"
    benchmark:
        "logs/bwa_index/{sample}-" + config['organism'] + ".benchmark.tsv"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input} > {log} 2>&1
        touch {output}
        """

rule bwa_mem:
    input:
        ref="references/{sample}/" + config['organism'] "-references.fasta",
        reads='preprocess/samtools/view/{sample}.fastq',
        bwa_index_done="mapped/{sample}/bwa_index.done"
    output:
        "mapped/{sample}/{sample}.bam"
    threads:
        config['bwa_mem']['threads']
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "logs/bwa_mem/{sample}.benchmark.tsv"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        (bwa mem \
          -t {threads} \
          {input.ref} \
          {input.reads} \
        | samtools sort --threads {threads} -o {output}) \
          2> {log}
        """

rule samtools_flagstat:
    input: "mapped/{tool}/{sample}.bam"
    output: "mapped/{tool}/{sample}.bam.flagstat"
    wrapper:
        "0.27.1/bio/samtools/flagstat"

rule samtools_depth:
    input: "mapped/{tool}/{sample}.bam"
    output: "mapped/{tool}/{sample}-depth.tsv"
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "logs/samtools_depth/{sample}-{tool}.benchmark.tsv"
    shell:
        """
        samtools depth -d 0 {input} > {output}
        """
