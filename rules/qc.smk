

rule fastqc:
    input:
        'preprocess/fastqs/{sample}.fastq'
    output:
        html='qc/fastqc/{sample}.html',
        zip='qc/fastqc/{sample}_fastqc.zip'
    params: ''
    log:
        'logs/fastqc/{sample}.log'
    wrapper:
        '0.27.1/bio/fastqc'


rule quast_assembly_assessment:
    input:
        expand('assembly/spades/{sample}/contigs.fasta', sample=samples.index)
    output:
        directory('qc/quast/'),
    threads: 16
    log:
        'logs/quast.log'
    conda:
        '../envs/quast.yaml'
    shell:
        'quast.py -L -o {output} --threads {threads} {input} &> {log}'


rule multiqc:
    input:
        'qc/quast/',
        expand('qc/fastqc/{sample}_fastqc.zip', sample=samples.index),
        expand('mapping/{sample}/{sample}.bam', sample=samples.index),
        expand('preprocess/samtools/depth/{sample}.tsv', sample=samples.index),
        expand('preprocess/samtools/flagstat/{sample}.flagstat', sample=samples.index),
        expand('preprocess/samtools/idxstats/{sample}.tsv', sample=samples.index),
        expand('variant_calling/{sample}-filtered.vcf', sample=samples.index),
        expand('variant_calling/snpeff/{sample}.csv', sample=samples.index)
    output:
        'qc/multiqc.html'
    params:
        ''  # Optional: extra parameters for multiqc.
    log:
        'logs/multiqc.log'
    wrapper:
        '0.27.1/bio/multiqc'
