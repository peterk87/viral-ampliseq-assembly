rule prepare_msa_input:
    input:
        expand('consensus/{sample}.fasta', sample=samples.index),
        'references/' + config['organism'] + '.fasta'
    output:
        'msa/samples-pre-aln.fasta'
    shell:
        'cat {input} > {output}'

rule mafft_msa:
    input:
        'msa/samples-pre-aln.fasta'
    output:
        'msa/alignment.fasta'
    log: 
        'logs/mafft.log'
    conda:
        '../envs/mafft.yaml'
    threads: 56
    shell:
        'mafft --globalpair --thread {threads} {input} > {output} 2> >(tee -a {log} >&2)'

