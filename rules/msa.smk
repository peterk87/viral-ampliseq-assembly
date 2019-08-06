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
    params:
        '--auto'
    log: 
        'logs/mafft.log'
    threads: -1
    shell:
        'mafft --thread {threads} {params} {input} > {output} 2> >(tee -a {log} >&2)'

