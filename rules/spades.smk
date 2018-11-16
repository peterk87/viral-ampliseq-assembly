"""
Rule for assembly with SPAdes
"""


rule spades_assembly:
    input:
        forward=get_fastq_forward,
        reverse=get_fastq_reverse
    output:
        contigs='assembly/spades/{sample}/contigs.fasta',
        assembly_graph='assembly/spades/{sample}/contigs.gfa',
        scaffolds='assembly/spades/{sample}/scaffolds.fasta',
        dir='assembly/spades/{sample}'
    threads: config['spades']['threads']
    params:
        tmp=config['spades']['tmp'],
        careful='--careful' if config['spades']['careful'] else '',
        meta='--meta' if config['spades']['meta'] else '',
        args=config['spades']['args']
    log: 
        'logs/spades/{sample}.log'
    benchmark:
        'benchmark/spades/{sample}.tsv'
    conda:
        '../envs/spades.yaml'
    shell:
        """
        spades.py \
        -1 {input.forward} \
        -2 {input.reverse} \
        -o {output.dir} \
        --tmp-dir {params.tmp}\
        --threads {threads} \
        {params.careful} \
        {params.args} \
        > {log} 2>&1
        """
