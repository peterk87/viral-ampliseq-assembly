"""
Rule for assembly with SPAdes
"""


rule spades_assembly:
    input:
        'preprocess/trimmed_fastqs/{sample}.fastq' if config['trim'] else 'preprocess/fastqs/{sample}.fastq'
    output:
        contigs='assembly/spades/{sample}/contigs.fasta',
        scaffolds='assembly/spades/{sample}/scaffolds.fasta'
    threads: config['spades'].get('threads', 16)
    params:
        tmp=config['spades'].get('tmp', '/tmp/viral-ampliseq-assembly-spades'),
        careful='--careful' if config['spades'].get('careful', None) else '',
        args=config['spades'].get('args', '')
    log: 
        'logs/spades/{sample}.log'
    benchmark:
        'benchmarks/spades/{sample}.tsv'
    shell:
        '''
        OUTDIR=$(dirname {output.contigs})
        spades.py --iontorrent -s {input} -o $OUTDIR \
          --tmp-dir {params.tmp} --threads {threads} {params.careful} \
          {params.args} &> {log}
        '''

rule rename_denovo_assembly:
    input:
        'assembly/spades/{sample}/contigs.fasta'
    output:
        report('assembly/spades-{sample}.fasta',
               caption='../report/results/assembly.rst',
               category='De Novo Assembly')
    shell:
        '''
        cp {input} {output}
        '''
