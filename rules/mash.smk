"""
Mash screen references genomes for user-specified organism against 
"""


rule mash_sketch_reference_genomes:
    input:
        'references/' + config['organism'] + '.fasta'
    output:
        'references/' + config['organism'] + '.msh'
    params:
        kmer_size=config['mash'].get('kmer_size', 16),
        n_sketches=config['mash'].get('n_sketches', 10000)
    log:
        'logs/mash/' + config['organism'] + '-mash_sketch_reference_genomes.log'
    benchmark:
        'benchmarks/mash_sketch/' + config['organism'] + '.tsv'
    conda:
        '../envs/mash.yaml'
    shell:
        'mash sketch -i -k {params.kmer_size} -s {params.n_sketches} -o {output} {input} > {log} 2>&1'


rule mash_screen_reads_vs_references:
    input:
        reads='preprocess/fastqs/{sample}.fastq',
        sketches='references/' + config['organism'] + '.msh'
    output:
        'preprocess/mash/{sample}-screen_references.tsv'
    threads:
        config['mash']['threads']
    log:
        'logs/mash/{sample}-mash_screen_reads_vs_references.log'
    benchmark:
        'benchmarks/mash_screen/{sample}.tsv'
    conda:
        '../envs/mash.yaml'
    shell:
        'mash screen -p {threads} {input.sketches} {input.reads} > {output} 2> {log}'


rule sort_mash_screen:
    input:
        'preprocess/mash/{sample}-screen_references.tsv'
    output:
        'preprocess/mash/{sample}-screen_references-sorted.tsv'
    shell:
        '''
        echo -ne "identity\tmatching_sketches\tmultiplicity\tpvalue\tmatch\tmatch_comment\n" > {output}
        sort -grk1 {input} >> {output}
        '''


rule top_reference_by_mash:
    input:
        mash='preprocess/mash/{sample}-screen_references-sorted.tsv',
        genbank='references/' + config['organism'] + '.genbank'
    output:
        fasta='mapping/{sample}/reference.fasta',
        genbank='mapping/{sample}/reference.genbank',
        gff='mapping/{sample}/reference.gff'
    log:
        'logs/scripts/top_reference_by_mash/{sample}.log'
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/top_reference_by_mash.py'
