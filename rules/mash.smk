"""
Mash screen references genomes for user-specified organism against 
"""


rule mash_sketch_reference_genomes:
    input:
        'references/{sample}.fasta'
    output:
        'references/{sample}.msh'
    params:
        kmer_size=config['mash'].get('kmer_size', 16),
        n_sketches=config['mash'].get('n_sketches', 10000)
    log:
        'logs/mash/{sample}-mash_sketch_reference_genomes.log'
    benchmark:
        'benchmarks/mash_sketch/{sample}.tsv'
    conda:
        '../envs/mash.yaml'
    shell:
        'mash sketch -i -k {params.kmer_size} -s {params.n_sketches} -o {output} {input} > {log} 2>&1'


rule mash_screen_reads_vs_references:
    input:
        reads='preprocess/fastqs/{sample}.fastq',
        sketches='references/{sample}.msh'
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
        fasta='references/{sample}.fasta'
    output:
        'mapping/{sample}/reference.fasta'
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/top_reference_by_mash.py'
