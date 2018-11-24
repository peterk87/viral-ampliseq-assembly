rule download_reference_genomes:
    output:
        fasta='references/' + config['organism'] + '.fasta',
        genbank='references/' + config['organism'] + '.genbank'
    params:
        organism=config['organism']
    log:
        'logs/scripts/download_reference_genomes/' + config['organism'] + '.log'
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/download_reference_genomes.py'
