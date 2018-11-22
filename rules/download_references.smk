rule download_reference_genomes:
    input:
        get_bam_file
    output:
        fasta='references/{sample}.fasta',
        genbank='references/{sample}.genbank'
    params:
        organism=config['organism']
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/download_reference_genomes.py'
