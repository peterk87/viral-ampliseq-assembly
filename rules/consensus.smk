
rule replace_ambiguous_bases_in_reference:
    '''
    Replace ambiguous bases with N so that bcftools consensus won't have issues
    with the REF allele and reference sequence not matching up.
    '''
    input:
        'mapping/{sample}/reference.fasta'
    output:
        'consensus/{sample}/reference.fasta'
    log:
        'logs/scripts/replace_ambiguous_bases_in_reference/{sample}.log'
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/replace_ambiguous_bases_in_reference.py'


rule bcftools_convert:
    input:
        'variant_calling/snpeff/{sample}.vcf'
    output:
        'consensus/{sample}.vcf.gz'
    conda:
        '../envs/snippy.yaml'
    shell:
        'bcftools convert -Oz -o {output} {input}'

rule bcftools_index:
    input:
        'consensus/{sample}.vcf.gz'
    output:
        'consensus/{sample}.vcf.gz.csi'
    conda:
        '../envs/snippy.yaml'
    shell:
        '''
        bcftools index -f {input}
        touch {output}
        '''

rule bcftools_consensus:
    input:
        ref='consensus/{sample}/reference.fasta',
        vcfgz='consensus/{sample}.vcf.gz',
        vcfgz_index='consensus/{sample}.vcf.gz.csi'
    output:
        'consensus/{sample}-raw.fasta'
    conda:
        '../envs/snippy.yaml'
    shell:
        'bcftools consensus -f {input.ref} -o {output} {input.vcfgz}'

rule fix_consensus:
    input:
        fasta='consensus/{sample}-raw.fasta',
        depth='mapping/{sample}/{sample}-depth.tsv'
    output:
        fasta='consensus/{sample}.fasta'
    params:
        unmapped_char='-'
    conda:
        '../envs/python_biopython.yaml'
    script:
        '../scripts/fix_zero_coverage_depth_in_consensus.py'
