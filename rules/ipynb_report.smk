rule create_ipynb_report:
    input:
        depth='mapping/{sample}/{sample}-depth.tsv',
        assembly='assembly/spades-{sample}.fasta'
    output:
        ipynb='notebooks/{sample}-{tool}.ipynb',
        html='notebooks/{sample}-{tool}.html'
    script:
        '../scripts/ipynb_report.py'
