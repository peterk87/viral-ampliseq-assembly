
if config['fast_tree']:
    rule clearcut:
        input:
            'msa/alignment.fasta'
        output:
            newick='msa/alignment.fasta.treefile'
        benchmark:
            'benchmarks/clearcut.tsv'
        conda:
            '../envs/clearcut.yaml'
        shell:
            'clearcut --alignment --DNA --in={input} --out={output}'
else:
    rule iqtree:
        input:
            'msa/alignment.fasta'
        output:
            report='msa/alignment.fasta.iqtree',
            newick='msa/alignment.fasta.treefile'
        params:
            model='GTR+F+R3'
        log:
            'logs/iqtree.log'
        benchmark:
            'benchmarks/iqtree.tsv'
        conda:
            '../envs/iqtree.yaml'
        shell:
            'iqtree -s {input} -nt 8 -bb 1000 -m {params.model} -redo | tee {log}'

rule render_tree:
    input:
        newick='msa/alignment.fasta.treefile',
        references='references/' + config['organism'] + '.genbank'
    output:
        rooted_newick=report('phylogeny/rooted_tree.newick',
                             caption='../report/results/phylogeny_newick.rst',
                             category='Phylogeny'),
        tree_png='phylogeny/rooted_tree.png',
        tree_svg='phylogeny/rooted_tree.svg',
        metadata_tsv='phylogeny/genome-metadata.tsv'
    params:
        samples=samples.index
    conda:
        '../envs/python_ete3.yaml'
    script:
        '../scripts/render_tree.py'
