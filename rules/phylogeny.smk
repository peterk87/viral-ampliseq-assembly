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
        rooted_newick='phylogeny/rooted_tree.newick',
        tree_png='phylogeny/rooted_tree.png',
        tree_svg='phylogeny/rooted_tree.svg'
    params:
        samples=samples.index,
        tree_img_width=2000,
        tree_img_height=2000
    conda:
        '../envs/python_ete3.yaml'
    script:
        '../scripts/render_tree.py'
