
if config['fast_tree']:
    rule clearcut:
        input:
            'msa/alignment.fasta'
        output:
            newick='msa/alignment.fasta.treefile'
        benchmark:
            'benchmarks/clearcut.tsv'
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
        shell:
            'iqtree -s {input} -nt 8 -bb 1000 -m {params.model} -redo | tee {log}'

rule render_tree:
    input:
        newick='msa/alignment.fasta.treefile',
        references='references/' + config['organism'] + '.genbank'
    output:
        tree_html=report('phylogeny/tree.html', 
            caption='../report/results/phylogeny_newick.rst',
            category='Phylogenetic tree'),
        metadata_tsv='phylogeny/genome-metadata.tsv'
    log:
        'logs/shiptv.log'
    benchmark:
        'benchmarks/shiptv.tsv'
    shell:
        'shiptv -r {input.references} -n {input.newick} -o {output.tree_html} -m {output.metadata_tsv}'
