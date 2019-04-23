#!/usr/bin/env python3
from typing import List, Dict, Union

import pandas as pd
from ete3 import Tree, faces, TreeStyle, NodeStyle, TreeNode
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# INPUT
input_references_genbank = snakemake.input.references
input_newick = snakemake.input.newick
# OUTPUT
output_rooted_newick = snakemake.output.rooted_newick
output_tree_png = snakemake.output.tree_png
output_tree_svg = snakemake.output.tree_svg
output_metadata_tsv = snakemake.output.metadata_tsv
# PARAMETERS
sample_names: List[str] = [x for x in snakemake.params.samples]


def trunc(s: str, w: int = 50) -> str:
    """Truncate a string `s` to length `w` and add '...' if longer than `w`"""
    n = len(s)
    if n > w:
        return s[:(w-3)] + '...'
    else:
        return s

def genbank_source_metadata(rec: SeqRecord) -> Dict[str, str]:
    """Get source feature metadata dictionary for a SeqRecord"""
    return {k: v[0] if v is not None and len(v) == 1 else v for k,v in rec.features[0].qualifiers.items()}

def highlight_user_samples(tree: TreeNode, sample_names: List[str]) -> None:
    node_style = NodeStyle()
    node_style["bgcolor"] = "LightSteelBlue"
    for leaf_name in sample_names:
        node = tree.get_leaves_by_name(leaf_name)[0]
        node.set_style(node_style)


id_to_rec = {r.id: r for r in SeqIO.parse(input_references_genbank, 'genbank')}
df_metadata = pd.DataFrame({gid: genbank_source_metadata(rec) for gid, rec in id_to_rec.items()}).transpose()
if 'isolate' in df_metadata and 'strain' in df_metadata:
    df_metadata['strain'] = df_metadata['isolate'].combine_first(df_metadata['strain'])

metadata_columns = '''
strain
note
country
collection_date
host
isolation_source
'''.strip().split('\n')
# only use columns present in the reference genome metadata
metadata_columns = [x for x in metadata_columns if x in df_metadata]

# Read phylogenetic tree newick file using ete3
tree = Tree(newick=input_newick)
# Calculate the midpoint node
midpoint_node = tree.get_midpoint_outgroup()
# Set midpoint as output
tree.set_outgroup(midpoint_node)
# Output re-rooted tree
tree.write(outfile=output_rooted_newick)

# Save phylogeny ordered genome metadata table
df_metadata_out = df_metadata.loc[list(tree.iter_leaf_names()), metadata_columns]
df_metadata_out.to_csv(output_metadata_tsv, sep='\t')

def tree_style_layout(node):
    if node.is_leaf() and node.name in df_metadata.index:
        for i, col_name in enumerate(metadata_columns):
            genome_info = df_metadata.loc[node.name, col_name]
            if pd.isnull(genome_info): continue
            col_face = faces.TextFace(trunc(genome_info, w=50), fsize=10)
            col_face.margin_top = 1
            col_face.margin_bottom = 1
            col_face.border.margin = 1
            col_face.margin_right = 5
            # Note that this faces is added in "aligned" mode
            faces.add_face_to_node(col_face, node, column=i, aligned=True)
        # Sets the style of leaf nodes
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"


highlight_user_samples(tree, sample_names)

tree_style = TreeStyle()
tree_style.show_leaf_name = True
tree_style.show_branch_support = True
tree_style.layout_fn = tree_style_layout
tree.render(output_tree_png, tree_style=tree_style)
tree.render(output_tree_svg, tree_style=tree_style)
