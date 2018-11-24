import logging
from collections import namedtuple

import pandas as pd
from Bio import SeqIO, Entrez
from BCBio import GFF

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'

logging.basicConfig(filename=str(snakemake.log), format=LOG_FORMAT, level=logging.INFO)

df_mash = pd.read_table(snakemake.input.mash).set_index('match', drop=False)
# ensure that Mash screen results table is sorted by identity in descending order
df_mash.sort_values('identity', ascending=False, inplace=True)

top_mash_genome = df_mash.match[0]

logging.info(f'Top reference genome by Mash screen was {top_mash_genome} with Mash result {df_mash.loc[top_mash_genome,].to_dict()}. Parsing "{snakemake.input.genbank}" for id="{top_mash_genome}"')

for rec in SeqIO.parse(snakemake.input.genbank, format='genbank'):
    if rec.id == top_mash_genome:
        logging.info(f'Found {top_mash_genome}. "{rec}".')
        logging.info(f'Writing to FASTA to "{snakemake.output.fasta}"')
        logging.info(f'GenBank to "{snakemake.output.genbank}"')
        logging.info(f'GFF to "{snakemake.output.gff}"')
        SeqIO.write([rec], snakemake.output.fasta, 'fasta')
        SeqIO.write([rec], snakemake.output.genbank, 'genbank')
        GFF.write([rec], open(snakemake.output.gff, 'w'))
        break
