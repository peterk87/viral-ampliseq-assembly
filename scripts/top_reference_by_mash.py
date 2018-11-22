import logging
from collections import namedtuple

import pandas as pd
from Bio import SeqIO, Entrez

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)


df_mash = pd.read_table(snakemake.input.mash).set_index('match', drop=False)
# ensure that Mash screen results table is sorted by identity in descending order
df_mash.sort_values('identity', ascending=False, inplace=True)

top_mash_genome = df_mash.match[0]

logging.info(f'Top reference genome by Mash screen was {top_mash_genome} with Mash result {df_mash.loc[top_mash_genome,].to_dict()}. Parsing "{snakemake.input.fasta}" for id="{top_mash_genome}"')

for rec in SeqIO.parse(snakemake.input.fasta, format='fasta'):
    if rec.id == top_mash_genome:
        logging.info(f'Found {top_mash_genome}. {rec}. Writing to {snakemake.output[0]}.')
        SeqIO.write([rec], snakemake.output[0], 'fasta')
        break
