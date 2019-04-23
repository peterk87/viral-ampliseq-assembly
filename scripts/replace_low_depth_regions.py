import logging

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'
logging.basicConfig(filename=str(snakemake.log), format=LOG_FORMAT, level=logging.INFO)

input_fasta = snakemake.input.fasta
input_depth = snakemake.input.depth
output_fasta = snakemake.output[0]

params_low_coverage = snakemake.params.low_coverage

def replace_low_depth(rec, df_depth):
    replaced: int = 0
    mutable_seq = rec.seq.tomutable()
    for idx, nt in enumerate(rec.seq):
        coverage = df_depth.query(f'genome == "{rec.id}" and position == {idx + 1}').coverage[0]
        if coverage == params_low_coverage:
            mutable_seq[idx] = 'N'
            replaced += 1
    if replaced > 0:
        logging.info(f'Found and replaced {replaced} bases in "{rec.description}".')
        rec.seq = mutable_seq.toseq()
    return rec


df = pd.read_table(input_depth, names='genome position coverage'.split()).set_index(['genome', 'position'])
recs = [replace_ambiguous(r, df) for r in SeqIO.parse(input_fasta, format='fasta')]

SeqIO.write(recs, output_fasta, 'fasta')
