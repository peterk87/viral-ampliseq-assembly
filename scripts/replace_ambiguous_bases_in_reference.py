import logging

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'
logging.basicConfig(filename=str(snakemake.log), format=LOG_FORMAT, level=logging.INFO)

good_nts = set('ACGTacgtNn')

input_fasta = snakemake.input[0]

def replace_ambiguous(rec, repl='N'):
    replaced: int = 0
    mutable_seq = rec.seq.tomutable()
    for idx, nt in enumerate(rec.seq):
        if nt not in good_nts:
            logging.info(f'Found ambigous base in {rec.id} "{nt}" at position {idx + 1}. Replacing with "{repl}"')
            mutable_seq[idx] = repl
            replaced += 1
    if replaced > 0:
        logging.info(f'Found and replaced {replaced} bases in "{rec.description}" with "{repl}"')
        rec.seq = mutable_seq.toseq()
    return rec


recs = [replace_ambiguous(r) for r in SeqIO.parse(input_fasta, format='fasta')]

SeqIO.write(recs, snakemake.output[0], 'fasta')
