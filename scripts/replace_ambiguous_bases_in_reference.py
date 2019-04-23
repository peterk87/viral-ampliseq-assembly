import logging

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'
logging.basicConfig(filename=str(snakemake.log), format=LOG_FORMAT, level=logging.INFO)

good_nts = set('ACGT')

# IUPAC nucleotide code   Base
# A   Adenine
# C   Cytosine
# G   Guanine
# T (or U)    Thymine (or Uracil)
# R   A or G
# Y   C or T
# S   G or C
# W   A or T
# K   G or T
# M   A or C
# B   C or G or T
# D   A or G or T
# H   A or C or T
# V   A or C or G
# N   any base
# . or -  gap
# https://www.bioinformatics.org/sms/iupac.html

ambigous_repl = set('RYSWKMBDHV')

input_fasta = snakemake.input[0]
output_fasta = snakemake.output[0]


def replace_ambiguous(rec):
    replaced: int = 0
    mutable_seq = rec.seq.tomutable()
    for idx, nt in enumerate(rec.seq):
        nt = nt.upper()
        if nt not in good_nts and nt in ambigous_repl:
            logging.info(f'Found ambigous base in {rec.id} "{nt}" at position {idx + 1}. Replacing with "N"')
            mutable_seq[idx] = 'N'
            replaced += 1
        elif nt not in good_nts and nt not in ambigous_repl:
            logging.warning(f'Character "{nt}" at index {idx} in "{rec.id}" not in ambigous base substitution dictionary! Replacing with "N"')
            mutable_seq[idx] = 'N'
    if replaced > 0:
        logging.info(f'Found and replaced {replaced} bases in "{rec.description}".')
        rec.seq = mutable_seq.toseq()
    return rec

recs = [replace_ambiguous(r) for r in SeqIO.parse(input_fasta, format='fasta')]

SeqIO.write(recs, output_fasta, 'fasta')
