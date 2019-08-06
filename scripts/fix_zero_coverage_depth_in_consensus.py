import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sample_name = snakemake.wildcards.sample
consensus_fasta = snakemake.input.fasta
depth_file = snakemake.input.depth
unmapped_char = snakemake.params.unmapped_char

seq_rec = SeqIO.read(consensus_fasta, format='fasta')
df = pd.read_table(depth_file, names='genome position coverage'.split())
mutable_seq = seq_rec.seq.tomutable()
for position in df[df.coverage == 0].position:
    mutable_seq[position - 1] = unmapped_char

consensus_seq_rec = SeqRecord(id=sample_name,
                              seq=mutable_seq.toseq(),
                              description=f'consensus_from_ref="{seq_rec.description}"')
SeqIO.write([consensus_seq_rec], snakemake.output.fasta, 'fasta')
