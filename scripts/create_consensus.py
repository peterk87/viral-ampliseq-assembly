from typing import List, Dict

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# INPUT
input_ref_fasta = snakemake.input.ref
input_vcf_tsv = snakemake.input.vcf_tsv
input_depth_tsv = snakemake.input.depth
# OUTPUT
output_consensus_fasta = snakemake.output.fasta
# WILDCARDS
sample_name = snakemake.wildcards.sample
# PARAMS
min_coverage: int = snakemake.params.get('min_coverage', 5)
unmapped_char: str = snakemake.params.get('unmapped_char', '-')

assert len(unmapped_char) == 1, '"unmapped_char" must be a str of length of 1, e.g. "-"/single dash character.'


def replace_low_depth_positions(input_ref_fasta: str, 
                                input_depth_tsv: str, 
                                min_coverage: int = 5, 
                                unmapped_char: str = '-') -> SeqRecord:
    """Replace low coverage depth positions with a specified character.

    Args:
        input_ref_fasta: FASTA file path
        input_depth_tsv: Samtools depth tab-delimited file path
        min_coverage: Minimum coverage depth threshold, below which positions are changed to the `unmapped_char` (default: 5)
        unmapped_char: Character to substitute at low coverage positions (default: "-")

    Returns:
        SeqRecord with low coverage depth positions substituted with `unmapped_char` 
    """
    seq_rec = SeqIO.read(input_ref_fasta, format='fasta')
    df = pd.read_table(input_depth_tsv, names='genome position coverage'.split())
    low_coverage_positions: pd.Series = df[df.coverage < min_coverage].position
    if low_coverage_positions.size == 0:
        return seq_rec
    mutable_seq = seq_rec.seq.tomutable()
    for position in low_coverage_positions:
        mutable_seq[position - 1] = unmapped_char
    seq_rec.seq = mutable_seq.toseq()
    return seq_rec


def parse_vcf_evidence(s: str) -> Dict[str, int]:
    """Parse VCF evidence column value into dict of variant sequence to frequency.
    
    >>> parse_vcf_evidence('CGA,CCGA:48,38 CTGA:0')
    {'CGA': 48, 'CCGA': 38, 'CTGA': 0}
    
    >>> parse_vcf_evidence('CCGA:48 CTGA:0')
    {'CCGA': 48, 'CTGA': 0}
    
    """
    alt_evid, ref_evid = s.split(' ')
    ref, ref_count = ref_evid.split(':')
    alts, counts = alt_evid.split(':')
    alt_count = {a: int(c) for a, c in zip(alts.split(','), counts.split(','))}
    alt_count.update({ref: int(ref_count)})
    return alt_count


def top_var(ref_var: str, var_to_count: Dict[str, int]) -> str:
    """Pick the most frequently observed variant with same length as the reference variant.
    
    >>> top_var('CTGA', {'CGA': 48, 'CCGA': 38, 'CTGA': 0})
    'CCGA'
    """
    max_count = -1
    max_var = ''
    ref_len = len(ref_var)
    for var, count in var_to_count.items():
        if count > max_count and len(var) == ref_len:
            max_count = count
            max_var = var
    return max_var


def consensus_segment(seq: str, 
                      curr_position: int, 
                      ref_var: str, 
                      top_variant: str, 
                      prev_position: int = 0) -> (str, int):
    return (seq[prev_position:(curr_position - 1)] + top_variant), (curr_position + len(ref_var) - 1)


def create_cons_seq(seq: str, df_vcf: pd.DataFrame) -> str:
    """Create consensus sequence given a reference sequence and a table of a Snippy vcf_to_tab

    For each variant, the ALT or REF allele that has the greatest support and 
    the same length as the REF allele is selected. Given the high error rate for
    IonTorrent data, shorter variants are ignored as they are likely due to 
    homopolymer sequencing error. 

    Args:
        seq: Reference sequence
        df_vcf: DataFrame of Snippy vcf_to_tab output from VCF file
    
    Returns:
        Consensus sequence favouring variants of the same length as the REF variant.
    """
    segs = []
    prev_position = 0
    for idx, curr_var in df_vcf.iterrows():
        segment, prev_position = consensus_segment(seq=seq, 
                                                   curr_position=curr_var.POS, 
                                                   ref_var=curr_var.REF, 
                                                   top_variant=top_var(curr_var.REF, parse_vcf_evidence(curr_var.EVIDENCE)), 
                                                   prev_position=prev_position)
        segs.append(segment)
    # append the rest of the reference sequence
    segs.append(seq[prev_position:])
    return ''.join(segs)



ref_seq_record: SeqRecord = replace_low_depth_positions(input_ref_fasta=input_ref_fasta, 
                                                        input_depth_tsv=input_depth_tsv,
                                                        min_coverage=min_coverage,
                                                        unmapped_char=unmapped_char)
df_vcf_tsv = pd.read_table(input_vcf_tsv)
consensus_seq: str = create_cons_seq(str(ref_seq_record.seq), df_vcf_tsv)

with open(output_consensus_fasta, 'w') as f:
    f.write(f'>{sample_name} consensus_from_ref="{ref_seq_record.description}"\n{consensus_seq}\n')
