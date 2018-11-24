import logging
from collections import namedtuple

import pandas as pd
from Bio import SeqIO, Entrez

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'
logging.basicConfig(filename=str(snakemake.log), format=LOG_FORMAT, level=logging.INFO)


# Organism corresponds to the config organism parameter
# fmd, csf, ebola, zika
NcbiOrganismInfo = namedtuple('NcbiOrganismInfo', ['taxid', 'tax_rank', 'min_genome_size', 'max_genome_size'])

organism_query_info = dict(
    csf=NcbiOrganismInfo(taxid=11096,
                         tax_rank='species',
                         min_genome_size=11000,
                         max_genome_size=13000),
    fmd=NcbiOrganismInfo(taxid=12110,
                         tax_rank='species',
                         min_genome_size=7000,
                         max_genome_size=8500),
    ebola=NcbiOrganismInfo(taxid=186536,
                           tax_rank='genus',
                           min_genome_size=18000,
                           max_genome_size=20000),
    zika=NcbiOrganismInfo(taxid=64320,
                          tax_rank='species',
                          min_genome_size=10000,
                          max_genome_size=11500))

logging.error(f'snakemake log {snakemake.log} {snakemake.log.__dict__}')

organism = snakemake.params.get('organism', None)
if organism is None or organism not in organism_query_info:
    raise Exception('An organism [one of fmd, csf, ebola, zika] needs to be specified!')


Entrez.email = snakemake.params.get('email', 'my.email@canada.ca')
logging.info(f'NCBI Entrez API needs an email. Using "{Entrez.email}".')

def get_ncbi_query(organism: str) -> str:
    org_info = organism_query_info[organism]
    return (f'txid{org_info.taxid}[Organism:noexp] AND biomol_genomic[PROP] '
            f'AND {org_info.min_genome_size}:{org_info.max_genome_size}[Sequence Length] '
            f'AND "complete genome"[Title]')


ncbi_query_string = get_ncbi_query(organism)
logging.info(f'Querying NCBI Nucleotide database for organism "{organism}" with query string: {ncbi_query_string}')
with Entrez.esearch(db='nucleotide', 
                    term=ncbi_query_string, 
                    retmax=10000) as esearch_handle:
    esearch_result = Entrez.read(esearch_handle)

assert 'IdList' in esearch_result
assert len(esearch_result['IdList']) > 0

logging.info(f'Fetching {len(esearch_result["IdList"])} for organism "{organism}" reference genome sequences from NCBI Nucleotide database.')
with Entrez.efetch(db='nucleotide', 
                   id=esearch_result['IdList'], 
                   rettype='gb', 
                   retmode='text') as efetch_handle:
    recs = [x for x in SeqIO.parse(efetch_handle, 'genbank')]
logging.info(f'Fetched {len(recs)} sequence records for reference genomes from NCBI Nucleotide database.')
SeqIO.write(sequences=recs,
            handle=snakemake.output.fasta,
            format='fasta')

SeqIO.write(sequences=recs,
            handle=snakemake.output.genbank,
            format='genbank')

logging.info(f'Output reference genome sequences to FASTA file "{snakemake.output.fasta}" and GenBank file "{snakemake.output.genbank}"')
