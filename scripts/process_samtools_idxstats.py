import pandas as pd

columns = 'genome length n_mapped n_unmapped'.split()
df = pd.read_table(snakemake.input[0], 
                   names=columns)\
       .set_index('genome', drop=False)
df.sort_values(['n_mapped', 'length'], ascending=[False,False], inplace=True)
# write sorted table of mapping stats with headers to a file
df.to_csv(snakemake.output.sorted, sep='\t', index=False)
# write accession of topped mapped genome to a file
with open(snakemake.output.top_mapped, 'wt') as f:
    f.write(df.genome[0])
