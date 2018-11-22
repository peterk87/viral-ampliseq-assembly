import pandas as pd

df_idxstats = pd.read_table(snakemake.input.idxstats).set_index('genome', drop=False)

depth_columns = 'genome position coverage'.split()
df_depth = pd.read_table(snakemake.input.depth, 
                         names=depth_columns)\
             .set_index('genome', drop=False)
df_depth_summary = df_depth.copy()
df_depth_summary['is_mapped'] = df_depth_summary['coverage'] > 0
genome_grouped = df_depth_summary.drop(columns=['position', ]).groupby('genome')
df_mapped_summary = pd.DataFrame(dict(
    mapped_positions=genome_grouped.is_mapped.sum(),
    mean_coverage=genome_grouped.coverage.mean(),
    median_coverage=genome_grouped.coverage.median(),
    max_coverage=genome_grouped.coverage.max(),
    ))
df_mapped_summary = pd.merge(df_mapped_summary, 
                             df_idxstats, 
                             on='genome', 
                             how='left')\
                      .set_index('genome', drop=False)
df_mapped_summary['p_mapped_positions'] = df_mapped_summary['mapped_positions'] / df_mapped_summary['length']
df_mapped_summary.sort_values('p_mapped_positions', ascending=False, inplace=True)
df_mapped_summary.to_csv(snakemake.output.genome_extent, sep='\t')

position_grouped = df_depth.drop(columns='genome').groupby('position')
df_position_max = position_grouped.max()
df_position_max.to_csv(snakemake.output.extent, sep='\t')
