rule samtools_index_bam:
    input:
        get_bam_file
    output:
        'preprocess/samtools/index/{sample}.done'
    threads: config['samtools']['threads']
    conda:
        '../envs/bwa.yaml'
    benchmark:
        'benchmark/samtools_index_bam/{sample}.tsv'
    shell:
        """
        samtools index -@ {threads} {input}
        touch {output}
        """


rule samtools_flagstat:
    input:
        get_bam_file
    output:
        'preprocess/samtools/flagstat/{sample}.flagstat'
    threads: config['samtools']['threads']
    conda:
        '../envs/bwa.yaml'
    shell:
        """
        samtools flagstat -@ {threads} {input} > {output}
        """


rule samtools_idxstats:
    """
    Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index. 
    The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 

    """
    input:
        bam=get_bam_file,
        bai_done='preprocess/samtools/index/{sample}.done'
    output:
        'preprocess/samtools/idxstats/{sample}.tsv'
    conda:
        '../envs/bwa.yaml'
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """


rule process_idxstat:
    input:
        'preprocess/samtools/idxstats/{sample}.tsv'
    output:
        sorted='preprocess/samtools/idxstats/{sample}-sorted.tsv',
        top_mapped='preprocess/samtools/idxstats/{sample}-top_mapped.txt'
    run:
        import pandas as pd

        df = pd.read_table(snakemake.input, names='genome length n_mapped n_unmapped'.split()).set_index('genome', drop=False)
        df.sort_values(['n_mapped', 'length'], ascending=[False,False], inplace=True)
        # write sorted table of mapping stats with headers to a file
        df.to_csv(snakemake.output.sorted, sep='\t', index=False)
        # write accession of topped mapped genome to a file
        with open(snakemake.output.top_mapped, 'wt') as f:
            f.write(df.genome[0])


rule samtools_depth:
    """
    Calculate depth for all positions including those with zero coverage and no
    limit on coverage.
    """
    input:
        get_bam_file
    output:
        'preprocess/samtools/depth/{sample}.tsv'
    conda:
        '../envs/bwa.yaml'
    benchmark:
        'benchmark/samtools_depth/{sample}.tsv'
    shell:
        """
        samtools depth -a -d 0 {input} > {output}
        """


rule process_samtools_depth:
    input:
        'preprocess/samtools/depth/{sample}.tsv'
    output:
        genome_extent='preprocess/samtools/depth/{sample}-genome_extent.tsv',
        extent='preprocess/samtools/depth/{sample}-extent.tsv'
    run:
        import pandas as pd

        df = pd.read_table(input[0], names='genome position coverage'.split())
        df = df.set_index('genome', drop=False)
        df_depth = df.copy()
        df_depth['is_mapped'] = df_depth['coverage'] > 0
        g = df_depth.drop(columns=['position', ]).groupby('genome')
        df_mapped_summary = pd.DataFrame(dict(mapped_positions=g.is_mapped.sum(),
                                  mean_coverage=g.coverage.mean()))
        df_mapped_summary.sort_values('mapped_positions', ascending=False, inplace=True)
        df_mapped_summary.to_csv(output.genome_extent, sep='\t')
    
        position_grouped = df.drop(columns='genome').groupby('position')
        df_position_max = position_grouped.max()
        df_position_max.to_csv(output.extent, sep='\t')

rule bam_to_fastq:
    input:
        get_bam_file
    output:
        'preprocess/fastqs/{sample}.fastq'
    conda:
        '../envs/bwa.yaml'
    shell:
        "samtools bam2fq {input} > {output}"
