"""
Variant calling rules
"""

rule freebayes:
    """
    Variant calling with FreeBayes
    """
    input:
        ref='references/{sample}/reference-no_ambig.fasta',
        samples='mapping/{sample}/{sample}.bam'
    output:
        'variant_calling/{sample}.vcf'
    log:
        'logs/freebayes/{sample}.log'
    params:
        extra='-p 2 -P 0 -C ' + str(config['freebayes'].get('minobs', 5)) + \
              ' --min-repeat-entropy 1.5 --strict-vcf ' + \
              ' -q ' + str(config['freebayes'].get('basequal', 13)) + \
              ' -m ' + str(config['freebayes'].get('mapqual', 60)) + \
              ' --min-coverage ' + str(config['freebayes'].get('mincov', 5)) + \
              ' -F 0.05'
    threads: 1
    wrapper:
        '0.27.1/bio/freebayes'


rule filter_vcf:
    """
    Filter VCF using bcftools
    """
    input:
        'variant_calling/{sample}.vcf'
    output:
        'variant_calling/{sample}-filtered.vcf'
    params:
        minqual=config['vcf_filtering'].get('minqual', 100),
        mincov=config['freebayes'].get('mincov', 5),
        alt_ref_ratio=config['vcf_filtering'].get('alt_ref_ratio', 1.1)
    shell:
        '''
        bcftools view \
        --include 'QUAL>={params.minqual} && FMT/DP>={params.mincov} && (FMT/AO)/(FMT/RO)>={params.alt_ref_ratio}' \
        {input} > {output}
        '''


rule snpeff_build:
    """
    Build a snpEff genome database for the reference genome for the current 
    sample.
    """
    input:
        gff='references/{sample}/reference.gff',
        fasta='references/{sample}/reference-no_ambig.fasta'
    output:
        config='variant_calling/snpeff/{sample}/snpeff.config'
    log:
        'logs/snpeff_build/{sample}.log'
    shell:
        '''
        touch {output.config}
        touch {log}
        INPUT_GFF=$(realpath {input.gff})
        INPUT_FASTA=$(realpath {input.fasta})
        OUTPUT_CONFIG=$(realpath {output.config})
        LOGFILE=$(realpath {log})

        BINDIR=$(dirname $(which snpEff))
        WORKDIR=$(dirname {output.config})
        ORIGINAL_CONFIG=$BINDIR/../etc/snpeff.config
        test -r $ORIGINAL_CONFIG

        # Colours!
        bold() {{ echo -e "\\e[1m$@\\e[0m"; }}
        yellow() {{ echo -e "\\e[33m$@\\e[0m"; }}
        bold_yellow() {{ bold $(yellow $@); }}
        green() {{ echo -e "\\e[32m$@\\e[0m"; }}
        bold_green() {{ bold $(green $@); }}
        red() {{ echo -e "\\e[31m$@\\e[0m"; }}
        blue() {{ echo -e "\\e[34m$@\\e[0m"; }}
        bold_blue() {{ bold $(blue $@); }}
        purple() {{ echo -e "\\e[35m$@\\e[0m"; }}
        bold_purple() {{ bold $(purple $@); }}

        bold_yellow "Starting GFF and config file preprocessing" | tee $LOGFILE
        green "Copying $(blue $ORIGINAL_CONFIG) to $(purple $OUTPUT_CONFIG)" | tee -a $LOGFILE
        cp $ORIGINAL_CONFIG $OUTPUT_CONFIG
        echo -e "Get reference genome accession from GFF file" | tee -a $LOGFILE
        ACCESSION=$(awk '/^#/{{f=1;next}} f{{ print $1 }}' {input.gff} | uniq -)
        echo -e "Ref genome accession is $(red $ACCESSION)" | tee -a $LOGFILE
        yellow "Adding entry for $(red $ACCESSION) to $(purple $OUTPUT_CONFIG)" | tee -a $LOGFILE
        echo -e "ref.genome : The Reference" | tee -a $OUTPUT_CONFIG
        echo -e "\\tref.chromosome : $ACCESSION" | tee -a $OUTPUT_CONFIG
        echo -e "\\tref.$ACCESSION.codonTable : Standard" | tee -a $OUTPUT_CONFIG

        blue "Printing last 3 lines of $(purple $OUTPUT_CONFIG)" | tee -a $LOGFILE
        tail -n3 $OUTPUT_CONFIG | tee -a $LOGFILE
        pushd $WORKDIR > /dev/null
        bold_yellow "Changed dir to snpEff build work dir $(bold_blue $PWD)" | tee -a $LOGFILE
        mkdir -p ref
        echo "Copying reference genome GFF and FASTA to ref/" | tee -a $LOGFILE
        echo "Stripping away UTR features since they cause issues for snpEff" | tee -a $LOGFILE
        grep -Pv "^${{ACCESSION}}\\tfeature\\t[^\\s]*UTR\\t" $INPUT_GFF > ref/genes.gff
        green "Original GFF lineno=$(bold_purple $(wc -l $INPUT_GFF))" | tee -a $LOGFILE
        green "snpEff GFF   lineno=$(bold_purple $(wc -l $(realpath ref/genes.gff)))" | tee -a $LOGFILE
        cp $INPUT_FASTA ref/sequences.fa
        bold_yellow "Finished GFF and config file preprocessing" | tee -a $LOGFILE

        bold_green "Running snpEff build..." | tee -a $LOGFILE
        snpEff build -v -c $OUTPUT_CONFIG -dataDir . -gff3 ref | tee -a $LOGFILE
        # go back to original work dir
        popd > /dev/null
        '''


rule snpeff:
    """
    Find the effects of SNPs in the sample genome compared to an annotated 
    reference genome.
    """
    input:
        vcf='variant_calling/{sample}-filtered.vcf',
        snpeff_config='variant_calling/snpeff/{sample}/snpeff.config'
    output:
        vcf='variant_calling/snpeff/{sample}.vcf',
        htmlstats=report('variant_calling/snpeff/{sample}.html',
                         caption='../report/results/snpeff.rst',
                         category='Variant Calling'),
        csvstats='variant_calling/snpeff/{sample}.csv'
    params:
        extra='-Xmx4g' # optional parameters (e.g., max memory 4g)
    shell:
        '''
        touch {output.vcf}
        touch {output.htmlstats}
        touch {output.csvstats}
        INPUT_VCF=$(realpath {input.vcf})
        INPUT_SNPEFF_CONFIG=$(realpath {input.snpeff_config})
        OUTPUT_VCF=$(realpath {output.vcf})
        OUTPUT_HTMLSTATS=$(realpath {output.htmlstats})
        OUTPUT_CSVSTATS=$(realpath {output.csvstats})
        WORKDIR=$(dirname {input.snpeff_config})
        pushd $WORKDIR > /dev/null
        snpEff ann {params.extra} \
        -noLog -csvStats $OUTPUT_CSVSTATS \
        -htmlStats $OUTPUT_HTMLSTATS \
        -c $INPUT_SNPEFF_CONFIG \
        -dataDir . \
        ref $INPUT_VCF > $OUTPUT_VCF
        popd > /dev/null
        '''


rule vcf_to_tab:
    input:
        gff='references/{sample}/reference.gff',
        fasta='references/{sample}/reference-no_ambig.fasta',
        vcf='variant_calling/snpeff/{sample}.vcf'
    output:
        report('variant_calling/{sample}-vcf.tsv', 
               caption='../report/results/variants_table.rst', 
               category='Variant Calling')
    log:
        'logs/vcf_to_tab/{sample}.log'
    shell:
        'snippy-vcf_to_tab -gff {input.gff} --ref {input.fasta} --vcf {input.vcf} > {output} 2> {log}'
