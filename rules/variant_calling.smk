
rule freebayes:
    """
    Variant calling with FreeBayes
    """
    input:
        ref='mapping/{sample}/reference.fasta',
        samples='mapping/{sample}/{sample}.bam'
    output:
        'variant_calling/{sample}.vcf'
    log:
        'logs/freebayes/{sample}.log'
    params:
        extra='-p 2 -P 0 -C ' + str(config['freebayes'].get('mincov', 10)) + \
              ' --min-repeat-entropy 1.5 --strict-vcf ' + \
              ' -q ' + str(config['freebayes'].get('basequal', 13)) + \
              ' -m ' + str(config['freebayes'].get('mapqual', 60)) + \
              ' --min-coverage ' + str(config['freebayes'].get('mincov', 10)) + \
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
        minqual=80, # Minimum read mapping quality to consider
        mincov=10, # Minimum site depth to for calling alleles
        minfrac=0 # Minumum proportion for variant evidence (0=AUTO)
    conda:
        '../envs/snippy.yaml'
    shell:
        '''
        bcftools view \
        --include 'QUAL>={params.minqual} && FMT/DP>={params.mincov} && (FMT/AO)/(FMT/DP)>={params.minfrac}' \
        {input} > {output}
        '''


rule snpeff_build:
    """
    Build a snpEff genome database for the reference genome for the current 
    sample.
    """
    input:
        gff='mapping/{sample}/reference.gff',
        fasta='mapping/{sample}/reference.fasta'
    output:
        config='variant_calling/snpeff/{sample}/snpeff.config'
    log:
        'logs/snpeff_build/{sample}.log'
    conda:
        '../envs/snippy.yaml'
    shell:
        '''
        INPUT_GFF=$(realpath {input.gff})
        INPUT_FASTA=$(realpath {input.fasta})
        touch {output.config}
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
        
        bold_yellow "Starting GFF and config file preprocessing" #> $LOGFILE
        
        green "Copying $(blue $ORIGINAL_CONFIG) to $(purple $OUTPUT_CONFIG)" #>> $LOGFILE
        cp $ORIGINAL_CONFIG $OUTPUT_CONFIG
        
        echo -e "Get reference genome accession from GFF file" #>> $LOGFILE
        ACCESSION=$(awk '/^#/{{f=1;next}} f{{ print $1 }}' {input.gff} | uniq -)
        echo -e "Ref genome accession is $(red $ACCESSION)" #>> $LOGFILE
        yellow "Adding entry for $(red $ACCESSION) to $(purple $OUTPUT_CONFIG)" #>> $LOGFILE
        echo -e "ref.genome : The Reference" #>> $OUTPUT_CONFIG
        echo -e "\\tref.chromosome : $ACCESSION" #>> $OUTPUT_CONFIG
        echo -e "\\tref.$ACCESSION.codonTable : Standard" #>> $OUTPUT_CONFIG

        blue "Printing last 3 lines of $(purple $OUTPUT_CONFIG)" #>> $LOGFILE
        tail -n3 $OUTPUT_CONFIG #>> $LOGFILE
        pushd $WORKDIR
        bold_yellow "Changed dir to snpEff build work dir $(bold_blue $PWD)" #>> $LOGFILE
        mkdir -p ref
        echo "Copying reference genome GFF and FASTA to ref/" #>> $LOGFILE
        echo "Stripping away UTR features since they cause issues for snpEff" #>> $LOGFILE
        grep -Pv "^${{ACCESSION}}\\tfeature\\t[^\\s]*UTR\\t" $INPUT_GFF > ref/genes.gff
        green "Original GFF lineno=$(bold_purple $(wc -l $INPUT_GFF))" #>> $LOGFILE
        green "snpEff GFF   lineno=$(bold_purple $(wc -l $(realpath ref/genes.gff)))" #>> $LOGFILE
        cp $INPUT_FASTA ref/sequences.fa
        bold_yellow "Finished GFF and config file preprocessing" #>> $LOGFILE
        
        bold_green "Running snpEff build..." #>> $LOGFILE
        snpEff build -v -c $OUTPUT_CONFIG -dataDir . -gff3 ref #&>> $LOGFILE
        # go back to original work dir
        popd
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
        htmlstats='variant_calling/snpeff/{sample}.html',
        csvstats='variant_calling/snpeff/{sample}.csv'
    params:
        extra='-Xmx4g' # optional parameters (e.g., max memory 4g)
    conda: 
        '../envs/snippy.yaml'
    shell:
        '''
        INPUT_VCF=$(realpath {input.vcf})
        INPUT_SNPEFF_CONFIG=$(realpath {input.snpeff_config})
        OUTPUT_VCF=$(realpath {output.vcf})
        OUTPUT_HTMLSTATS=$(realpath {output.htmlstats})
        OUTPUT_CSVSTATS=$(realpath {output.csvstats})
        WORKDIR=$(dirname {input.snpeff_config})
        pushd $WORKDIR
        snpEff ann {params.extra} \
        -noLog -csvStats $OUTPUT_CSVSTATS \
        -htmlStats $OUTPUT_HTMLSTATS \
        -c $INPUT_SNPEFF_CONFIG \
        -dataDir . \
        ref $INPUT_VCF > $OUTPUT_VCF
        popd
        '''


rule vcf_to_tab:
    input:
        gff='mapping/{sample}/reference.gff',
        fasta='mapping/{sample}/reference.fasta',
        vcf='variant_calling/snpeff/{sample}.vcf'
    output:
        'variant_calling/{sample}-vcf.tsv'
    log:
        'logs/vcf_to_tab/{sample}.log'
    conda:
        '../envs/snippy.yaml'
    shell:
        'snippy-vcf_to_tab -gff {input.gff} --ref {input.fasta} --vcf {input.vcf} > {output} 2> {log}'
