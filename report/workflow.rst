viral-ampliseq-assembly (version 1.0.0): Analysis and assembly of viral AmpliSeq IonTorrent sequence data. 


**Preprocessing**

- Duplicate reads were removed using Picard_
{% if snakemake.config['trim'] %}
- Reads were trimmed with Trimmomatic_ prior to SPAdes_ assembly
{% endif %}
- BAM file stats computed using Samtools_ (coverage depth, extent, extent per genome, # of reads mapped)


**Reference Genome Selection**

- Downloading of all {{ snakemake.config['organism']|upper }} virus genomes from `NCBI Entrez API`_ using BioPython_ 
- Mash_ screen of deduplicated reads against all reference genomes with sketch size of {{ snakemake.config['mash']['n_sketches'] }} and sketch k-mer size of {{ snakemake.config['mash']['kmer_size'] }}, sorting by Mash screen identity to find top reference genome for read mapping and variant calling


**Read Mapping & Variant Calling**

- Read mapping with `BWA MEM`_
- Removal of duplicate reads with Samtools_
- Variant calling with FreeBayes_
    - Minimum coverage of {{ snakemake.config['freebayes']['mincov'] }}X
    - Minimum mapping quality of {{ snakemake.config['freebayes']['mapqual'] }}
    - Minimum PHRED score base quality of {{ snakemake.config['freebayes']['basequal'] }}
- SnpEff_ was used to predict and report variant effects using reference genome annotation


**De Novo Assembly**

- SPAdes_ de novo assembly of {% if snakemake.config['trim'] %} trimmed {% endif %} deduplicated reads.
- QUAST_ quality assessment of assemblies


**Quality Control**

MultiQC_ interactive report of FastQC_, Samtools_, QUAST_, SnpEff_


.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _SPAdes: http://cab.spbu.ru/software/spades/
.. _Mash: https://mash.readthedocs.io/en/latest/
.. _NCBI Entrez API: https://www.ncbi.nlm.nih.gov/books/NBK25501/
.. _BioPython: https://biopython.org/
.. _FreeBayes: https://github.com/ekg/freebayes
.. _QUAST: http://quast.sourceforge.net/
.. _BWA MEM: http://bio-bwa.sourceforge.net/
.. _Picard: https://broadinstitute.github.io/picard
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Samtools: http://samtools.sourceforge.net/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
