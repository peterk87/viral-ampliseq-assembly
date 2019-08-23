# Snakemake workflow: viral-ampliseq-assembly

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.4-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/peterk87/viral-ampliseq-assembly.svg?branch=master)](https://travis-ci.com/peterk87/viral-ampliseq-assembly)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3359)

[Snakemake][] workflow for analysis and assembly of viral genomes such as Classical Swine Fever Virus ([CSFV][]) from IonTorrent AmpliSeq data.

## Overview

- Preprocessing

    - Duplicate reads were removed using [Picard][]
    - Reads were trimmed with [Trimmomatic][] prior to [SPAdes][] assembly
    - BAM file stats computed using [Samtools][] (coverage depth, extent, extent per genome, # of reads mapped)

- Reference Genome Selection

    - Downloading of all Classical swine fever virus ([CSFV][]) (or FMDV, Ebola, Zika) virus genomes from [NCBI Entrez API][] using [BioPython][]
    - [Mash][] screen of deduplicated reads against all reference genomes with sketch size of 10000 and sketch k-mer size of 16, sorting by Mash screen identity to find top reference genome for read mapping and variant calling

- Read Mapping & Variant Calling

    - Read mapping with [BWA MEM][]
    - Removal of duplicate reads with [Samtools][]
    - Variant calling with [FreeBayes][]
    - [SnpEff][] was used to predict and report variant effects using reference genome annotation

- De Novo Assembly
    - [SPAdes][] de novo assembly of trimmed deduplicated reads.
    - [QUAST][] quality assessment of assemblies

- Quality Control

    - [MultiQC][] interactive report of [FastQC][], [Samtools][], [QUAST][], [SnpEff][]

- Phylogenetic Tree

    - Phylogenetic tree constructed with [IQ-TREE][] (or [Clearcut][] if a quick and dirty tree is okay)
    - Interactive HTML phylogenetic tree visualization with [PhyloCanvas][] using [shiptv][]

## Authors

* Peter Kruczkiewicz (@peterk87)

## Usage

### Step 0: Install pre-requisites

Running this workflow with [Singularity][] is recommended, but you can use [Conda][] if you prefer. The Singularity image will come with all the dependencies bundled together in a single file. 

#### Install [Singularity][] (recommended)

Follow the instructions for installing Singularity [here](https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-start)

#### Setup and activate the [Conda] environment if not using [Singularity][] (optional)

Install [Conda] if you haven't already following [these instructions](https://conda.io/en/latest/miniconda.html) and setup the [BioConda channel](https://bioconda.github.io/user/install.html#set-up-channels).

Download or `git clone` this repo

```bash
git clone https://github.com/peterk87/viral-ampliseq-assembly.git
cd viral-ampliseq-assembly
# create a conda environment named "viral-ampliseq-assembly-1.0.0"
conda env create -f environment.yml
conda activate viral-ampliseq-assembly-1.0.0
# install snakemake into this env
conda install -y snakemake
# run Snakemake on the test directory
snakemake --directory test/
```

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/peterk87/viral-ampliseq-assembly/releases).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).


### Step 2: Configure workflow

Create an analysis directory, copy and modify the example `config.yaml` and `samples.tsv` files to suit your needs. 

e.g.

```
mkdir ~/my-ampliseq-analysis
cp viral-ampliseq-assembly/config.yaml ~/my-ampliseq-analysis/
cp viral-ampliseq-assembly/samples.tsv ~/my-ampliseq-analysis/
```

Edit your `config.yaml` as needed.

Add sample entries to your `samples.tsv`

```
sample  bam_file
Sample1 bams/Sample1.bam
Sample2 bams/Sample2.bam
Sample3 bams/Sample3.bam
... <more sample entries>
```
where `bam_file` can be the relative or absolute path to a sample's BAM file.


#### [IQ-TREE] maximum-likelihood or [Clearcut] RNJ tree

In your `config.yaml` the `fast_tree` parameter controls which method (ML or RNJ) is used for phylogenetic tree construction.

If you want a quick and dirty tree, set

```yaml
fast_tree: true
```

in your `config.yaml` to generate a Relaxed Neighbor Joining (RNJ) tree.

Otherwise, if you want a high accuracy phylogenetic tree and are willing to wait for it, then set

```yaml
fast_tree: false
```

to use [IQ-TREE] to generate a maximum-likelihood phylogenetic tree with 1000 ultrafast bootstraps (UFBoot) (see [Minh et al., 2016](http://dx.doi.org/10.1093/molbev/mst024) for more info on UFBoot).

### Step 3: Execute workflow

*If you do not have [Singularity] installed then remove the `--use-singularity` flag*

Test your configuration by performing a dry-run via

    snakemake --use-singularity -n

Execute the workflow locally via

    snakemake --use-singularity --cores $N

using `$N` cores.

#### Cluster execution

*Note: You may need to install the `drmaa` Python library (`pip install drmaa`)*

You can execute the workflow on a SLURM/DRMAA cluster environment with

    snakemake --directory test --use-singularity --drmaa " -c 4 -p YourClusterQueueName --mem=4096 " -j 8 -w 60

This will run the workflow on the test data in the `test/` directory with 4 CPUs and 4G memory per job and 8 jobs at once (`-j 8`) while waiting 60 seconds for output files to appear on the shared filesystem (`-w 60`).

The cluster partition or queue to schedule jobs to is specified with `-p YourClusterQueueName`.

The above will run each rule or job with 4 CPUs and 4GB memory each, which may be way more than needed or not enough so you could create a YAML (or JSON) file to specify default and specific resource requirements for some steps:

Example `cluster-config.yaml`:

```yaml
__default__:
    cpu: 1
    partition: YourClusterQueueName
    memory: 1024
samtools_index_bam_initial:
    cpu: 32
    memory: 16384
spades_assembly:
    cpu: 32
    memory: 16384
bwa_mem:
    cpu: 32
    memory: 4096
mafft_msa:
    cpu: 32
    memory: 4096
iqtree:
    cpu: 8
    memory: 4096
snpeff:
    memory: 4096
```

With the `cluster-config.yaml`, run the workflow in a cluster environment via

    snakemake --directory test --use-singularity --cluster-config cluster-config.yaml --drmaa " -c {cluster.cpu} -p {cluster.partition} --mem={cluster.memory} " -j 8 -w 60
    

With the above command and `cluster-config.yaml`, by default, a rule or step in the workflow will only use 1 CPU and request 1G of memory, while the rules like `iqtree` or `spades_assembly` will request more CPUs and memory from the SLURM/DRMAA scheduler. 


See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

## Testing

Tests cases are in the subfolder `test`. They should be executed via continuous integration with Travis CI.

## Output

If you were to copy the files in `test` (`samples.tsv`, `bam/` and `config.yaml`) to a new directory `my-analysis-directory` and run the workflow on that directory, i.e.


```bash
snakemake --directory my-analysis-directory/ # other args
```

The contents of `my-analysis-directory` should look like:

```bash
my-analysis-directory
├── phylogeny # Phylogenetic Tree Output
│   ├── genome-metadata.tsv
│   └── tree.html
├── config.yaml # INPUT: Workflow Execution Config File 
├── qc # Quality Control Output
│   ├── multiqc.html # MultiQC report file
│   ├── fastqc # FastQC Output
│   │   ├── Sample1.html
│   │   └── Sample1_fastqc.zip
│   ├── multiqc_data
│   │   ├── [Text files]
│   └── quast # QUAST Output
│       ├── report.tex
│       ├── icarus_viewers
│       │   └── contig_size_viewer.html
│       ├── report.html
│       ├── basic_stats
│       │   ├── [QUAST PDFs]
│       ├── icarus.html
│       ├── transposed_report.tex
│       ├── quast.log
│       ├── report.pdf
│       ├── report.txt
│       ├── .snakemake_timestamp
│       ├── report.tsv
│       ├── transposed_report.tsv
│       └── transposed_report.txt
├── variant_calling # Variant Calling Output
│   ├── Sample1-filtered.vcf # Filtered variants for Sample1 in VCF format
│   ├── Sample1.vcf # Unfiltered variants for Sample1 in VCF format
│   ├── snpeff # SnpEff Output
│   │   ├── Sample1
│   │   │   ├── [SnpEff specific files]
│   │   ├── Sample1.vcf
│   │   ├── Sample1.csv
│   │   ├── Sample1.html # SnpEff report for Sample1
│   │   └── Sample1.genes.txt
│   └── Sample1-vcf.tsv # SnpEff annotated variants in a tab-delimited table
├── mapping # Read Mapping Output
│   └── Sample1 # Read mapping output and summary files for Sample1
│       ├── Sample1-extent.tsv
│       ├── Sample1-genome_extent.tsv
│       ├── Sample1-idxstats.tsv
│       ├── Sample1.bam
│       ├── Sample1-depth.tsv
│       ├── Sample1-idxstats-sorted.tsv
│       ├── Sample1-idxstats-top_mapped.txt
│       └── Sample1.bam.bai
├── bam # Input directory with Sample1 BAM file specified in config.yaml
│   └── a.bam
├── consensus # Consensus Sequence Output
│   └── Sample1.fasta # Consensus sequence for Sample1 from reference mapping and variant calling
├── logs # Log files for various tools
│   ├── <tool name>
│   │   └── Sample1.log
├── samples.tsv # INPUT: tab-delimited table with 2 fields: "sample" and "bam_file"
├── references # Reference Genomes Downloaded From NCBI
│   ├── Sample1 # Top Reference Genome
│   │   ├── reference.gff
│   │   ├── reference-no_ambig.fasta.bwt
│   │   ├── reference-no_ambig.fasta.pac
│   │   ├── reference.genbank
│   │   ├── reference-no_ambig.fasta.amb
│   │   ├── reference-no_ambig.fasta.ann
│   │   ├── reference-no_ambig.fasta
│   │   ├── reference-no_ambig.fasta.sa
│   │   ├── reference.fasta
│   │   └── reference-no_ambig.fasta.fai
│   ├── csf.msh # Mash sketch database from "csf.fasta"
│   ├── csf.genbank # CSFV genomes downloaded from NCBI in GenBank format
│   └── csf.fasta # CSFV genomes downloaded from NCBI in FASTA format
├── assembly # Assembly Output
│   ├── spades # SPAdes assembly outputs for each input sample
│   │   └── Sample1 # SPAdes assembly output for Sample1
│   │       ├── before_rr.fasta
│   │       ├── params.txt
│   │       ├── contigs.paths
│   │       ├── input_dataset.yaml
│   │       ├── <SPAdes specific output directories>
│   │       ├── scaffolds.paths
│   │       ├── contigs.fasta
│   │       ├── spades.log
│   │       ├── assembly_graph.fastg
│   │       ├── dataset.info
│   │       ├── scaffolds.fasta
│   │       └── assembly_graph_with_scaffolds.gfa
│   └── spades-Sample1.fasta
├── benchmarks # Benchmark runtime info for tools in workflow
│   ├── <benchmark tab-delimited files for various tools in workflow>
├── msa # Multiple sequence alignment (MSA) output and IQ-TREE/Clearcut phylogenetic tree
│   ├── alignment.fasta
│   ├── samples-pre-aln.fasta
│   └── alignment.fasta.treefile
└── preprocess # Preprocessing Output of Input BAM Files 
    ├── samtools # Initial BAM file stats output
    │   ├── depth
    │   │   ├── Sample1-extent.tsv
    │   │   ├── Sample1-genome_extent.tsv
    │   │   └── Sample1.tsv
    │   ├── flagstat
    │   │   └── Sample1.flagstat
    │   ├── index
    │   │   └── Sample1.done
    │   └── idxstats
    │       ├── Sample1-top_mapped.txt
    │       ├── Sample1.tsv
    │       └── Sample1-sorted.tsv
    ├── fastqs # Deduplicated reads in FASTQ format
    │   └── Sample1.fastq
    ├── mash # Mash Screen results
    │   ├── Sample1-screen_references-sorted.tsv
    │   └── Sample1-screen_references.tsv
    ├── trimmed_fastqs # Trimmomatic trimmed reads
    │   └── Sample1.fastq
    └── dedup # Deduplicated BAM files
        ├── Sample1.bam
        ├── Sample1.metrics.txt
        └── Sample1.bam.bai
```


[BioPython]: https://biopython.org/
[BWA MEM]: https://github.com/lh3/bwa
[Clearcut]: http://bioinformatics.hungry.com/clearcut/
[Conda]: https://conda.io/en/latest/
[CSFV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=11096
[FastQC]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[FreeBayes]: https://github.com/ekg/freebayes
[IQ-TREE]: http://www.iqtree.org/
[Mash]: https://mash.readthedocs.io/en/latest/
[MultiQC]: https://multiqc.info/
[NCBI Entrez API]: https://www.ncbi.nlm.nih.gov/books/NBK25501/
[PhyloCanvas]: http://phylocanvas.org/
[Picard]: https://broadinstitute.github.io/picard/
[QUAST]: http://quast.sourceforge.net/quast.html
[Samtools]: https://samtools.github.io/
[shiptv]: https://github.com/peterk87/shiptv
[Singularity]: https://sylabs.io/docs/
[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[SnpEff]: http://snpeff.sourceforge.net/SnpEff.html
[SPAdes]: http://cab.spbu.ru/software/spades/
[Trimmomatic]: http://www.usadellab.org/cms/?page=trimmomatic
