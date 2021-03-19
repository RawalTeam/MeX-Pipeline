# MeX Pipeline

### Pre-requisites
1. Conda (<a href="https://docs.conda.io/en/latest/miniconda.html">Miniconda</a> or <a href="https://www.anaconda.com/products/individual">Anaconda</a>)
2. Linux
3. Git

### Getting Started
Creating conda environment
```bash
git clone https://github.com/RawalTeam/MeX-Pipeline.git
cd  Mex-Pipeline
conda env create -f envs/mex.yaml --name mex
```

Installing additional external dependencies
```bash
conda activate mex
python install_deps.py
```

Running MeX Pipeline
```bash
conda activate mex
python mex.py --fq1 path/of/fastq1 --genome path/of/fasta --te path/of/fasta --outdir path/of/output_folder -p 2
```

Help
```bash
conda activate mex
python mex.py -h
```

### Components of MeX Pipeline
* <a href="https://github.com/OpenGene/fastp">FASTp</a>\
A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with a multithreading supported to afford high performance.

* <a href="https://github.com/s-andrews/FastQC">FASTQc</a>\
FastQC is a program designed to spot potential problems in high througput sequencing datasets. It runs a set of analyses on one or more raw sequence files in fastq or bam format and produces a report which summarises the results.

* <a href="https://github.com/bergmanlab/ngs_te_mapper">ngs_te_mapper</a>\
ngs_te_mapper is an R implementation of the method for detecting transposable element (TE) insertions from next-generation sequencing (NGS).
  
* <a href="https://www.ensembl.org/vep">Ensembl Variant Effect Predictor (VEP)</a>\
VEP determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.

### Inputs
#### Required
1. FASTq file 1 (--fq1, -1)\
Either the Read1 FASTQ file from a paired-end sequencing, or the FASTQ file from an unpaired sequencing.
   
2. Genome FASTA file (--genome, -g)\
The genome sequence of the reference genome in FASTA format.

3. TE FASTA file (--te, -te)\
A FASTA file containing a consensus sequence for each family.
   
#### Optional
1. FASTq file 2 (--fq2, -2)\
The Read2 FASTQ file from a paired-end sequencing run.
   
### Outputs
    --- /path/of/outdir
        |_ logs* (various log files)
        |_ outputs
            |_ fastp*
            |_ fastqc*
            |_ ngs_te_mapper*
            |_ vep* (created only if non-reference TE found)
        |_ config.json (internal configuration file)
        |_ Snakefile (snakemake file)
        |_ workflow.html (snakemake report)
    
    * Is a directory

### Usages
```
usage: mex.py [-h] -1 FQ1 -g GENOME -te TE -O OUTDIR [-2 FQ2] [-p PROCESSES]
              [--force]

optional arguments:
  -h, --help            show this help message and exit
  -1 FQ1, --fq1 FQ1     FASTQ Read 1 *
  -g GENOME, --genome GENOME
                        Genome FASTA *
  -te TE, --te TE       TE FASTA *
  -O OUTDIR, --outdir OUTDIR
                        Output Directory *
  -2 FQ2, --fq2 FQ2     FASTQ Read 2
  -p PROCESSES, --processes PROCESSES
                        Number of processes for multiprocessing
  --force               Rerun entire MeX pipeline

```