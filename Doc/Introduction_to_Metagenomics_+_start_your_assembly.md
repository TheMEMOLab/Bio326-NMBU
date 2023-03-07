# Prokaryota dry lab: 01 Metagenomic Reconstruction
#### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples

Hello! Welcome to the metagenomic dry lab session! Here we will take the raw basecalled reads from the nanopore sequenator and try to reconstruct the microbial genomes that they come from.


## Filter raw reads

Your raw reads from the prokaryotic sequencing session reside in "/mnt/courses/BIO326/PROK/data/metagenomic_assembly".

There is one file named "raw_reads.fastq.gz"

Fastq is a raw read format containing a base quality score for each position along each read.
We will filter these reads using filtlong.
By specifying `--min_length 1000` and `--keep_percent 90` we keep only the reads that live up to these requirements.

Before we get started, create a directory in your work dir named metaG and copy the file containing the raw reads.

Then create a slurm-script with the following contents:

```
#!/bin/bash
#SBATCH --job-name=filtlong
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem=8G
#SBATCH --partition=normal

filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    output/samtools/mapped.sorted.bam.fastq.gz \
    | gzip \
    > output/filtlong/output.fastq.gz

```





