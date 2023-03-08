# Prokaryota dry lab: 01 Metagenomic Reconstruction
#### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples

Hello! Welcome to the metagenomic dry lab session! Here we will take the raw basecalled reads from the nanopore sequenator and try to reconstruct the microbial genomes that they come from.


### Filter raw reads

Your raw reads from the prokaryotic sequencing session reside in "/mnt/courses/BIO326/PROK/data/metagenomic_assembly".

There is one file named "raw_reads.fastq.gz"

Fastq is a raw read format containing a base quality score for each position along each read.
We will filter these reads using filtlong.
By specifying `--min_length 1000` and `--keep_percent 90` we keep only the reads that live up to these requirements.

Before we get started, create a directory in your work dir named metaG and copy the file containing the raw reads.

Then create a slurm-script with the following contents:

Create a file named 01_filter-filtlong.sh and submit the job with sbatch:

```
#!/bin/bash
#SBATCH --job-name=filter-filtlong
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem=8G
#SBATCH --partition=normal

INPUT="output/samtools/mapped.sorted.bam.fastq.gz"
OUTPUT="output/filtlong/output.fastq.gz"

filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    $INPUT \
    | gzip \
    > $OUTPUT

```


### Assemble reads into a draft assembly

Then we will use flye in "meta" mode. Flye builds contiguous sequences by overlapping each read.
You can read more about how to configure flye here: https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md

Create a file named 02_assemble-flye.sh and submit the job with sbatch:

```
#!/bin/bash
#SBATCH --job-name=assemble-flye
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G
#SBATCH --partition=normal

INPUT="output/filtlong/output.fastq.gz"
OUTPUT="output/flye" # Note: this is a directory, not a file.

flye \
    --nano-raw $INPUT \
    --meta \
    --threads $SLURM_NPROCS \
    --out-dir $OUTPUT \
    --iterations 2
    
```


### Polishing with Racon and Medaka

Create a file named 03_polish-racon.sh and submit the job with sbatch:


```
#!/bin/bash
#SBATCH --job-name=polish1-racon
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G
#SBATCH --partition=normal



>&2 echo "Round 1"
>&2 echo "Mapping minimap2 ..."
minimap2 \
    -x map-ont \
    -t 8 \
    output/flye/assembly.fasta \
    output/filtlong/output.fastq.gz \
    > output/racon_art/minimap2_round1.paf

>&2 echo "Correcting Racon ..."
racon             -t 8             output/filtlong/output.fastq.gz             output/racon_art/minimap2_round1.paf             output/flye/assembly.fasta > output/racon_art/racon_round1.fna


>&2 echo "Round 2"
>&2 echo "Mapping minimap2 ..."
minimap2             -x map-ont             -t 8             output/racon_art/racon_round1.fna             output/filtlong/output.fastq.gz > output/racon_art/minimap2_round2.paf

>&2 echo "Correcting Racon ..."
racon             -t 8             output/filtlong/output.fastq.gz             output/racon_art/minimap2_round2.paf             output/racon_art/racon_round1.fna > output/racon_art/racon_round2.fna
```














