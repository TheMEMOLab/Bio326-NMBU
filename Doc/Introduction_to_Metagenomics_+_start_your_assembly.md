# Prokaryota dry lab: 01 Metagenomic Reconstruction
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples

Hello! Welcome to the metagenomic dry lab session! Here we will take the raw basecalled reads from the nanopore sequenator and try to reconstruct the microbial genomes that they come from.

#### Checking access to conda environment
There is a conda environment in /mnt/courses/BIO326/PROK/condaenv that we will use in each of the slurm batch scripts.

You can activate it temporarily and check that all software is ready to rock and roll:

```bash
conda activate /mnt/courses/BIO326/PROK/condaenv
flye --version
#> 2.9.1-b1780
filtlong --version
#> Filtlong v0.2.1
flye --version
#> 2.9.1-b1780
minimap2 --version
#> 2.24-r1122
racon --version
#> 1.5.0
medaka --version
#> medaka 1.7.2
metabat2 2> >(grep -oE "\(version.+$")
#> (version 2:2.15 (Bioconda); 2022-11-19T22:34:33)
conda deactivate

```



## Quality control and filtering of the raw reads 🛂

Your raw reads from the prokaryotic sequencing session reside in "/mnt/courses/BIO326/PROK/data/metagenomic_assembly/".

There is one file named "raw_reads_nanopore.fastq.gz"

Fastq is a raw read format containing a base quality scores for nucleotide position along each read.
Sometimes when we sequence, we see a lot of low quality reads that we want to get rid of, because they mostly contain noise that confuse the downstream analysis.

By specifying `--min_length 1000` and `--keep_percent 90` we keep only the reads that live up to these requirements.


📝 Create a file named 01_filter-filtlong.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=filter-filtlong
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem=8G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define IO
in="/mnt/courses/BIO326/PROK/data/metagenomic_assembly/raw_reads_nanopore.fastq.gz"
out="output/filtlong/output.fastq.gz"

# Make sure that the output directory exists
mkdir --parents $(dirname $out)


filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    $in \
    | gzip \
    > $out


```

Now check your output/filtlong/ directory. There should be a compressed fastq output file.

```bash
tree output/filtlong/
#> output/filtlong/
#> └── output.fastq.gz
#> 
#> 0 directories, 1 file
```

We can visualize the quality of these reads with nanoplot.

📝 Create a file named 01b_visualize-nanoplot.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=visualize-nanoplot
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task 2
#SBATCH --mem=8G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define IO
in="output/filtlong/output.fastq.gz"


NanoPlot \
    --threads $SLURM_NPROCS \
    --fastq $in \
    --plots hex dot 


```

Now check your output/filtlong/ directory. There should be a compressed fastq output file.

```bash
tree output/filtlong/
#> output/filtlong/
#> └── output.fastq.gz
#> 
#> 0 directories, 1 file
```



## Assemble reads into a draft assembly 🏗

Then we will use flye in "meta" mode. Flye builds contiguous sequences by overlapping each read.
You can read more about how to configure flye here: https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md

📝 Create a file named 02_assemble-flye.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=assemble-flye
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define IO
in="output/filtlong/output.fastq.gz"
out="output/flye" # Note: this is a directory, not a file.

# Make sure that the output directory exists
#mkdir --parents $out


flye \
    --nano-raw $in \
    --meta \
    --threads $SLURM_NPROCS \
    --out-dir $out \
    --iterations 2


```


## Polishing with Racon and Medaka ✨

### Racon

📝 Create a file named 03_polish-racon.sh with the following contents, and submit the job with sbatch:


```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=polish1-racon
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define IO
in_draft_assembly="output/flye/assembly.fasta"
in_reads="output/filtlong/output.fastq.gz"
out_polished_assembly="output/racon/racon_round2.fna"

# Make sure that the output directory exists
mkdir --parents $(dirname $out_polished_assembly)


>&2 echo "Round 1"
>&2 echo "Mapping minimap2 ..."
minimap2 \
    -x map-ont \
    -t 8 \
    $in_draft_assembly \
    $in_reads \
    > output/racon/minimap2_round1.paf

>&2 echo "Correcting Racon ..."
racon \
    -t 8 \
    $in_reads \
    output/racon/minimap2_round1.paf \
    $in_draft_assembly > output/racon/racon_round1.fna


>&2 echo "Round 2"
>&2 echo "Mapping minimap2 ..."
minimap2 \
    -x map-ont \
    -t 8 \
    output/racon/racon_round1.fna \
    $in_reads > output/racon/minimap2_round2.paf

>&2 echo "Correcting Racon ..."
racon \
    -t 8 \
    $in_reads \
    output/racon/minimap2_round2.paf \
    output/racon/racon_round1.fna > $out_polished_assembly


```



### Medaka


📝 Create a file named 04_polish-medaka.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=polish1-racon
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define IO
in_assembly="output/racon_art/racon_round2.fna"
in_reads="output/filtlong/output.fastq.gz"
OUTPUT="output/medaka_art"

# Make sure that the output directory exists
#mkdir --parents $OUTPUT


medaka_consensus \
    -t $SLURM_NPROCS \
    -d $in_assembly \
    -i $in_reads \
    -o $OUTPUT \
    -m r1041_e82_260bps_hac_g632


```

## Binning with Metabat2 🗑️🗑️




📝 Create a file named 05_bin-metabat.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=bin-metabat
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=8G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define IO
in_assembly="output/medaka_art/consensus.fasta"
out_bins="output/metabat2/bin"

# Make sure that the output directory exists
mkdir --parents $(dirname $out_bins)


metabat2 \
    --numThreads 4 \
    --inFile $in_assembly \
    --outFile $out_bins \
    --minClsSize 1000000


```












