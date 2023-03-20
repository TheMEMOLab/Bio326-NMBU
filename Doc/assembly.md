# Prokaryota dry lab: 01 Metagenomic Assembly
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Date??

Hello! Welcome to the metagenomic dry lab session! Here we will take the raw basecalled reads from the nanopore sequenator and try to reconstruct the microbial genomes that they come from.

#### Checking access to conda environment
We created a conda environment in /mnt/courses/BIO326/PROK/condaenv that contains all the software that we'll use to complete the metagenomic assembly.

You can activate it temporarily and check that all software is ready to rock and roll:

```bash
conda activate /mnt/courses/BIO326/PROK/condaenv
flye --version
#> 2.9.1-b1780
assembly-stats -v
#> Version: 1.0.1
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

We already merged the files from the three sequencing runs and performed rudimentary quality control. Everything looks great! Your reads reside in "/mnt/courses/BIO326/PROK/data/metagenomic_assembly/".

There is one file named "nanopore_reads_filtered.fastq.gz"

**Task**: Copy this file into your preferred working directory, so we can start using it in the next step.


## Assemble reads into a draft assembly 🏗

The sequenced reads represents fragments of biological genomic sequences - These we want to reconstruct inside the computer. It works a bit by solving a puzzle: Each piece is compared to many other pieces, until the whole picture can be put together. In the same way, we overlap the reads with one another, until we get long continuous sequences a.k.a. "contigs". These contigs can then be put together into "scaffolds" that represent larger parts of the original biological chromosomes and plasmids that harbor the prokaryotic cells.

Here we will use the Flye assembler (https://github.com/fenderglass/Flye/). It takes in reads from from genomic sequencing, and puts out a long draft assembly that contains sequence scaffolds from all the species that are present in the samples we sequenced.

📝 Create a file named 01_assemble-flye.sh with the following contents, and submit the job with sbatch: Remember to change the "in" variable to point to your copy of the nanopore_reads.fastq.gz file.

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=assemble-flye
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task 8
#SBATCH --mem=30G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in="<path to nanopore reads>"
out="results/flye" # Note: this is a directory, not a file.

mkdir results


flye --meta --nano-hq $in --threads $SLURM_CPUS_PER_TASK --out-dir $out --iterations 2


```


**Bonus points**: If you look closely at the flye program call in the sbatch script above, you can see that we're passing the "--meta" argument to flye. By investigating the Flye documentation (https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md), can you explain briefly what the "meta" mode does, and argue why we want to use it in this setting?

---

When Flye finishes, you will see a lot of output files in the results/flye directory. The main result file that we'll continue with is the results/flye/assembly.fasta file.

```bash
ls -lh results/flye/
#> total 480M
#>    2 Mar 15 15:00 00-assembly/
#>    3 Mar 15 15:49 10-consensus/
#>    6 Mar 15 16:13 20-repeat/
#>    7 Mar 15 16:15 30-contigger/
#>   11 Mar 16 04:17 40-polishing/
#> 233M Mar 16 04:18 assembly.fasta
#> 224M Mar 16 04:18 assembly_graph.gfa
#> 3.1M Mar 16 04:17 assembly_graph.gv
#> 384K Mar 16 04:18 assembly_info.txt
#>  20M Mar 16 04:18 flye.log
#>   92 Mar 16 04:17 params.json
```


You can investigate some basic statistics of the Flye assembly using the assembly-stats software:

```bash
assembly-stats results/flye/assembly.fasta 
#> stats for output/flye/assembly.fasta
#> sum = 239436989, n = 10266, ave = 23323.30, largest = 1194178
#> N50 = 33581, n = 1668
#> N60 = 25598, n = 2487
#> N70 = 19593, n = 3565
#> N80 = 15137, n = 4951
#> N90 = 10904, n = 6804
#> N100 = 114, n = 10266
#> N_count = 0
#> Gaps = 0
```

As the algorithms used in Flye are not deterministic, your output may vary slightly.