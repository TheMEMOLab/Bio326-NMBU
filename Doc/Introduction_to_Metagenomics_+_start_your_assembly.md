# Prokaryota dry lab: 01 Metagenomic Reconstruction
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples

Hello! Welcome to the metagenomic dry lab session! Here we will take the raw basecalled reads from the nanopore sequenator and try to reconstruct the microbial genomes that they come from.

#### Checking access to conda environment
There is a conda environment in /mnt/courses/BIO326/PROK/condaenv that we will use in each of the slurm batch scripts.

You can activate it temporarily and check that all software is ready to rock and roll:

```bash
conda activate /mnt/courses/BIO326/PROK/condaenv
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




## Polishing with Racon and Medaka ✨

The sequenced reads have an error rate around 1/100. This is not a big deal, as the assembly algorithm in flye is tolerant and perfectly happy with this. The problem for us though is that some of the positions in the draft assembly might then be representing these errors, rather than the actual biological sequence of the organism. Luckily there is a process referred to as *polishing* that can reduce the number of erroneous positions by overlapping the reads and calculating the probabilites for error. This process is reminiscent of extracting a *consensus* sequence from a set of aligned genomes.

Here we are going to apply several rounds of racon, and a final round of medaka. Note that Flye also has an internal polisher that we have already applied. As it turns out, Racon and Medaka have better performance, so we'll apply those as well.

### Racon 🦝

In order to run Racon we need to give it our draft assembly and our reads. First we'll map the reads to the draft with minimap2, then racon sifts through the alignment in the .paf file, and outputs a corrected assembly. This process is repeated for a total of two rounds of polishing.

If you want to know more about how to set up Racon, you can read about it here: https://github.com/isovic/racon#usage

📝 Create a file named 02a_polish1-racon.sh with the following contents, and submit the job with sbatch:


```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=polish1-racon
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_draft_assembly="results/flye/assembly.fasta"
in_reads="results/filtlong/output.fastq.gz"
out_polished_assembly="results/racon/racon_round2.fna"

# Make sure that the output directory exists
mkdir results/racon/



# Mapping minimap2 round 1
minimap2 -x map-ont -t $SLURM_CPUS_PER_TASK $in_draft_assembly $in_reads > results/racon/minimap2_round1.paf

# Correcting Racon round 1
racon -t $SLURM_CPUS_PER_TASK $in_reads results/racon/minimap2_round1.paf $in_draft_assembly > results/racon/racon_round1.fna



# Mapping minimap2 round 2
minimap2 -x map-ont -t $SLURM_CPUS_PER_TASK results/racon/racon_round1.fna $in_reads > results/racon/minimap2_round2.paf

# Correcting Racon round 2
racon -t $SLURM_CPUS_PER_TASK $in_reads results/racon/minimap2_round2.paf results/racon/racon_round1.fna > $out_polished_assembly


```

### Medaka 🐟

From the polished output of two rounds of Racon, we have an assembly that is quite good, but we can make it even better. Here we'll apply one round of medaka polishing. Medaka works much the same way but uses a different internal algorithm. In this case, we're telling medaka which exact sequencing platform we used for sequencing the reads - This is to let Medaka know about the specifics of the errors that this specific platform creates??.

If curious, you can read more about how to set up Medaka here: https://github.com/nanoporetech/medaka#usage


📝 Create a file named 02b_polish2-medaka.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=polish2-medaka
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=4G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="results/racon/racon_round2.fna"
in_reads="results/filtlong/output.fastq.gz"
out="results/medaka"


medaka_consensus -t $SLURM_CPUS_PER_TASK -d $in_assembly -i $in_reads -o $out -m r1041_e82_400bps_sup_g615


```




## Binning with Metabat2 🦇🗑️


### Calculating contig depths

📝 Create a file named 03a_depth-minimap.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=depth-minimap
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=8G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="results/medaka/consensus.fasta"
in_reads="results/filtlong/output.fastq.gz"
out_alignment="results/contig_depths/bam_for_depths.bam"
out_depth="results/contig_depths/depth.tsv"

# Make sure that the output directory exists
mkdir results/contig_depths/


# Map reads to polished assembly and sort the alignment
minimap2 -ax map-ont --sam-hit-only -t $SLURM_CPUS_PER_TASK $in_assembly $in_reads | samtools sort -@ $SLURM_CPUS_PER_TASK -o $out_alignment

# Calculate depths of above alignment
jgi_summarize_bam_contig_depths --outputDepth $out_depth $out_alignment


```















### Binning 

📝 Create a file named 03b_bin-metabat.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=bin-metabat
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=8G

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="results/medaka/consensus.fasta"
in_depth="results/contig_depths/depth.tsv"
out_bins="results/metabat2/bin"

# Make sure that the output directory exists
mkdir --parents $(dirname $out_bins)


metabat2  --numThreads $SLURM_CPUS_PER_TASK  --inFile $in_assembly  --outFile $out_bins  --abdFile $in_depth  --minClsSize 1000000


```












