# Prokaryota dry lab: 03 Binning
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Friday 14th of April 2023


## Binning with Metabat2 🦇🗑️


So far we created an assembly containing all contigs (continuous sequences) from each of the organisms in the microbial community that we sequenced from the cattle rumen.

Presently, the assembly consists of thousands of contigs, each coming from a single species. By grouping together the contigs from each species present in our sample, we can create what is referred to as a MAG, short for Metagenome Assembled Genome.

Popular binning algorithms like the one used in Metabat2 utilize contig depth as a tell tale to link the contigs together that come from the same species. This is done by mapping the reads onto the assembly and then it counts the read depth of each contig. The smart thing here is that contigs coming from the same species will have similar depth, because the abundance is linked. Another vital contig statistic that binners use is the GC-content. Each species has its own intrinsic GC-content, and by grouping contigs further on GC-content, we might get good representatives for distinct species in our sample.

So, here we will first calculate the depth of each contig.

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
source activate /mnt/courses/BIO326/PROK/checkm2_condaenv

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


--- 

When the binning has completed, we can check the total length and number of contigs in each bin with assembly-stats.




### Bin QC

Quality control is something we can do on all levels of our microbiological workflow. We have been continuously running assembly-stats to follow how each tool has reshaped our assembly. Now the time has come to look into some more qualitative statistics with CheckM2 (https://github.com/chklovski/CheckM2/)

📝 Create a file named 03c_qc-checkm2.sh with the following contents, and submit the job with sbatch:


```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=qc-checkm2
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 8
#SBATCH --mem=16G

# Activate the conda environment. Note that we're using a separate conda environment for this software.
source activate /mnt/courses/BIO326/PROK/checkm2

# Define paths
in_dir="results/metabat2/"
out_dir="results/checkm2"


checkm2 predict --threads $SLURM_CPUS_PER_TASK --input $in_dir --output-directory $out_dir --extension .fa --force

```

