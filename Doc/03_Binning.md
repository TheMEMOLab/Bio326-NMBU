# Prokaryota dry lab: 03 Binning
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Friday 14th of April 2023



## Checkpoint

Today we will need the final polished assembly from medaka as well as the original reads that we produced with filtlong. The binning scripts that we will run today expect these input files to be present in the $SCRATCH/prok subdirectories. 

You can check that your files are present and non-empty with these commands:

```bash 
cd $SCRATCH/prok
du -h results/filtlong/output.fastq.gz 
du -h results/medaka/consensus.fasta
#> 4.1G    results/filtlong/output.fastq.gz
#> 224M    results/medaka/consensus.fasta
```

If you get a "No such file or directory" error, you can copy our premade files into your directory and continue with the binning exercise below. <ins>Be aware</ins> that this action will possibly overwrite your own files if they do indeed exist.


```bash
# Only run this if you didn't previously successfully run filtlong and medaka
mkdir -p $SCRATCH/prok/results/filtlong/
cp /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/filtlong/output.fastq.gz $SCRATCH/prok/results/filtlong/
mkdir -p $SCRATCH/prok/results/medaka/
cp /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/medaka/consensus.fasta $SCRATCH/prok/results/medaka/
cd $SCRATCH/prok
```




## Binning with Metabat2 🦇🗑️


So far we created an assembly containing all contigs (continuous sequences) from each of the organisms in the microbial community that we sequenced from the cattle rumen.

Presently, the assembly consists of thousands of contigs, each coming from a single species. By grouping together the contigs from each species present in our sample, we can create what is referred to as a MAG, short for Metagenome Assembled Genome.

Popular binning algorithms like the ones used in Metabat2 utilize contig depth as a tell tale to link the individual contigs together that come from the same species. This is done by mapping the original reads onto the assembly and then counting the read depth of each contig. The smart thing here is that contigs coming from the same species will have similar depth. Another vital contig statistic that binners use is the GC-content. Each species has its own intrinsic GC-content, and by grouping contigs further on GC-content -in this case by counting the tetranucleotide frequency- we might get good representatives for distinct species in our sample. If our bins live up to our requirements, we can refer to them as MAGs.

So, here we will first calculate the depth of each contig.

### Calculating contig depths

📝 Create a file named 03a_depth-minimap.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=depth-minimap
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 8
#SBATCH --mem=12G
#SBATCH --output slurm-%j-%x.out.log
#SBATCH -p smallmem,hugemem

# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="$SCRATCH/prok/results/medaka/consensus.fasta"
in_reads="$SCRATCH/prok/results/filtlong/output.fastq.gz"
out_alignment="$SCRATCH/prok/results/contig_depths/bam_for_depths.bam"
out_depth="$SCRATCH/prok/results/contig_depths/depth.tsv"

# Make sure that the output directory exists
mkdir --parents $SCRATCH/prok/results/contig_depths/


# Map reads to polished assembly and sort the alignment
minimap2 -ax map-ont --sam-hit-only -t $SLURM_CPUS_PER_TASK $in_assembly $in_reads | samtools sort -@ $SLURM_CPUS_PER_TASK -o $out_alignment

# Calculate depths of above alignment
jgi_summarize_bam_contig_depths --outputDepth $out_depth $out_alignment


```

Calculating the depth of each contig should take around 30 minutes.


### Binning 

📝 Create a file named 03b_bin-metabat.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=bin-metabat
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem=8G
#SBATCH --output slurm-%j-%x.out.log
#SBATCH -p smallmem,hugemem


# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="$SCRATCH/prok/results/medaka/consensus.fasta"
in_depth="$SCRATCH/prok/results/contig_depths/depth.tsv"
out_bins="$SCRATCH/prok/results/metabat2/bin"

# Make sure that the output directory exists
mkdir --parents $(dirname $out_bins)


metabat2  --numThreads $SLURM_CPUS_PER_TASK  --inFile $in_assembly  --outFile $out_bins  --abdFile $in_depth  --minClsSize 1000000


```

Binning should take a few minutes.

--- 

When the binning has completed, we can check the total length and number of contigs in each bin with assembly-stats.

The binner automatically decides how many bins it thinks are present. If you run ls on its directory, you should see some 30 bins.

```bash
ls results/metabat2/
#> bin.1.fa   bin.12.fa  bin.15.fa  bin.18.fa  bin.20.fa  bin.23.fa  bin.26.fa  bin.29.fa  bin.31.fa  bin.34.fa  bin.5.fa  bin.8.fa
#> bin.10.fa  bin.13.fa  bin.16.fa  bin.19.fa  bin.21.fa  bin.24.fa  bin.27.fa  bin.3.fa   bin.32.fa  bin.35.fa  bin.6.fa  bin.9.fa
#> bin.11.fa  bin.14.fa  bin.17.fa  bin.2.fa   bin.22.fa  bin.25.fa  bin.28.fa  bin.30.fa  bin.33.fa  bin.4.fa   bin.7.fa
```




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
in_dir="$SCRATCH/prok/results/metabat2/"
out_dir="$SCRATCH/prok/results/checkm2"


checkm2 predict --threads $SLURM_CPUS_PER_TASK --input $in_dir --output-directory $out_dir --extension .fa --force

```




## Qualitative Bin Analysis and Annotation

We will use a genomes-to-report pipeline named assemblycomparator2. You can read more about it [here](https://github.com/cmkobel/assemblycomparator2).

You can install a shortcut to run the pipeline, by calling this instruction in your terminal:

```bash

echo """

export SNAKEMAKE_CONDA_PREFIX=/mnt/users/cako/prok23/assemblycomparator2/conda_base
export ASSCOM2_BASE=/mnt/courses/BIO326/PROK/assemblycomparator2
alias assemblycomparator2='/mnt/orion/opt/conda/miniconda3/bin/conda run \
    --live-stream \
    --prefix /net/fs-2/scale/OrionStore/Courses/BIO326/PROK/assemblycomparator2/ac2 \
    snakemake \
        --snakefile /mnt/courses/BIO326/PROK/assemblycomparator2/snakefile \
        --profile /mnt/courses/BIO326/PROK/assemblycomparator2/profiles/slurm-nmbu-orion/ \
        --configfile /mnt/courses/BIO326/PROK/assemblycomparator2/config.yaml'

""" >> ~/.bashrc && source ~/.bashrc

```


Assemblycomparator2 is then installed and set up specifically to use the slurm/sbatch system on Orion, so there is no need to create or launch any shell (.sh) scripts. 



First, enter the directory where your bins reside.

```bash 

cd $SCRATCH/prok/results/metabat2/
ls *.fa
#> rumen1.fa 
#> rumen2.fa
#> ...
#> rumenN.fa
```

Then, launch assemblycomparator2 with this command: (We will start it with an ampersand (&) at the end, to fork the process and let it continue running even if you disconnect your laptop from the network.)

```bash 

assemblycomparator2 --until assembly_stats sequence_lengths prokka busco checkm2 kraken2 gtdbtk & 

```

The "--until" argument lets the pipeline know to only run the specified analyses. In this case we're running the ones that are relevant for comparing bins or MAGs.

You will se a lot of output in your terminal. This is because assemblycomparator2 runs all independent analysis jobs at the same time. Some jobs run for each bin, and others run a comparison across all bins in a single job. It will likely take several hours for the complete pipeline to finish, but some of the faster jobs (like sequence_lengths, busco and kraken2) might finish after just 20 minutes. 

It is a good idea to open a new tab in your terminal window, and log in with another session. Then you can let assemblycomparator run in the first login window while you browse the results as they finish -real time- in the second.

Log into Orion in a second tab and surveil the output from assemblycomparator by running "tree" on the newly created "results_ac2" directory where assemblycomparator outputs its results.

```bash
tree -L 3
#> ???
```

The "-L 3" argument lets tree know to stop listing files after hitting a depth level of 2 in the directory. This is to avoid overflowing your terminal window.

---

When all of the jobs in the pipeline have finished, an report document will reside in "results_ac2/report_metabat2.html". You should download this file to your computer and open it with a web browser.

All platforms (windows, mac, linux) should be able to work with the FileZilla client (https://filezilla-project.org/).

If you're having trouble with the pipeline report document - Maybe it won't download or isn't created in the first place - You can download our report from here: [report_metabat2_orion_carl.html.zip](https://github.com/TheMEMOLab/Bio326-NMBU/files/11211669/report_metabat2_orion_carl.html.zip)

