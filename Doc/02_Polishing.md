# Prokaryota dry lab: 02 Polishing
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Wednesday 12th of April 2023


## A note about the problems that we had last time 

Last time we ran into some issues with filtlong because the locale settings on your user accounts are not set up correctly. This is of course not your fault, but we think that we solved the issue, and we just need you to run a few commands. Hopefully we will not get into more trouble.

Once logged into orion, please copy and paste these commands.
```bash
echo "export LC_ALL=C; unset LANGUAGE" >> ~/.bash_profile # Fixes filtlong
echo """ 
if [ -f ~/.bashrc ]; then
  . ~/.bashrc
fi
""" >> ~/.bash_profile
export PS1="\[\e[32m\]\u\[\e[m\]\[\e[33m\]@\[\e[m\]\[\e[34m\]\h\[\e[m\] \[\e[35m\]\W\[\e[m\] \[\e[33m\]\\$\[\e[m\] " >> ~/.bash_profile # Makes your terminal always show your current working directory
exit

```

After you run these commands on orion, it should automatically log you out, and you can log in again. When you log in again, you should see your terminal having new colors:



<img width="299" alt="Screenshot 2023-03-22 at 15 11 12" src="https://user-images.githubusercontent.com/5913696/226930948-b2531d29-8810-48a0-bdee-756d3ff31ae4.png">





## Recap

Last time we filtered our raw nanopore reads with Filtlong, and then we created an assembly using Flye. Ideally, you should have a directory in your $SCRATCH directory named prok, containing the results for both Filtlong and Flye. 

For the analysis today we need the assembly from flye. We will account on it to be in the directory that is specified in the polishing scripts, and before we start, you should check that it is in fact present and is not an empty file.

```bash
ls -lh $SCRATCH/prok/results/flye/assembly.fasta
#> -rw-rw-r-- 1 cako nobody 233M Mar 22 15:16 /mnt/SCRATCH/cako/prok/results/flye/assembly.fasta

```

If you do not see that your file either does not exist, or does not have a size comparable to the 233M (Megabytes) that is shown in the output above, you can copy it from the demo directory using the following command:

```bash
cp -v /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/* $SCRATCH/prok/results/flye/
#> cp: omitting directory '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/00-assembly'
#> cp: omitting directory '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/10-consensus'
#> cp: omitting directory '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/20-repeat'
#> cp: omitting directory '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/30-contigger'
#> cp: omitting directory '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/40-polishing'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly.fasta' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly.fasta'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_graph.gfa' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_graph.gfa'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_graph.gv' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_graph.gv'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_info.txt' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_info.txt'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/flye.log' -> '/mnt/SCRATCH/cako/prok/results/flye/flye.log'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/params.json' -> '/mnt/SCRATCH/cako/prok/results/flye/params.json'
```

Don't worry about the "cp: omitting directory" warnings for now, these subdirectories are not necessary for polishing.





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
#SBATCH --output slurm-%j-%x.out.log


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


So far we created a draft assembly, containing all contigs (continuous sequences) from each of the organisms in the microbial community that we sequenced (cattle rumen). In this *draft assembly* we then proceeded to remove the potential sequencing errors in order to create a polished *assembly*.

Now comes the task of extracting each of the several species present in our sample. This process in referred to as *binning*. This is because we effectively classify each contig in the assembly as coming from one of several sources (species). The way the binning algorithms typically work is that they look into the abundance of each of the contigs. This is done by mapping the reads onto the assembly (like we did in the racon polishing step, and then counting the read depth of each contig). The smart thing here is that contigs coming from the same species will have similar depth, because the abundance is linked. The binning algorithm uses this information to classify each contig ??


So, here we will first calculate the depth of each contig

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












