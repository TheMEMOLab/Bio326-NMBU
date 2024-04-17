# Prokaryota dry lab: 03 Binning
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Wednesday 17th of April 2024



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

<table><tr><td>

#### (In case of missing data from last time)

If you get a "No such file or directory" error, you can copy our premade files into your directory and continue with the binning exercise below. <ins>Be aware</ins> that this action will possibly overwrite your own files if they do indeed exist.


```bash
# Only run this if you didn't previously successfully run filtlong and medaka
mkdir -p $SCRATCH/prok/results/filtlong/
cp /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/filtlong/output.fastq.gz $SCRATCH/prok/results/filtlong/
mkdir -p $SCRATCH/prok/results/medaka/
cp /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/medaka/consensus.fasta $SCRATCH/prok/results/medaka/
cd $SCRATCH/prok
```

</td></tr></table>




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
#SBATCH --mem=24G
#SBATCH --output slurm-%j-%x.out.log
#SBATCH -p smallmem,hugemem,hugemem-avx2

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

When the job finishes, you can see the contig depths in the file below:
```bash
cat $SCRATCH/prok/results/contig_depths/depth.tsv | column -t | less -S
#> contigName    contigLen    totalAvgDepth  bam_for_depths.bam  bam_for_depths.bam-var
#> contig_6249   42425        4.15245        4.15245             16.7985
#> contig_746    10384        0              0                   0
#> contig_850    13519        0              0                   0
#> contig_4906   9255         0.37441        0.37441             0.234247
#> contig_2344   21841        5.9964         5.9964              19.6254
#> contig_6063   40379        6.67968        6.67968             15.6048
#> contig_2018   15964        2.46883        2.46883             3.64768
#> contig_10853  15213        0.907256       0.907256            0.633342
#> contig_5190   17779        1.36417        1.36417             1.6962
#> contig_9671   10145        0.109855       0.109855            0.0978116
#> contig_1104   105858       3.80065        3.80065             6.97997
#> contig_7858   31531        6.38456        6.38456             40.9298
#> contig_9245   16570        1.68544        1.68544             2.59145
#> contig_2602   21434        1.71077        1.71077             2.06949
#> contig_7752   14634        0.405482       0.405482            0.241066
#> ...
#> (press q to exit)
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
#SBATCH --output slurm-%j-%x.out.log
#SBATCH -p smallmem,hugemem,hugemem-avx2


# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths

bin_prefix="bin"
in_assembly="$SCRATCH/prok/results/medaka/consensus.fasta"
in_depth="$SCRATCH/prok/results/contig_depths/depth.tsv"
out_bins="$SCRATCH/prok/results/metabat2/${bin_prefix}"

# Make sure that the output directory exists
mkdir --parents $(dirname $out_bins)


metabat2  --numThreads $SLURM_CPUS_PER_TASK  --inFile $in_assembly  --outFile $out_bins  --abdFile $in_depth  --minClsSize 1000000



```

Binning should take a few minutes.

--- 



The binner automatically decides how many bins it thinks are present. If you run ls on its directory, you should see some 30 bins.

```bash
ls $SCRATCH/prok/results/metabat2/
#> bin.1.fa   bin.12.fa  bin.15.fa  bin.18.fa  bin.20.fa  bin.23.fa  bin.26.fa  bin.29.fa  bin.31.fa  bin.34.fa  bin.5.fa  bin.8.fa
#> bin.10.fa  bin.13.fa  bin.16.fa  bin.19.fa  bin.21.fa  bin.24.fa  bin.27.fa  bin.3.fa   bin.32.fa  bin.35.fa  bin.6.fa  bin.9.fa
#> bin.11.fa  bin.14.fa  bin.17.fa  bin.2.fa   bin.22.fa  bin.25.fa  bin.28.fa  bin.30.fa  bin.33.fa  bin.4.fa   bin.7.fa
```

Depending on destiny, you might see more or less than 30 files. Don't worry, this is normal. Metabat2 uses a non-deterministic algorithm, and if you rerun it, you might get a different number of bins.


We can check the total length and number of contigs in each bin with assembly-stats.
```bash 
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats -t $SCRATCH/prok/results/metabat2/*.fa | column -t | less -S
#> filename                    total_length  number  mean_length  longest  shortest  N_count  Gaps  N50      N50n  N70      N70n  N90      N90n
#> results/metabat2/bin.1.fa   1368180       34      40240.59     212034   5187      0        0     61619    7     35298    13    19171    23
#> results/metabat2/bin.10.fa  2335282       31      75331.68     304038   3069      0        0     125057   6     90932    10    45038    18
#> results/metabat2/bin.11.fa  7434991       306     24297.36     80231    5294      0        0     27951    86    20575    148   14098    235
#> results/metabat2/bin.12.fa  2814547       123     22882.50     90800    5402      0        0     29328    30    19676    54    11538    92
#> results/metabat2/bin.13.fa  1298639       2       649319.50    1187215  111424    0        0     1187215  1     1187215  1     1187215  1
#> results/metabat2/bin.14.fa  1150712       52      22129.08     115632   3682      0        0     28276    13    17582    23    12284    39
#> results/metabat2/bin.15.fa  2054121       15      136941.40    434392   21264     0        0     173722   4     149421   7     84010    11
#> results/metabat2/bin.16.fa  7176178       273     26286.37     148334   3715      0        0     36670    64    22000    114   13444    199
#> results/metabat2/bin.17.fa  2806080       23      122003.48    410382   7343      0        0     223668   5     113842   9     54987    15
#> ... 
#> (press q to exit)
```

---

## Qualitative Bin Analysis and Annotation

Now, our bins are created. These bins represent individual bacterial and archaeal genomes of the the microbiome in the rumen material that we have sequenced. 

We will use a genomes-to-report pipeline named Assemblycomparator2 to characterize the functional potential of our bins. You can read more about this pipeline [here]([https://github.com/cmkobel/assemblycomparator2](https://github.com/cmkobel/assemblycomparator2?tab=readme-ov-file#what-analyses-does-it-do)).

<img width="150" alt="snakemake logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/7188e748-9d37-43ae-a5d5-100e9560df1f">

Assemblycomparator2 is based on [Snakemake](https://snakemake.github.io/). 

![image](https://github.com/TheMEMOLab/Bio326-NMBU/assets/5913696/bee98c41-6d96-4985-8bea-5784cb008d8c)
_Figure showing the jobs that assemblycomparator2 computes on the bins, and in which order. This is a directed acyclic graph._


Assemblycomparator2 is already installed on Orion, and can be called after activating the following conda environment.

```bash
source activate /mnt/courses/BIO326/PROK/conda_ac26
#> (/mnt/courses/BIO326/PROK/conda_ac26) [bio326-2024-n@login: ~]$ 

# Check that you have the newest version.
asscom2 --version
#> assemblycomparator2 v2.6.1

```

The pipeline is then ready to rock and roll. Assemblycomparator2 is configured to use the Slurm workload manager behind-the-scenes, so there is no need to create or launch any more shell (.sh) scripts. 

To get going analyzing our freshly made genomic bins, first enter the directory where your bins reside like so.

```bash 

cd $SCRATCH/prok/results/metabat2/
ls *.fa
#> bin.1.fa 
#> bin.2.fa
#> ...
#> bin.N.fa
```

Then, launch the pipeline with this command: (We will call it with an ampersand (&) at the end, to "fork" the process and let it continue running even if you disconnect your laptop from the network.)

```bash 

asscom2 --unlock # Clears previously started, failed runs
asscom2 --until report assembly_stats sequence_lengths prokka busco checkm2 mashtree dbcan & 

```

<table><tr><td>
  
#### (How to stop a "forked" process)
If you regret starting the pipeline, you can -in the same terminal window- press "fg" on your keyboard and then hit enter, followed by ctrl+c once to stop the process. "fg" brings the forked process to the foreground, and ctrl+c interrupts the pipeline

</td></tr></table>

The "--until" argument lets the pipeline know to only run the specified analyses. In this case we're running the ones that are relevant for comparing bins or MAGs.

You will se a lot of output in your terminal. This is because the pipeline runs all independent analysis jobs at the same time. Some jobs run for each bin, and others run a comparison across all bins in a single job. It will likely take several hours for the complete pipeline to finish, but some of the faster jobs (like sequence_lengths and busco) might finish after just 20 minutes. 

It is a good idea to open a new tab in your terminal window, and log in with another session. Then you can let assemblycomparator run in the first login window while you browse the results as they finish -real time- in the second.

Log into Orion in a second tab and surveil the output from assemblycomparator by running "tree" on the newly created "results_ac2" directory where assemblycomparator outputs its results.

```bash
tree -L 3 $SCRATCH/prok/results/metabat2/results_ac2
#> ???
```

The "-L <N>" argument lets tree know to stop listing files after hitting a depth level of N in the directory. This is to avoid overflowing your terminal window.

---

When all of the jobs in the pipeline have finished, an report document will reside in "results_ac2/report_metabat2.html". You should download this file to your computer and open it with a web browser. On all platforms (windows, mac, linux), it should be possible to use  the FileZilla client to download this .html file (https://filezilla-project.org/). You can also use `rsync` on the command line (see here: https://orion.nmbu.no/en/Copydata).

If you're having trouble with the pipeline report document - Maybe it won't download or isn't created in the first place - Ask for help. Or download our demonstration report from here: [report_metabat2_orion_carl.html.zip](https://github.com/TheMEMOLab/Bio326-NMBU/files/11211669/report_metabat2_orion_carl.html.zip)

