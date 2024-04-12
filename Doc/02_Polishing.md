# Prokaryota dry lab: 02 Polishing
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Friday 12th of April 2024



## Recap 

Last time we filtered our raw nanopore reads with Filtlong, and then we created an assembly using Flye. Ideally, you should have a directory in your personal $SCRATCH directory named prok, containing the results for both Filtlong and Flye. 

For the analysis today we need the assembly from Flye. We will account on it to be in the directory that is specified in the polishing scripts, and before we start, you should check that it is in fact present and is not an empty file.

First let's go into your $SCRATCH/prok/ directory and list its contents.

```bash
mkdir -p $SCRATCH/prok/ && cd $SCRATCH/prok/ && pwd
#> /mnt/SCRATCH/cako/prok
```

Then we list the contents of this directory. Don't expect your files to look exactly like these.

```bash
ls
#> 01a_filter-filtlong.sh  raw_reads_nanopore.fastq.gz  slurm-11939921.out
#> 01b_assemble-flye.sh    results                      slurm-11939933.out
```


For the polishing exercise today we require the assembly from Flye to be present in a predefined directory. You can make sure that the output from Flye is present by running the `ls` command:


```bash
ls -lh $SCRATCH/prok/results/flye/
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

<table><tr><td>


If you see that your file either does not exist or does not have a size comparable to around 233M (megabytes) which is exemplified in the ls output above, you can copy the assembly file we made for you in the demo directory using the following command.


```bash
mkdir -p $SCRATCH/prok/results/flye/
cp -v /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly* $SCRATCH/prok/results/flye/
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly.fasta' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly.fasta'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_graph.gfa' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_graph.gfa'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_graph.gv' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_graph.gv'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_info.txt' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_info.txt'
```

We also need the reads from filtlong, these you can get from here:

```bash
mkdir -p $SCRATCH/prok/results/filtlong/
cp /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/filtlong/output.fastq.gz $SCRATCH/prok/results/filtlong/
```

</td></tr></table>

Now we're ready to continue with polishing.

---

## Polishing with Racon and Medaka ✨

A basic model of how polishing works is that the polisher stacks all relevant reads on top of the genome and decides for each position whether the present nucleotide letter is the best representative for that position, or not. There are several sources of variation that make draft assemblies _polishable_. The main sources are multistrain variation from closely related species as well as incorporation of sequencing errors during the sequencing process. Ideally, assemblers would be perfect, and we wouldn't have to perform polishing. But because of some noise or artefacts that are present in our data, we might make our genomes more truthful to their biological origin by performing these polishing steps.

Genome polishing is reminiscent of generating a consensus genome. Consensus genome creation is a term used in reference mapping. This is why you may incidentally see the term _consensus_ being used in the tools that we're gonna run.

Here we are going to apply several rounds of racon, and a final round of medaka. Note that Flye also has an internal polisher that we have already applied (if you look closely at the program call for Flye, it says "--iterations 2" ...). As it turns out, Racon and Medaka have better performance, so we'll apply those as well.

Deciding in which order and how many iterations to apply each polishing tool typically is decided with trial and error (empirical).


### Racon 🦝

In order to run Racon we need to give it our draft assembly and our reads. First we'll map the reads to the draft with minimap2, then racon sifts through the alignment in the .paf file, and outputs a corrected assembly. This process is repeated for a total of two rounds of polishing.

If you want to know more about how to set up Racon, you can read about it here: https://github.com/isovic/racon#usage

📝 Create a file named 02a_polish1-racon.sh with the following contents, and submit the job with sbatch:


```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=polish1-racon
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task 8
#SBATCH --mem=16G
#SBATCH --output slurm-%j-%x.out.log
#SBATCH -p smallmem,hugemem,hugemem-avx2


# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_draft_assembly="$SCRATCH/prok/results/flye/assembly.fasta"
in_reads="$SCRATCH/prok/results/filtlong/output.fastq.gz"
out_polished_assembly="$SCRATCH/prok/results/racon/racon_round2.fna"

# Make sure that the output directory exists
mkdir $SCRATCH/prok/results/racon/


# Mapping minimap2 round 1
minimap2 -x map-ont -t $SLURM_CPUS_PER_TASK $in_draft_assembly $in_reads > $SCRATCH/prok/results/racon/minimap2_round1.paf

# Correcting Racon round 1
racon -t $SLURM_CPUS_PER_TASK $in_reads $SCRATCH/prok/results/racon/minimap2_round1.paf $in_draft_assembly > $SCRATCH/prok/results/racon/racon_round1.fna


# Mapping minimap2 round 2
minimap2 -x map-ont -t $SLURM_CPUS_PER_TASK $SCRATCH/prok/results/racon/racon_round1.fna $in_reads > $SCRATCH/prok/results/racon/minimap2_round2.paf

# Correcting Racon round 2
racon -t $SLURM_CPUS_PER_TASK $in_reads $SCRATCH/prok/results/racon/minimap2_round2.paf $SCRATCH/prok/results/racon/racon_round1.fna > $out_polished_assembly


```

**NB! You can copy the script from ```/mnt/courses/BIO326/PROK/scripts/BIO326_24/02a_polish1-racon.sh```**


<table><tr><td>

#### Assembly-stats progress check on Racon results

Let's first check the assembly-stats on the original Flye assembly.

```bash
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats $SCRATCH/prok/results/flye/assembly.fasta
#> stats for /mnt/SCRATCH/bio326-2023-19/prok/results/flye/assembly.fasta
#> sum = 239251443, n = 10266, ave = 23305.23, largest = 1190509
#> N50 = 33554, n = 1668
#> N60 = 25596, n = 2486
#> N70 = 19567, n = 3564
#> N80 = 15104, n = 4951
#> N90 = 10900, n = 6804
#> N100 = 114, n = 10266
#> N_count = 0
#> Gaps = 0

```

Now, wait until Racon finishes, also run assembly-stats on that one, and compare to Flye.

(If you want to see the Racon results before your job finishes, you can copy a premade file from "/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/racon/racon_round2.fna")

```bash
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats $SCRATCH/prok/results/racon/racon_round2.fna
#> stats for /mnt/SCRATCH/bio326-2023-19/prok/results/racon/racon_round2.fna
#> sum = 234221872, n = 9546, ave = 24536.13, largest = 1187476
#> N50 = 33613, n = 1620
#> N60 = 25811, n = 2417
#> N70 = 19643, n = 3465
#> N80 = 15334, n = 4811
#> N90 = 11080, n = 6602
#> N100 = 463, n = 9546
#> N_count = 0
#> Gaps = 0

```

What changes do you see before and after running Racon?
  - Is the average contig length longer?
  - Are there fewer contigs?

</td></tr></table>


### Medaka 🐟

From the polished output of two rounds of Racon, we have an assembly that is quite good, but we can make it even better. Here we'll apply one round of medaka polishing. Medaka works much the same way but uses a different internal algorithm. In this case, we're telling medaka which exact sequencing platform we used for sequencing the reads - This is to let Medaka know about the specifics of the errors that this specific platform creates.

If curious, you can read more about how to set up Medaka here: https://github.com/nanoporetech/medaka#usage

📝 Create a file named 02b_polish2-medaka.sh with the following contents, and submit the job with sbatch:

```bash
#!/bin/bash

# Define slurm parameters
#SBATCH --job-name=polish2-medaka
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task 16
#SBATCH --mem=16G
#SBATCH --output slurm-%j-%x.out.log
#SBATCH -p smallmem,hugemem,hugemem-avx2


# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="$SCRATCH/prok/results/racon/racon_round2.fna"
in_reads="$SCRATCH/prok/results/filtlong/output.fastq.gz"
out="$SCRATCH/prok/results/medaka"


medaka_consensus -t $SLURM_CPUS_PER_TASK -d $in_assembly -i $in_reads -o $out -m r1041_e82_400bps_sup_g615


```

**NB! You can copy the script from ```/mnt/courses/BIO326/PROK/scripts/BIO326_24/02b_polish2-medaka.sh```**

Wait for the medaka job to finish. It could take several hours.

<table><tr><td>


#### Assembly-stats progress check on Medaka results

We should also check how Medaka improves our assembly.


```bash
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats $SCRATCH/prok/results/medaka/consensus.fasta
#> stats for /mnt/SCRATCH/bio326-2023-19/prok/results/medaka/consensus.fasta
#> sum = 234434401, n = 9546, ave = 24558.39, largest = 1187215
#> N50 = 33599, n = 1621
#> N60 = 25856, n = 2417
#> N70 = 19651, n = 3464
#> N80 = 15333, n = 4812
#> N90 = 11103, n = 6602
#> N100 = 461, n = 9546
#> N_count = 0
#> Gaps = 0
```

Again, do you see that your assembly is improved?

</td></tr></table>

---

### Summary 

At this point you will have taken your draft assembly from Flye and polished it with both Racon (2x) and Medaka (1x). Now, we can refer to our assembly without using the term _draft_.   


### Feedback

If you have any issues running any of the commands in this tutorial, please write an email to the instructors.
