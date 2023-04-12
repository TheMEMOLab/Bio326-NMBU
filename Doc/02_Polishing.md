# Prokaryota dry lab: 02 Polishing
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Wednesday 12th of April 2023


## A note about the problems that we had last time 

Last time we ran into some issues with filtlong because the locale settings on your user accounts were not set up correctly from the beginning. This is of course not at your fault, but we think that we solved the issue, and we just need you to run a few commands to fix it. Hopefully we will not get into more trouble.

Once you are logged into orion, please copy and paste these commands:

```bash
echo """ 
export LC_ALL=C; unset LANGUAGE # Fixes a bug in filtlong
if [ -f ~/.bashrc ]; then 
  . ~/.bashrc
fi
export PS1=\"\\[\\e[32m\\]\\u\\[\\e[m\\]\\[\\e[33m\\]@\\[\\e[m\\]\\[\\e[34m\\]\\h\\[\\e[m\\] \\[\\e[35m\\]\\W\\[\\e[m\\] \\[\\e[33m\\]\\\\$\\[\\e[m\\] \"
""" >> ~/.bash_profile # Sets up a bash profile for your terminal.
exit

```

After you run these commands on orion, it should automatically log you out, and you can log in again. When you log in again, you should see your terminal having new colors:




<img width="346" alt="Screenshot 2023-03-22 at 15 41 45" src="https://user-images.githubusercontent.com/5913696/226940379-933b380a-e6a8-4af9-873c-202c799e435b.png">


---

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


For the polishing exercise today we require the assembly from Flye to be present in a predefined directory. You can make sure that this file is present by running the `ls` command:

```bash
ls -lh $SCRATCH/prok/results/flye/assembly.fasta
#> -rw-rw-r-- 1 cako nobody 233M Mar 22 15:16 /mnt/SCRATCH/cako/prok/results/flye/assembly.fasta
```

If you see that your file either does not exist or does not have a size comparable to around 233M (megabytes) which is exemplified in the ls output above, you can copy the assembly file we made for you in the demo directory using the following command:


```bash
cp -v /mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly* $SCRATCH/prok/results/flye/
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly.fasta' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly.fasta'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_graph.gfa' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_graph.gfa'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_graph.gv' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_graph.gv'
#> '/mnt/courses/BIO326/PROK/data/metagenomic_assembly/demo/results/flye/assembly_info.txt' -> '/mnt/SCRATCH/cako/prok/results/flye/assembly_info.txt'
```

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

Wait for the Racon job to finish. It should take less than an hour.


<table><tr><td>

#### Assembly-stats progress check on Racon results

Let's check how Racon has improved our Flye assembly.


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

If you compare the assembly-stats statistics before and after running Racon, what changes do you see in your assembly?
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
#SBATCH --cpus-per-task 8
#SBATCH --mem=16G
#SBATCH --output slurm-%j-%x.out.log


# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="$SCRATCH/prok/results/racon/racon_round2.fna"
in_reads="$SCRATCH/prok/results/filtlong/output.fastq.gz"
out="$SCRATCH/prok/results/medaka"


medaka_consensus -t $SLURM_CPUS_PER_TASK -d $in_assembly -i $in_reads -o $out -m r1041_e82_400bps_sup_g615


```

Wait for the medaka job to finish. It could take several hours.

<table><tr><td>


#### Assembly-stats progress check on Medaka results

We should also check how Medaka improves our assembly.


```bash
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats $SCRATCH/prok/results/medaka/consensus.fasta
#> ??
```

Do you see fewer contigs now, than what you had with your initial Flye assembly?

</td></tr></table>

---

### Summary 

At this point you will have taken your draft assembly from Flye and polished it with both Racon (2x) and Medaka (1x). Now, we can refer to our assembly without using the term _draft_.   
