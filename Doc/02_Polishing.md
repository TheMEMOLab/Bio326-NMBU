# Prokaryota dry lab: 02 Polishing
### Based on ONT (Oxford Nanopore Technologies) long read sequencing of cow rumen samples
#### Wednesday 12th of April 2023


## A note about the problems that we had last time 

Last time we ran into some issues with filtlong because the locale settings on your user accounts are not set up correctly. This is of course not your fault, but we think that we solved the issue, and we just need you to run a few commands. Hopefully we will not get into more trouble.

Once you are logged into orion, please copy and paste these commands.
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

For the analysis today we need the assembly from flye. We will account on it to be in the directory that is specified in the polishing scripts, and before we start, you should check that it is in fact present and is not an empty file.

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



For the polishing exercise today we require the assembly from Flye to be present in a specific directory. You can make sure that this file is present by running the `ls` command:

```bash
ls -lh $SCRATCH/prok/results/flye/assembly.fasta
#> -rw-rw-r-- 1 cako nobody 233M Mar 22 15:16 /mnt/SCRATCH/cako/prok/results/flye/assembly.fasta

```


If see that your file either does not exist, or does not have a size comparable to around 233M (Megabytes) that is shown in the output above, you can copy the assembly file we made for you in the demo directory using the following command:


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
#SBATCH --output slurm-%x-%j.out.log


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

#### Assembly-stats progress check on Racon results

Before we continue we should check how the Racon polisher has changed our assembly. 


```bash
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats $SCRATCH/prok/results/racon/racon_round2.fna
??
```

If you compare the assembly-stats statistics before and after running Racon, what changes do you see in your assembly?


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
#SBATCH --output slurm-%x-%j.out.log


# Activate the conda environment
source activate /mnt/courses/BIO326/PROK/condaenv

# Define paths
in_assembly="results/racon/racon_round2.fna"
in_reads="results/filtlong/output.fastq.gz"
out="results/medaka"


medaka_consensus -t $SLURM_CPUS_PER_TASK -d $in_assembly -i $in_reads -o $out -m r1041_e82_400bps_sup_g615


```



#### Assembly-stats progress check on Medaka results

We should also check how Medaka enhances our assembly.


```bash
/mnt/courses/BIO326/PROK/condaenv/bin/assembly-stats $SCRATCH/prok/results/medaka/consensus.fasta
??
```

Do you see fewer contigs now, that what you had with your initial Flye assembly?


---

#### Summary 

At this point you have taken your draft assembly from Flye and polished it with both Racon (2x) and Medaka (1x). Now, we can refer to our assembly without using the term _draft_. 