## Submit a BLAST job to compare sequences but using a SLURM script.

Let's use the following SLURM script to run the BLAST as we did in the interactive job.


```bash
#!/bin/bash

#####################SLURM flags
## Job name:
#SBATCH --job-name=MyFirstBlastp
## Wall time limit:
#SBATCH --time=00:10:00
## Other parameters:
#SBATCH --cpus-per-task 2  #Asking for 2 CPUS
#SBATCH --mem=2G     #2G of RAM
#SBATCH --nodes 1 #Only one node 
#SBATCH --partition=smallmem  #partition to run
#SBATCH --output=slurm-%x_%A.out  #slurm standar output file
#########################################

######Everything below this are the job instructions######

#Loading modules
module purge #This remove any module loaded 
module load Miniconda3 && eval "$(conda shell.bash hook)"

##Activating conda environments
conda activate $COURSES/BIO326/BestPracticesOrion/BLASTConda
echo "Working with this $CONDA_PREFIX environmet ..."

##Useful lines to know where and when the job starts
##Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID\_$SLURM_ARRAY_TASK_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation

cd $TMPDIR

#Check if $USER exists in $TMPDIR

if [[ -d $USER ]]
        then
                echo "$USER exists on $TMPDIR"
        else
                mkdir $USER
fi


echo "copying files to" $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID

##Enter to the $TMPDIR/$USER

cd $TMPDIR/$USER

##Create a work directory and enter to it

mkdir work.dir.of.$SLURM_JOB_ID 
cd work.dir.of.$SLURM_JOB_ID

##Copy the fasta files form the $SCRATCH dir

echo "Copy data ..." ##Legend to know what the job is doing 

cp /mnt/courses/BIO326/BestPracticesOrion/BLASTExample/*.fa* .

###Create a protein blast database ##

echo "Making database" ##Legend to know what the job is doing

time makeblastdb \
-dbtype prot \
-in Bacteroides51.faa 

###Run BLASTp##

echo "Running BLAST" ##Legend to know what the job is doing

time blastp \
-query amylase.Bgramini.fasta \
-db Bacteroides51.faa -dbsize 1000000000 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads $SLURM_CPUS_ON_NODE \
-out amylase.Bgramini.fasta.blastp.out

###Copy results to the $SCRATCH##

echo "Copy data to the $SCRATCH ..." ##Legend to know what the job is doing
mkdir $SCRATCH/MyFirstBlastp.dir

cp *fasta.blastp.out $SCRATCH/MyFirstBlastp.dir  

echo "Your files are in $SCRATCH/MyFirstBlastp.dir"
ls $SCRATCH/MyFirstBlastp.dir

###Remove the work.directory

cd $TMPDIR/$USER

rm -rf work.dir.of.*

echo "I am done at" ##Legend to know what the job is doing
date
```

**You can copy this script to your $SCRATCH or $HOME directory from ```$COURSES/BIO326/BestPracticesOrion/myfirstblast.SLURM.sh**

```
[bio326-21-0@login bio326-21-0]$ cp $COURSES/BIO326/BestPracticesOrion/myfirstblast.SLURM.sh .
```

### Running the Job by sbatch

The way that SLURM takes the bash script and submit a job is by using the SLURM **sbatch** comand following by the script we want to run:

```bash
sbatch myfirstblast.SLURM.sh
Submitted batch job 11818971
```
The job now is in the **queue** to run.



```bash
ls
```

We can check into this file:


```bash
more slurm-MyFirstBlastp_11818971
```
```console
Working with this /mnt/courses/BIO326/BestPracticesOrion/BLASTConda environmet ...
Hello bio326-2023-19
my submit directory is:
/net/fs-2/scale/OrionStore/Home/bio326-2023-19
this is the job:
11818971_
I am running on:
cn-12
I am running with:
2 cpus
Today is:
Tue Feb 21 21:50:36 CET 2023
bio326-2023-19 exists on /home/work
copying files to /home/work/bio326-2023-19/tmpDir_of.11818971
Copy data ...
Making database


Building a new DB, current time: 02/21/2023 21:50:37
New DB name:   /home/work/bio326-2023-19/work.dir.of.11818971/Bacteroides51.faa
New DB title:  Bacteroides51.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 4630 sequences in 0.178206 seconds.



real    0m0.273s
user    0m0.213s
sys     0m0.018s
Running BLAST
Warning: [blastp] Examining 5 or more matches is recommended

real    0m0.199s
user    0m0.184s
sys     0m0.025s
Copy data to the /mnt/SCRATCH/bio326-2023-19 ...
Your files are in /mnt/SCRATCH/bio326-2023-19/MyFirstBlastp.dir
amylase.Bgramini.fasta.blastp.out
I am done at
Tue Feb 21 21:50:37 CET 2023
```

As you can see it seems the Job ran smoothly and produced the result:

```bash
ls /mnt/SCRATCH/bio326-2023-19/MyFirstBlastp.dir 
```
```console
amylase.Bgramini.fasta.blastp.out
```
