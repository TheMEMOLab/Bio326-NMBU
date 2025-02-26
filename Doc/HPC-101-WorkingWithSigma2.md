# Working with SAGA Cluster Sigma2-NIRS

This workshop was doing by help of [Sigma2-NIRS](https://documentation.sigma2.no/index.html).

### What is this?

This document is intended to be a quick reference guide on the basic usage of the SAGA Sigma2 HPC cluster. For a complete reference please refer to the full documentation of [SAGA](https://documentation.sigma2.no/hpc_machines/saga.html).

## Login into SAGA

For login open a Command-line interface (CLI) or Terminal  and type something like this. 

```bash
ssh auve@saga.sigma2.no
```
> [!Important]
> Rember to change your user name

This will ask for your password, the one you got from Sigma2 by text. Type it.  
*Even you don't see anything the password is typed.*

After the first attempt of login you will get something like:

```
Welcome to saga.sigma2.no

Documentation:      https://documentation.sigma2.no/
Support email:      support@nris.no
Request resources:  https://www.sigma2.no/apply-e-infrastructure-resources

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Latest news from:       https://opslog.sigma2.no

  o 2025-02-17: Trouble accessing the puhuri portal via MyAccessID - LUMI
  o 2025-02-10: Filesystem problems on Betzy
  o 2025-02-10: Fram queuing system is down.
  o 2025-02-10: GPFS upgrade on the NIRD service platform.
  o 2025-02-03: Filesystem update on Betzy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NOTE: The $USERWORK autocleanup period is at least 21 days and up to 42 days if
sufficient storage is available.

Project directories have backup as described in
https://documentation.sigma2.no/files_storage/backup.html

Last login: Wed Feb 12 14:25:51 2025 from 128.39.239.134
(BASICS)[auve@login-3: ~]$

```
### SAGA main configuration 

The supercomputer is named after the goddess in Norse mythology associated with wisdom. Saga is also a term for the Icelandic epic prose literature. The supercomputer, placed at NTNU in Trondheim is designed to run both sequential and parallel workloads. It was made available to users right before the start of the 2019.2 period (October 2019).

Saga is provided by Hewlett Packard Enterprise and has a computational capacity of approximately 140 million CPU hours a year.

Let's take a look into this figure: 

![Cluster](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/cluster.png)


**All users have access to the $HOME, so please DO NOT USE THE $HOME FOR STORAGE OF LARGE FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and SLURM scripts**

## SAGA HPC nodes:

**Technical details**

Main components

* 200 standard compute nodes with 40 cores and 192 GiB memory each

* 120 standard compute nodes with 52 cores and 192 GiB memory each

* 28 medium memory compute nodes, with 40 cores and 384 GiB of memory each

* 8 big memory nodes, with 3 TiB and 64 cores each

* 2 huge memory nodes, with 6 TiB and 64 cores each

* 8 GPU nodes, with 4 NVIDIA P100 GPUs and 2 CPUs with 24 cores and 384 GiB memory each

* 8 GPU nodes, with 4 NVIDIA A100 GPUs and 1 CPU with 32 cores and 1 TiB memory each

* 8 login and service nodes with 256 cores in total

* 6.5 PB high metadata performance BeeGFS scratch file system

## Running jobs in SAGA

The HPC clusters are resources that are shared between many users, and to ensure fair use everyone must do their computations by submitting jobs through a queue system (batch system) that will execute the applications on the available resources. In our case Slurm is used as workload manager and job scheduler.

When you log in to a cluster, you are logged in to a login node shared by all users. The login nodes are meant for logging in, copying files, editing, compiling, running short tests (no more than a couple of minutes), submitting jobs, checking job status, etc. If you are unsure about the basic interaction with Unix-like systems, here is a good resource to start with. Jobs started via Slurm run on the compute nodes.

>[!Warning]
> **NEVER RUN A JOB IN THE LOGIN NODE!!! THE LOGIN NODE IS ONLY FOR LOOKING AND MANAGING FILES, INSTALLING SOFTWARE AND WRITE SCRIPTS** 

There are two ways for doing this:
* Interactive Job (via SLURM)
* Schedule a Job (via SLURM)

### What is SLURM? 

[Slurm](https://slurm.schedmd.com/) is an open source and highly scalable cluster management and job scheduling system for large and small Linux clusters. As a cluster workload manager, Slurm has three key functions.

- First, it allocates access to resources (compute nodes) to users for some duration of time so they can perform work
- Second, it provides a framework for starting, executing, and monitoring work (normally a parallel job) on the set of allocated nodes
- Finally, it arbitrates contention for resources by managing a queue of pending work (from Slurm overview)
It is important to know that:

**All Slurm commands start with letter “s" (e.g. sbatch, scancel, srun, etc...)**

**Resource allocation depends on your fairshare i.e. priority in the queue, so remember not to be "greedy" when you submit a job!!!**

### SAGA Resources and information through SLURM.

If we want to know the amount of CPU, RAM and other configuration in the cluster, we can use a set of tools (commands) SLURM provides to find available resources in SAGA.

For example we can display the amount of CPUs and RAM available in each node by:

```bash
freecores|sort -V
```

## File system in SAGA:

| Environment | Space Purpose | Backedup | Default Quota |Life span |
| ----------- | ------------- | -------- | ------------- | -------- |
| $HOME	| Personal user home space that is best for small files |	YES	|20G |	Active users |
| $USERWORK |	Temporary working directory for data analyses | 	YES	|NA|	28 days |
| $PROJECTS |	Shared disk space for research projects	| YES |	On demand	Project period |
| $LOCALSCRATCH |	Local to compute node	| NO | Varies (4.5 – 16 T)	| Until a job is completed |

>[!Warning]
> Although all users have access to the $HOME, please DO NOT USE THE $HOME FOR STORAGE any FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and configuration files.

### Where can I storage large files? 

>[!Important]
> During the BIO326 Course all data, scripts, results and so must be written in
```
/cluster/projects/nn9987k
```

Let's go there and create a folder where you can work:


```bash
cd /cluster/projects/nn9987k
mkdir $USER
tree /cluster/projects/nn9987k/$USER
cd $USER
pwd
```

### Queues for different usage:

SAGA offers different ques to run, users must select the queue that fits the most to their needs:

![queue](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/Capture.PNG)

All details of these queues can be found [here](https://documentation.sigma2.no/jobs/choosing_job_types.html)


## Running an interactive job to test programs and get used to working in the cluster

The easiest way to test software and look into huge files without messing the login node and other users, is by running an **interactive** job in SAGA. This means you can "book" a compute node and type your commands directly in that node. Let's run an interactive job by the following commands:

```bash
srun \
--account=nn9987k \
--partition=normal \
--gres=localscratch:10G \
--cpus-per-task 4 \
--nodes 1 \
--mem=10G \
--time=02:00:00 \
--pty bash \
-i
```
After you will see something like this:

```
srun: job 14011180 queued and waiting for resources
srun: job 14011180 has been allocated resources
```
Check the prompt now:

```
(BASICS)[auve@c2-41: ~]$
```
You can notice that now the prompt has changed and shows the node (computer) we are running on. In this case the node ```c2-41```. Also if this is not displayed we can take advantage of the many [SLURM_environment_variables](https://slurm.schedmd.com/pdfs/summary.pdf). These are dynamic values that SLURM uses to control the computers. For example, if you would like to know what is the node you are working on and the number of CPUs requested for this job you can print the values of that by using different SLURM variables and the command "echo" followed by the name of the variable:

```bash
echo $SLURM_NODELIST
echo $SLURM_CPUS_ON_NODE
```

In the interactive jobs we can run short parsing scripts, test software with a small datasets, etc. This is super helpful for debugging, testing software, moving data, etc.

### Temporary working directory $LOCALSCRATCH, faster and more efficient jobs

Generally any software can read (data) and write (results) from any partition of the cluster user has access to (i.e. $HOME, /cluster/projects/, etc), however, I/O (reading and writing) from those locations uses network-traffic resources resulting in a high inefficiency for heavy jobs (e.g mapping reads to large genomes/metagenomes or assembling genomes). Also if multiple users are running jobs in the same way the traffic in the network, even using the BeeGFS (1-10 Gbps), makes the jobs super slow. 
To avoid this we can take advantage of the **$LOCALSCRATCH**. This is a physical hard-drive allocated in each of the computer nodes. We can migrate the data to there for faster I/O. Often, quite some efficiency can be gained by doing this.

>[!Note]
>The amount of space required in the ```/localscratch``` directory is requested by the flag ```--gres=localscratch:<amount of memory Gb>```

Let's move to that directory and check if the amount of space requested matches with the size of the disk:

```bash
cd $LOCALSCRATCH
df -h .
```

```
Filesystem                       Size  Used Avail Use% Mounted on
/dev/mapper/xcatvg-localscratch   10G     0   10G   0% /localscratch
```
### Monitoring the job:

By using the SLRUM comand ```squeue``` we can check if our interactive job is sill running:

```bash
squeue -u $USER
```

### Example of an Interactive job by running BLAST.

**Problem: We want to detect the presence of an a-amylase sequence similar to Bacteroides gramini in a new Bacteroides genome (51).**

**Solution: we can use BLAST tool to align the sequence to all the predicted proteins in the Bacteroides51 genome and look for an ortholog of the a-amylase.**

Let's enter to that directory and then copy some fasta files from the ```/cluster/projects/nn9987k/BIO326-2025/```, this is a share directory we (teachers) will use to upload data for you. In this class we are using the files from ```/cluster/projects/nn9987k/BIO326-2025/HPC101/BLASTExample``` path.

First take a look of the data

```bash
tree /cluster/projects/nn9987k/BIO326-2025/HPC101/
```

```
/cluster/projects/nn9987k/BIO326-2025/HPC101/
└── BLASTExample
    ├── amylase.Bgramini.fasta
    └── Bacteroides51.faa

2 directories, 2 files
```
*Tip: Having more than one terminal open helps to look into multiple directories faster*

As you can see there are multiple files here. Let's copy the two fasta files **.faa and .fasta** into the $LOCALSCRATCH


```
cp /cluster/projects/nn9987k/BIO326-2025/HPC101/BLASTExample/*.* $LOCALSCRATCH
cd $LOCALSCRATCH/
ls
```

```
amylase.Bgramini.fasta  Bacteroides51.faa
```

```bash
less amylase.Bgramini.fasta 

```

```console
>WP_024997086.1 alpha-amylase [Bacteroides graminisolvens]
MKRYKYWFLLLIPFLIVACSGSDDPVIEPPVVLKEGLNYSPTAPDADQELTITFKAGSTSALYNYVGDVY
VHIGVIVDGSWKYVPAEWTENISKCKMTKTADNVWSVKLSPTVRQWFASGETSIQKLGIVIRNADGSKKG
LTDDAFVSVTDSKYKPFTPAAIKYATLPAGVKEGINIVNSSTVTLVLYDKDKSGNHKDYAHVIGDFNSWK
LTNDDKSQMNRDDAAGCWWITLSGLTGTKEYAFQYYVGTAAEGATRLADAYSRKILDPDNDSYISSTTYN
EDKTYPQGAEGIVSVFKTEPDTYTWKNTAFKMKDKDDLVIYEMLLRDFTASGDLNGAKAKLSYLKSLGVN
AIELMPVQEFDGNDSWGYNPCFFFALDKAYGTDKMYKEFIDACHGEGIAVIFDVVYNHATGSHPFAKLYW
NSATNKTSAQNPWFNVDAPHPYSVFHDFNHESPLVRAFVKRNLEFLLKEYKIDGFRFDLTKGFTQKSSTE
STASAYDATRIAILKDYNSTVKTVNPSAMMILEHFCDNAEEKELANDGMYLWRNMNYAYCESAMGLPGNS
DFSGLYDTSMPMGSLVGFMESHDEERMSFKQIAYGNYTFKTSLADRMKQLKVNTAFFLTVPGPKMIWQFG
ELGYDYSIEENGRTGKKPVKWEYYDDASRKALYDTYAKLMTLRNANTELFDTSALFSWQVKGNTNWLNGR
FLTLEGGGKKLVVAGNFTNQAGSYTVTFPHTGTWYNYMTGESVSVSATNQTISIPAHEFKLFVDFQSN
```

This is the sequence of an enzyme (a-amylase) of the bacteria Bacteroides fragilis, I would like to know if an homologue of this sequence is present in the set of sequences of **Bacteroides51.faa** (Bacteroides sp. from cockroaches). The easiest way is by doing a BLAST search. But is BLAST already installed?


```bash
blastp

```

```console
bash: blastp: command not found
```

It seems blastp is not installed as a default software in SAGA.

### Conda envrionment:

Conda is an open source package management system and environment management system that runs on Windows, macOS, and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language. You can read more about conda [here](https://docs.conda.io/en/latest/).

For **BIO-326** we will use different Conda environments previously installed in SAGA.

1. Activate the Module Anaconda to load all the conda basics

```bash
module load Anaconda3/2022.10
```

Then let's configure the interactive session with the conda global environments

```bash
eval "$(conda shell.bash hook)"
```

What we are doing here is being sure Conda is loaded and then export all the conda configurations to our shell.  
**NB! Remember that the aim of this course is not to be a Linux expert so do not worry if this is a bit cryptic for you :-)** 

If everything was OK, you should now see the ```base``` conda environment loaded and the prompt shows this:

```console
(base)[auve@c2-33: 14014597]$
```

Now we can **activate** the conda environment:

```bash
conda activate /cluster/projects/nn9987k/.share/conda_environments/BLAST/
```

```console
(BLAST)[auve@c2-33: 14014597]$ 
```

>[!Note]
>The prompt will change every time we activate different conda environments.

We can then run a blast experiment:

- Let's check blast is now running:

```bash
blastn
```

```console
BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified
Please refer to the BLAST+ user manual.
```

Now it is running...


* Index a database using ```makeblastdb``` and the molecule type prot:

```bash
makeblastdb -dbtype prot -in Bacteroides51.faa

```

```console
Building a new DB, current time: 02/20/2024 17:27:49
New DB name:   /home/work/bio326-2024-1/work.dir.of.14221763/Bacteroides51.faa
New DB title:  Bacteroides51.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 4630 sequences in 0.091763 seconds.

```

Check the results by ls:

```bash

ls

```

```console
amylase.Bgramini.fasta  Bacteroides51.faa.pdb  Bacteroides51.faa.pin  Bacteroides51.faa.pot  Bacteroides51.faa.ptf
Bacteroides51.faa       Bacteroides51.faa.phr  Bacteroides51.faa.pjs  Bacteroides51.faa.psq  Bacteroides51.faa.pto
```

And now lets run the BLAST. As we want to search for protein in a protein database the command we need to use is BLASTP:

```bash
blastp -query amylase.Bgramini.fasta -db Bacteroides51.faa -dbsize 1000000000 -max_target_seqs 1 -outfmt 6 -num_threads $SLURM_CPUS_ON_NODE -out amylase.Bgramini.fasta.blastp.out
```

Take a look into the results:

```bash
less amylase.Bgramini.fasta.blastp.out
```

```console
WP_024997086.1  D0T87_RS12665   57.772  772     301     13      8       763     28      790     0.0     908
```

It seems the amylase of *B. fragilis* has a match with the D0T87_RS12665 sequence of Bacteroides51. We can corroborate this by looking into the fasta file annotation header by doing something like this:

```bash
grep D0T87_RS12665 Bacteroides51.faa

```
We found the amylase!!!

>[!Important]
> Remember to copy the results to the ```/cluster/projects/nn9987k/$USER```

```bash

rsync -aLhv amylase.Bgramini.fasta.blastp.out /cluster/projects/nn9987k/$USER/results/

```

# sbatch scripts.

Most of the time you do not use the interactive way for submiting jobs into the cluster. To submit jobs, you need to write all the instructions you want the computer to execute. This is what a script is.

SLURM can use bash or computer scripting language (e.g. perl, python, etc) base script to read the instructions. The first line #!/bin/bash, called bash shebang, are reserved words that the computer needs to read and interpret in order to launch the program. The following lines, need to start with #SLURM and these are specific instructions SLURM uses to know how many resources and other parameters the job will use while is running. In our interactive job we used some of these SLURM parameters directly in the command line

```console
-c 2 --mem=2G -p smallmem,hugemem-avx2,test -t 01:00:00
```

SLURM uses a [bash](https://www.gnu.org/software/bash/) (computer language) base script to read the instructions. The first lines, are reserved words that SLURM needs to read inorder to launch the program:

```console
-p --partition <partition-name>       --pty <software-name/path>
--mem <memory>                        --gres <general-resources>
-n --ntasks <number of tasks>         -t --time <days-hours:minutes>
-N --nodes <number-of-nodes>          -A --account <account>
-c --cpus-per-task <number-of-cpus>   -L --licenses <license>
-w --nodelist <list-of-node-names>    -J --job-name <jobname>
```

We can indicate these options by using the ```#SBATCH``` word following by any of these flag (e.g -c 2 ; means 2 CPUs).



These parameters in a SLURM script must start with a #SBATCH string and will tell the SLURM schedule the following information. The batch script may contain options preceded with "#SBATCH" before any executable commands in the script. sbatch will stop processing further #SBATCH directives once the first non-comment non-whitespace line has been reached in the script.

- Number of nodes
- Desired number of processors or jobs
- Type of partition/queue you want to use (optional)
- Memory requirement (Optional)
- Length of time you want to run the job (Each partition has a default)
- Where to write output and error files
- Name for your job while running on HPC
- Email ID to get job status (Optional)



### Creating a SBATCH job template

Let's take a look on a basic SLURM template:


```bash
#!/bin/bash


#####Everything after this will be the instructions to SLURM###
#################################################################
## Job name:
#SBATCH --job-name=MySbatchScript  #Name of the job
#
###Account
#SBATCH --account=nn9987k #Particular account for each project
## Wall time limit:
#SBATCH --time=00:05:00  #Run for 5 minutes
#
##Partition
#SBATCH --partition=normal  #Partiion submitting the job

## Other parameters:
#SBATCH --cpus-per-task 2    #Number of cpus the job will use
#SBATCH --mem=1G             #Memory RAM
#SBATCH --nodes 1        #Number of computers
#SBATCH -o slurm-%x_%j.out    #Standar output message
#SBATCH -e slurm-%x_%j.err    #Standar error message
##############################################################

######Everything below this are the job instructions######

echo "Hello $USER this is my first JOB"
echo "I am running on the NODE $SLURM_NODELIST"
echo "I am running with $SLURM_CPUS_ON_NODE cpus"

echo "Starting $SLURM_JOB_ID at"
date

sleep 10 && echo "I slept for 10 seconds" > 10.txt

echo "Ending $SLURM_JOB_ID at"
date

```

>[!Note]
>A copy of this script is storage in ```/cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/myfirstsbatch.SLURM.sh```

We can then submit the job by:

1) go to the ```/cluster/projects/nn9987k/$USER``` folder"

```bash
cd /cluster/projects/nn9987k/$USER
```
2) Use the command sbatch to launch the job:

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/myfirstsbatch.SLURM.sh
```
## Monitoring the jobs by squeue

Users can check the status of the Job by the command ```squeue -u $USER``` 

```bash
squeue -u $USER
```


```console
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          14027935    normal MySbatch     auve  R       0:01      1 c5-45
```


Now it is runing (R). 

When the job starts it produces an out and and err file ```slurm-JOBNAME-$JOB_ID.out  slurm-JOBNAME-$JOB_ID.err``` These are all the log files and possible errors.

And after 10 seconds it will finish and we should ended up with 3 files.

```console
slurm-MySbatchScript_$JOBID.err  slurm-MySbatchScript_$JOBID$.out 10.txt
```

Let's check the output of the files:

```bash
more slurm-MySbatchScript_14027935.*
```

```console
::::::::::::::
/cluster/projects/nn9987k/auve/slurm-MySbatchScript_14027935.err
::::::::::::::
::::::::::::::
/cluster/projects/nn9987k/auve/slurm-MySbatchScript_14027935.out
::::::::::::::
Starting job 14027935 on c5-45 on saga at Thu Feb 20 19:35:52 CET 2025

Hello auve this is my first JOB
I am running on the NODE c5-45
I am running with 2 cpus
Starting 14027935 at
Thu Feb 20 07:35:52 PM CET 2025
Ending 14027935 at
Thu Feb 20 07:36:02 PM CET 2025

Job 14027935 consumed 0.0 billing hours from project nn9987k.

Submitted 2025-02-20T19:35:51; waited 1.0 seconds in the queue after becoming eligible to run.

Requested wallclock time: 5.0 minutes
Elapsed wallclock time:   10.0 seconds

Job exited normally.

Task and CPU statistics:
ID              CPUs  Tasks  CPU util                Start  Elapsed  Exit status
14027935           2            0.0 %  2025-02-20T19:35:52   10.0 s  0
14027935.batch     2      1     0.1 %  2025-02-20T19:35:52   10.0 s  0

Used CPU time:   0.0 CPU seconds
Unused CPU time: 20.0 CPU seconds

Memory statistics, in MiB:
ID               Alloc   Usage
14027935        1024.0
14027935.batch  1024.0     0.0

Job 14027935 completed at Thu Feb 20 19:36:02 CET 2025
```

This job also created a text (.txt) file names 10.txt, let's take a look:

```bash
more 10.txt
```

```console
I slept for 10 seconds
```

### Debunging errors during sbatch execution:

Let's run the following script ```/cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/mysecondbath.SLURM.sh ``` that launches a job to sleep for 20 seconds and then create a text file ```20.txt```

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/mysecondbath.SLURM.sh
```
Let's list the results:

```bash
ls -lrth
```

```console
total 3.0K
-rw-rw-r-- 1 auve nn9987k    0 Feb 20 19:35 slurm-MySbatchScript_14027935.err
-rw-rw-r-- 1 auve auve    1017 Feb 20 19:36 slurm-MySbatchScript_14027935.out
-rw-rw-r-- 1 auve nn9987k   23 Feb 20 19:36 10.txt
-rw-rw-r-- 1 auve nn9987k   76 Feb 20 19:55 slurm-MySbatchScript_14027942.err
-rw-rw-r-- 1 auve nn9987k    0 Feb 20 19:55 20.txt
-rw-rw-r-- 1 auve auve    1017 Feb 20 19:55 slurm-MySbatchScript_14027942.out
```

We can notice the ```20.txt``` file is empty (value 0 in the 4th colum), we can check then the ```slurm-MySbatchScript_$JOBID.err``` to debug:

```bash
more slurm-MySbatchScript_14027942.err
```

```console
/var/tmp/slurmd/job14225128/slurm_script: line 32: ech: command not found
```

It looks like the error is in line 32, let's take a look of the slurm script:

```bash
less +32 -N /cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/mysecondbath.SLURM.sh
```

```console
32 sleep 20 && ech "I sleept for 20 seconds" > 20.txt
```

Changing the error in that line will correct the code. Let's make a copy and then change that line. 

```bash
cp /cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/mysecondbath.SLURM.sh /cluster/projects/nn9987k/$USER
sed -i.back '32s/ech/echo/' mysecondbath.SLURM.sh

```

After changing, save the script, look for the line to be corrected and if it is correct resubmit as follow:

```bash
less +32 -N mysecondbath.SLURM.sh
```

```console
32 sleep 20 && echo "I sleept for 20 seconds" > 20.txt
```

Now the command ```echo``` is well spelled let's resubmit:

```bash
sbatch mysecondbath.SLURM.sh
```

Now the ```20.txt``` file is written and has info on it:

```bash
less 20.txt
```

```console
I sleept for 20 seconds
20.txt (END)
```

### Running BLAST as a SLRUM sbatch job.

Now that we know how to create a SLRUM sbatch job we can apply the same to submit the BLAST example. The following template will do the job:

```bash
#!/bin/bash


#####Everything after this will be the instructions to SLURM###
#################################################################
## Job name:
#SBATCH --job-name=BLASTP  #Name of the job
#
###Account
#SBATCH --account=nn9987k #Particular account for each project
## Wall time limit:
#SBATCH --time=00:10:00  #Run for 10 minutes
#
##Partition
#SBATCH --partition=normal  #Partiion submitting the job

## Other parameters:
#SBATCH --cpus-per-task 4    #Number of cpus the job will use
#SBATCH --mem=10G             #Memory RAM
#SBATCH --nodes 1        #Number of computers
#SBATCH -o slurm-%x_%j.out    #Standar output message
#SBATCH -e slurm-%x_%j.err    #Standar error message
##############################################################

######Everything below this are the job instructions######


##Activate conda environments

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2022.10


##Activate conda environments

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /cluster/projects/nn9987k/.share/conda_environments/BLAST
echo "I am workung with this" $CONDA_PREFIX

####Do some work:########

echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID\_$SLURM_ARRAY_TASK_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "I am working with this enviroment loaded"
echo $CONDA_PREFIX
echo "Today is:"
date

## Copying data to local node for faster computation

cp /cluster/projects/nn9987k/BIO326-2025/HPC101/BLASTExample/*.* $LOCALSCRATCH
cd $LOCALSCRATCH/
ls

#Index a database using makeblastdb and the molecule type prot:

echo "Indexing database..."
makeblastdb -dbtype prot -in Bacteroides51.faa

#run the BLAST. As we want to search for protein in a protein database the command we need to use is BLASTP:

echo "Running blast"

blastp \
-query amylase.Bgramini.fasta \
-db Bacteroides51.faa \
-dbsize 1000000000 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads $SLURM_CPUS_ON_NODE \
-out amylase.Bgramini.fasta.blastp.sbatch.out

#Remember to copy the results to the ```/cluster/projects/nn9987k/$USER```

echo "Copy results..."

rsync -aLhv amylase.Bgramini.fasta.blastp.sbatch.out /cluster/projects/nn9987k/$USER/results/

##Done

echo "Done Bye :-)"

date

```

We can submit the script by:

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/blast.SLRUM.sh
```

>[!Note]
> Discuss about the results, does this work?

## SFTP Guide: Transferring Files to and from the Server

This guide will show you how to use SFTP (Secure File Transfer Protocol) to upload and download files between your local computer and the remote server `login.saga.sigma2.no`. We will use the common project directory: `/cluster/project/nn9987k/$USER`.

## 1. Connecting to the Server

To establish an SFTP connection, use the following command:

```bash
sftp -P 12 <your_username>@login.saga.sigma2.no
```

Replace `<your_username>` with your actual username.

Once connected, you will see an `sftp>` prompt, indicating that you can start transferring files.

## 2. Uploading Files to the Server

### a) Upload a Single File

To upload a file from your local computer to the common project directory:

```bash
put /path/to/local/file /cluster/project/nn9987k/$USER/
```

Example:

```bash
put mydata.txt /cluster/project/nn9987k/$USER/
```

### b) Upload Multiple Files

You can also upload multiple files at once:

```bash
mput file1.txt file2.txt /cluster/project/nn9987k/$USER/
```

### c) Upload a Directory

To upload a directory and its contents, use:

```bash
put -r /path/to/local/directory /cluster/project/nn9987k/$USER/
```

## 3. Downloading Files from the Server

### a) Download a Single File

To download a file from the server to your local machine:

```bash
get /cluster/project/nn9987k/$USER/filename /path/to/local/destination/
```

Example:

```bash
get /cluster/project/nn9987k/$USER/mydata.txt ./
```

(This will save the file in the current local directory.)

### b) Download Multiple Files

To download multiple files:

```bash
mget /cluster/project/nn9987k/$USER/file1.txt /cluster/project/nn9987k/$USER/file2.txt .
```

### c) Download a Directory

To download an entire directory:

```bash
get -r /cluster/project/nn9987k/$USER/myfolder ./
```

## 4. Checking and Navigating Directories

### a) List Files

To see the files in the current remote directory:

```bash
ls
```

### b) Change Directory

To move into another directory on the server:

```bash
cd /cluster/project/nn9987k/$USER/
```

To check which directory you are currently in:

```bash
pwd
```

## 5. Exiting SFTP

When you are done, type:

```bash
bye
```

or

```bash
exit
```

This will close the SFTP session.

## Additional Tips

- Use `lcd` to change the local directory before downloading files:
  ```bash
  lcd /path/to/local/destination/
  ```
- Use `!` to run local commands without exiting SFTP (e.g., `!pwd` to check your local directory).

This guide covers the basics of using SFTP to transfer files between your computer and the remote server. Happy coding!

## Remarks

* Do not use the login node to run process (e.g. BLAST, SPADES, HMMER) **Do not get naked in the lobby**.
* Do not use the $HOME partition for storag.
* Always use the /cluster/projects/ to work.
* Use interactive jobs for testing and debugging.
* Use the $LOCALSCRATCH for faster computation.
* Monitoring your jobs by squeue.
* Use sbatch command to submit your "final" jobs scripts.

## Welcome to the world of HPC environment:
![Dali](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/DaliHPC.png)
