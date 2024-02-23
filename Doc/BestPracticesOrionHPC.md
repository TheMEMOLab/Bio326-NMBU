# Basic exercises for best practices in Orion
v4, February, 2024 Author: Arturo Vera-Ponce de Leon

### What is this?
This document is intended to be a quick reference guide on the basic usage of the NMBU-Orion HPC cluster. For a complete reference please referer to the full documentation of [Orion](https://orion.nmbu.no/)

**Login into orion** 

To login into Orion cluster we need two things:
- Establish a [VPN](https://nmbuhjelp.nmbu.no/tas/public/ssp/content/detail/knowledgeitem?unid=4c4870bb-f0fc-4c9f-b0f5-0c036951436f) connection (if not using a NMBU network)
- Use a secure-shell command [ssh](https://en.wikipedia.org/wiki/SSH_(Secure_Shell)) via the command line

For login open a Command-line interfase (CLI) or Terminal  and type something like this. 

```bash
ssh bio326-2024-1@login.orion.nmbu.no
```
*Remember to change to your username bio326-y-x*

This will ask for your password. Type it

*Even you don't see anything the password is typed*

After the first attempt of login you will get something like:

```console
bio326-2024-1@login.orion.nmbu.no's password:
You are required to change your password immediately (root enforced)
Last login: Tue Feb 13 13:37:59 2024 from 10.42.16.74
WARNING: Your password has expired.
You must change your password now and login again!
Changing password for user bio326-2024-1.
Changing password for bio326-2024-1.
(current) UNIX password:
```

Type the password again and then create a new password.

**!NB If is your first time using the CLI it is strongly suggested to use a simple password you don't forget!** 

```console
New password:
Retype new password:
passwd: all authentication tokens updated successfully.
Connection to login.orion.nmbu.no closed.
```

If this was done correctly, you will be logout and need to relogin:

```bash
ssh bio326-2024-1@login.orion.nmbu.no

```

```console
bio326-2024-1@login.orion.nmbu.no's password:
```

Then this message will apear:

```console
Last login: Tue Feb 20 14:47:50 2024 from 10.42.26.49


Welcome to the NMBU Orion compute cluster environment.

You have successfully logged in to a machine that grants you access to your home directory, enables script editing, file management, and job submission within the cluster environment. Please refrain from running any jobs on this machine, as they may be automatically terminated..




Your home directory quota is:
Volume          Type    ID                      Used    Quota   Limit   Grace
Home            USR     bio326-2024-1           16K     200G    300G    none


NEWS:
- 2023-07-12: We upgraded Slurm, JupyterHub, RStudio and R to the latest version. Please email us if you miss anything or notice any issues.

IMPORTANT:
- To learn more about Orion, please visit our introduction page: https://orion.nmbu.no/
- If your project requires additional CPU hours, we recommend applying for national infrastructure resources at https://www.sigma2.no/
- For optimal file management, please use filemanager.orion.nmbu.no and we strongly advise compressing your fastq, vcf, and other non-compressed files using tools like pigz.
- Please always consider using the $TMPDIR storage, which is local to each compute node and offers faster read and write speeds. More information about this and other Orion file systems can be found at https://orion.nmbu.no/en/filesystems and https://orion.nmbu.no/en/storagepolicy. If you are unsure how to use $TMPDIR, please contact us.

For any inquiries related to Orion, please contact us at orion-support@nmbu.no. You can also find us on Teams: https://bit.ly/orion-teams


Press Enter to acknowledge that you read the message above and continue.
```


After pressing enter you will see this:

```console
            ,        ,
            /(        )`
            \ \___   / |
            /- _  `-/  '
           (/\/ \ \   /\
           / /   | `    \
           O O   ) /    |
           `-^--'`<     '  Linux Version 3.10.0-1160.42.2.el7.x86_64
          (_.)  _  )   /  Compiled #1 SMP Tue Sep 7 14:49:57 UTC 2021
           `.___/`    /  Eight 2.19GHz Intel Intel(R) Xeon(R) Platinum 8352Y CPU @ 2.20GHz Processors, 32GB RAM
             `-----' /  35118 Bogomips Total
<----.     __ / __   \  Load Average 1.30, 1.21, 1.23
<----|====O)))==) \) /====  login
<----'    `--' `.__,' \
             |        |
              \       /       /\
         ______( (_  / \______/
       ,'  ,-----'   |
       `--{__________)
Hello, bio326-2024-1! Have a fantastic afternoon!

-bash-4.2$
```

**Now you are logged into the Orion login-node.**

This a very basic configuration of the [promt](https://en.wikibooks.org/wiki/Guide_to_Unix/Explanations/Shell_Prompt), so let's changing to be more colorful and personalized. Just copy the following command:

```bash

rsync -aPlvh --no-owner --no-group /mnt/courses/BIO326/BestPracticesOrion/.bash* . && source .bash_profile

```

You will see something like this after:

```console
[bio326-2024-1@login: ~]$

```


### Orion main configuration 

**The Orion HPC system currently consists of 1680 processor cores with more than 12 terabytes of RAM and 1 petabyte of storage accessible on a 10/1 Gbit network via NFS. The Orion HPC's computation nodes utilize various processors with Centos Linux 7.9 as an operating environment.**

Let's take a look into this figure: 

![Cluster](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/cluster.png)

**NEVER RUN A JOB IN THE LOGIN NODE!!! THE LOGIN NODE IS ONLY FOR LOOKING AND MANAGING FILES, INSTALL SOFTWARE AND WRITE SCRIPTS** 

How can I be sure of the number of CPUs and RAM of this "login" computer node and other nodes?

* CPUS: Use the command nproc

```bash
nproc
```


* RAM: We need to look for the "Total memory". All this info is allocated in the meminfo file at /proc directory. So we can use the grep command to look for this into the file.

```bash
grep MemTotal /proc/meminfo |awk '{print $1,$2/1000000 " GB"}'

```

```console
MemTotal: 32.7794 GB
```

As you can see, this computer is not well suitable for "heavy" computational work. So if we want to do some work (e.g. run BLAST or assembly a genome) we need to send this (job) into a compute node.

### The Orion HPC nodes:

**Orion HPC provides a supercomputing environment with thousands of cores to researchers and students to handle their computational problems. The following table summarizes the computing power capacity of Orion HPC nodes.**

|     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- |--- |
| Number of nodes | RAM(GB)\* | CPU type | Clock rate (GHz) | Cores \*\* | $TMPDIR(TB)\*\*\* | Nodes | Remark|
| 1   | 3000 | Xeon(R) Gold 6140 | 2.30 | 144 | 2   | cn-1    ||
| 7   | 192 | Xeon(R) CPU E5-2650 v2 | 2.60 | 32  | 4.5 | cn[4-11]    ||
| 1   | 256 | Xeon(R) CPU E5-2683 v4 | 2.10 | 64  | 4.5 | cn-13    ||
| 2   | 1000 | Xeon(R) CPU E7- 4870 | 2.40 | 80  | 18/8   | cn-[2-3]    ||
| 4   | 256 | EPYC 7302 16-Core | 3.0 | 64  | 16  | gn-[0-3] | 3 x NVIDIA Quadro RTX 8000|
| 2   | 256 | Xeon(R) CPU E5-2650 v2 | 2.60 | 32  | 11   |cn-12,15     ||
| 1   | 2000 | EPYC 7742 64-Core | 2.25 - 3.40 | 256 | 8   | cn-14|    |
| 2   | 1000 | EPYC 7702 64-Core | 2 - 3.35 | 256 | 2   |  cn-[16-17]   ||

There are two ways for doing this:
* Interactive Job (via SLURM)
* Schedule a Job (via SLURM)

### What is SLURM? 

[Slurm](https://slurm.schedmd.com/) is an open source and highly scalable cluster management and job scheduling system for large and small Linux clusters. As a cluster workload manager, Slurm has three key functions.

- First, it allocates access to resources (compute nodes) to users for some duration of time so they can perform work
- Second, it provides a framework for starting, executing, and monitoring work (normally a parallel job) on the set of allocated nodes
- Finally, it arbitrates contention for resources by managing a queue of pending work (from Slurm overview)
It is important to know that:

**All Slurm commands start with letter “s" (e.g sbatch, scancel, srun, etc...)**

**Resource allocation depends on your fairshare i.e. priority in the queue, so remember not to be "greedy" when you submit a job!!!**

### Orion Resources and information through SLURM.

If we want to know the amount of CPU, RAM and other configuration in the cluster, we can use a set of tools (commands) SLURM provide to find available resources in Orion.

For example we can display the Partition, No. of CPUs, Memmory ammount of each node (computer) in Orion using the following instructions:

```bash
sinfo -l -N
```

This will display this table:

```console
Tue Feb 20 15:50:09 2024
NODELIST   NODES    PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON
cn-1           1         test   allocated 144    4:18:2 309453        0      1 cpu_xeon none
cn-2           1       orion*       mixed 80     40:1:2 103186        0      1 cpu_xeon none
cn-2           1      hugemem       mixed 80     40:1:2 103186        0      1 cpu_xeon none
cn-3           1       orion*       mixed 80     40:1:2 967353        0      1 cpu_xeon none
cn-3           1      hugemem       mixed 80     40:1:2 967353        0      1 cpu_xeon none
cn-4           1      RStudio       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-5           1       orion*   allocated 32     32:1:1 193230        0      1 cpu_xeon none
cn-5           1     smallmem   allocated 32     32:1:1 193230        0      1 cpu_xeon none
cn-6           1       orion*   allocated 32     32:1:1 193230        0      1 cpu_xeon none
cn-6           1     smallmem   allocated 32     32:1:1 193230        0      1 cpu_xeon none
cn-7           1      RStudio       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-8           1       orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-8           1     smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-9           1       orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-9           1     smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-10          1       orion*   allocated 32     32:1:1 193230        0      1 cpu_xeon none
cn-10          1     smallmem   allocated 32     32:1:1 193230        0      1 cpu_xeon none
cn-11          1       orion*   allocated 64     2:16:2 257687        0      1 cpu_xeon none
cn-11          1 hugemem-avx2   allocated 64     2:16:2 257687        0      1 cpu_xeon none
cn-12          1       orion*       mixed 32      2:8:2 257738        0      1 cpu_xeon none
cn-12          1      hugemem       mixed 32      2:8:2 257738        0      1 cpu_xeon none
cn-14          1       orion*       mixed 256    2:64:2 205153        0      1 cpu_amd, none
cn-14          1 hugemem-avx2       mixed 256    2:64:2 205153        0      1 cpu_amd, none
cn-15          1       orion*   allocated 32      2:8:2 257738        0      1 cpu_xeon none
cn-15          1     smallmem   allocated 32      2:8:2 257738        0      1 cpu_xeon none
cn-16          1       orion*       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-16          1 hugemem-avx2       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-17          1 hugemem-avx2       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-17          1       orion*       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-18          1 hugemem-avx2       mixed 56     2:14:2 386700        0      1 cpu_amd, none
cn-18          1       orion*       mixed 56     2:14:2 386700        0      1 cpu_amd, none
cn-19          1 hugemem-avx2       mixed 56     2:14:2 386700        0      1 cpu_amd, none
cn-19          1        test2       mixed 56     2:14:2 386700        0      1 cpu_amd, none
gn-0           1          gpu       mixed 64     2:16:2 257710        0      1 cpu_amd, none
gn-1           1          gpu       mixed 64     2:16:2 257710        0      1 cpu_amd, none
gn-2           1          gpu       mixed 64     2:16:2 257710        0      1 cpu_amd, none
gn-3           1          gpu       mixed 64     2:16:2 257710        0      1 cpu_amd, none
```

In this case, the *State* column shows the status of the node. It means, how many resources can be allocated per node, in this example there are 4 different status: 

* allocated: The node has been allocated to one or more jobs.
* completing* : All jobs associated with this node are in the process of COMPLETING. This node state will be removed when all of the job's processes have terminated
* drained: The node is unavailable for use per system administrator request. 
* mixed: The node has some of its CPUs ALLOCATED while others are IDLE.

Summarizing, the only nodes that can accept jobs under the previous conditions are those with "MIXED" status. 

## File system in Orion:

| Environment | Space Purpose | Backedup | Default Quota |Life span |
| ----------- | ------------- | -------- | ------------- | -------- |
| $HOME	| Personal user home space that is best for small files |	YES	|200G - 300G* |	Active users |
| $SCRATCH or $SCRATCH_PROJECTS |	Temporary working directory for data analyses | 	YES	| 500G - 1T*|	180 days |
| $PROJECTS |	Shared disk space for research projects	| YES |	On demand	Project period (max 3 - 4 years) |
| $TMPDIR |	Local to compute node	| NO | Varies (4.5 – 16 T)	| Until a job is completed |

As an Orion user you can check the ammount of space used in different directories of the cluster. To check all the disks and mounted partitions in Orion we can run the following command:

```bash
df -h
```

As you can notice there are plenty of directories in Orion, but let's focus in the $HOME partition:

```bash
df -h .
```

*syntaxis: d(isk)f(ree) -h(uman redable) .(the directory I am now)*

```console
Filesystem                   Size  Used Avail Use% Mounted on
fs-2:/scale/OrionStore/Home   90T   76T   15T  85% /net/fs-2/scale/OrionStore/Home
```

**All users have access to the $HOME, so please DO NOT USE THE $HOME FOR STORAGE LARGE FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and SLURM scripts**

### Where can I storage large files? 

There are two dierectories designed to this:
* $SCRATCH
* $PROJECT

As a student of this course, we are using the $SCRATCH folder to keep our raw sequencing files and final results. This directory, in contrast to the $HOME and $PROJECT, is not backed up. **Remember to make a copy of your important files into another (permanent) location!!!** All students have a directory in the $SCRATCH: /mnt/SCRATCH/bio326-y-x ; where y is the year and x is the student number.

Let's move into that folder:

```bash
cd $SCRATCH
pwd
```
*For moving here we use the variable ```$SCRATCH```. A variable is a character string to which the user assign a value. The value assigned could be a number, text, filename, device, or any other type of data. It is important you get familiarized with this and other Unix/Linux terminology. This [link](https://orion.nmbu.no/en/LinuxTutorial) can help you to navigate through these terms.

```console
/mnt/SCRATCH/bio326-2024-1
```


### Queues for different RAM usage:

Orion offers three primary job queues, or partitions: hugemem-avx2, RStudio, smallmem and gpu. You are encouraged to choose the most appropriate partition.

|Partition 	|Job memory requirements| 	Available resources|
|-----------|-----------------------|----------------------|
|hugemem-avx2 	|> 100 GB RAM 	|6 x (>1 TB RAM, 80 cores)|
|RStudio  | ~100Gb RAM|2 x (~ 100 GB RAM,32 cores) |
|smallmem |	1-100 GB RAM 	|6 x (192 GB RAM, 32 cores)|
|gpu 	|1-100 GB RAM 	|4 x (256 GB RAM, 64 cores, 3 Quadro RTX 8000)|

To select the correct partition queue we can use the **SLURM** command: *-p \<partition>* and choosing the appropriate partition depending the amount of memory RAM the job will use to compute.
In general, there is a shorter job queue in smallmem than the hugemem queue. Also, be prepared to wait for your job to start if you need to use the hugemem queue.

### Resources in the cluster

Some times the cluster is super busy, a way we can check how many resources are available is by running the ```freecores``` command:

```bash
freecores
```

```console
cn-10   0 of 32 cores free,  38606 of 193230 MB free,  0.0 PEs free (ALLOCATED)
cn-11   0 of 64 cores free, 216727 of 257687 MB free,  0.0 PEs free (ALLOCATED)
cn-15   0 of 32 cores free,  70346 of 257738 MB free,  0.0 PEs free (ALLOCATED)
cn-6    0 of 32 cores free,  86734 of 193230 MB free,  0.0 PEs free (ALLOCATED)
cn-9    0 of 32 cores free, 107214 of 193230 MB free,  0.0 PEs free (ALLOCATED)
cn-12   4 of 32 cores free, 137930 of 257738 MB free,  4.0 PEs free (MIXED)
cn-14   4 of 256 cores free, 241105 of 2051537 MB free,  4.0 PEs free (MIXED)
cn-17   6 of 256 cores free, 140882 of 1031762 MB free,  6.0 PEs free (MIXED)
cn-3    8 of 80 cores free, 783321 of 967353 MB free,  8.0 PEs free (MIXED)
cn-5   10 of 32 cores free,  82638 of 193230 MB free, 10.0 PEs free (MIXED)
cn-2   12 of 80 cores free, 741045 of 1031861 MB free, 12.0 PEs free (MIXED)
cn-7   16 of 32 cores free, 129510 of 193230 MB free, 16.0 PEs free (MIXED)
cn-8   20 of 32 cores free, 119502 of 193230 MB free, 19.8 PEs free (MIXED)
gn-0   20 of 64 cores free,  33582 of 257710 MB free,  8.3 PEs free (MIXED)
cn-18  26 of 56 cores free,  41612 of 386700 MB free,  6.0 PEs free (MIXED)
cn-4   30 of 32 cores free, 152270 of 193230 MB free, 25.2 PEs free (MIXED)
gn-2   36 of 64 cores free,  81582 of 257710 MB free, 20.3 PEs free (MIXED)
gn-3   36 of 64 cores free,  81582 of 257710 MB free, 20.3 PEs free (MIXED)
cn-19  46 of 56 cores free, 233100 of 386700 MB free, 33.8 PEs free (MIXED)
gn-1   48 of 64 cores free, 155310 of 257710 MB free, 38.6 PEs free (MIXED)
cn-16  68 of 256 cores free, 228946 of 1031762 MB free, 56.8 PEs free (MIXED)
Total   390 of  1648 cores free,  3812 of  9808 GB free, 289.0 PEs free

```

Although this is useful, the information is not sorted and it is difficult to read, let's sort it a bit:

```bash
freecores |sort -V
```

```console
Total   390 of  1648 cores free,  3812 of  9808 GB free, 289.0 PEs free
cn-2   12 of 80 cores free, 741045 of 1031861 MB free, 12.0 PEs free (MIXED)
cn-3    8 of 80 cores free, 783321 of 967353 MB free,  8.0 PEs free (MIXED)
cn-4   30 of 32 cores free, 152270 of 193230 MB free, 25.2 PEs free (MIXED)
cn-5   10 of 32 cores free,  82638 of 193230 MB free, 10.0 PEs free (MIXED)
cn-6    0 of 32 cores free,  86734 of 193230 MB free,  0.0 PEs free (ALLOCATED)
cn-7   16 of 32 cores free, 129510 of 193230 MB free, 16.0 PEs free (MIXED)
cn-8   20 of 32 cores free, 119502 of 193230 MB free, 19.8 PEs free (MIXED)
cn-9    0 of 32 cores free, 107214 of 193230 MB free,  0.0 PEs free (ALLOCATED)
cn-10   0 of 32 cores free,  38606 of 193230 MB free,  0.0 PEs free (ALLOCATED)
cn-11   0 of 64 cores free, 216727 of 257687 MB free,  0.0 PEs free (ALLOCATED)
cn-12   4 of 32 cores free, 137930 of 257738 MB free,  4.0 PEs free (MIXED)
cn-14   4 of 256 cores free, 241105 of 2051537 MB free,  4.0 PEs free (MIXED)
cn-15   0 of 32 cores free,  70346 of 257738 MB free,  0.0 PEs free (ALLOCATED)
cn-16  68 of 256 cores free, 228946 of 1031762 MB free, 56.8 PEs free (MIXED)
cn-17   6 of 256 cores free, 140882 of 1031762 MB free,  6.0 PEs free (MIXED)
cn-18  26 of 56 cores free,  41612 of 386700 MB free,  6.0 PEs free (MIXED)
cn-19  46 of 56 cores free, 233100 of 386700 MB free, 33.8 PEs free (MIXED)
gn-0   20 of 64 cores free,  33582 of 257710 MB free,  8.3 PEs free (MIXED)
gn-1   48 of 64 cores free, 155310 of 257710 MB free, 38.6 PEs free (MIXED)
gn-2   36 of 64 cores free,  81582 of 257710 MB free, 20.3 PEs free (MIXED)
gn-3   36 of 64 cores free,  81582 of 257710 MB free, 20.3 PEs free (MIXED)

```

By using this we can plan our experiment and see how much resources are available at the moment.

## Running an interactive job to test programs and get used to working in the cluster

The easiest way to test software and look into huge files without messing the login node and other users, is by running an **interactive** job in Orion. This means you can "book" a compute node and type your commands directly in that node. Let's run an interactive job by the following commands:

```bash
qlogin -c 2 --mem=2G -p smallmem,hugemem-avx2,test -t 01:00:00
```

*Basic syntaxis the command:
 qlogin \<slurm-options> \<software-name/path>*
  
It might take a while to SLURM allocate the resources of this job. But as soon as it allocates the job a message like this will be displayed:

```console
salloc: Pending job allocation 14221700
salloc: job 14221700 queued and waiting for resources
salloc: job 14221700 has been allocated resources
salloc: Granted job allocation 14221700
```


You can notice that now the prompt has changed and shows the node (computer) we are running on. In this case the node "cn-x". Also if this is not displayed we can take advantage of the many [SLURM_environment_variables](https://slurm.schedmd.com/pdfs/summary.pdf). These are dynamic values that SLURM uses to control the computers. For example, if you would like to know what is the node I am working on and the number of CPUs requested for this job you can print the values of that by using different SLURM variables and the command "echo" follows by the name of the variable:

```bash
echo $SLURM_NODELIST
echo $SLURM_CPUS_ON_NODE
```

In the interactive jobs we can run short parsing scripts, test software with a small datasets, etc. This is super helful for debuging or testing software.

### Temporary working directory, faster and more efficient Jobs

Generaly any software can read (data) and write (results) from any partition of the cluster (i.e. $HOME, $SCRATCH, $PROJECT), however, I/O (reading and writing) from those locations uses a lot of network-trafic resources resulting in a high inefficenfy for heavy jobs (e.g mapping reads to large genomes/metagenomes or assembly genomes). Also if multiple users are running jobs in the same way the traffic in the network, even using the InfiniBand (1-10 Gbps), makes the jobs super slow. 
To avoid this we can take advantage of the **$TMPDIR** partition. This is a physical hard-drive allocated in each of the computer nodes. We can migrate the data to there for faster I/O. Often, quite some efficiency can be gained by doing this.

Let's take a look, first we need to check if our **$USER** exists in that **$TMPDIR**

```bash
echo $TMPDIR/$USER
/home/work/bio326-2024-1
```

This means the user **bio326-y-u** has a directory in the **$TMPDIR** (/home/work). Move to that directory:

```bash
cd $TMPDIR/$USER

```

It could happen that our $USER name is not in the $TMPDIR, if that is the case we can easily create this directory and go there:

```bash
 mkdir $TMPDIR/$USER
 cd $TMPDIR/$USER

```
Let's check how much space we have available here:

```bash
df -h .

```

```console
Filesystem               Size  Used Avail Use% Mounted on
/dev/mapper/centos-home  4.4T   33G  4.3T   1% /home
```

Now we need to create an other directory **a work directory** to copy data for executing some commands. We can use another SLURM variable, let's say the JOBID to be consistent.


```bash
mkdir work.dir.of.$SLURM_JOB_ID
ls
cd work.dir.of.$SLURM_JOB_ID
```

Then check we have been done this correctly:

```bash
pwd
```

```console
/mnt/SCRATCH/bio326-2024-1/work.dir.of.14221700
```

By using the $SLURM_JOB_ID we can further identify what job we are running.

### Monitoring the job:

By using the SLRUM comand ```squeue``` we can check if our interactive job is sill running:

```bash
squeue -u $USER
```

### Example of an Interactive job by running BLAST.

**Problem: We want to detect the presence of a a-amylase sequence similar to Bacteroides gramini in a new Bacteroides genome (51).**

**Solution: we can use BLAST tool to align the sequence to all the predicted proteins in the Bacteroides51 genome and look for an ortholg of the a-amylase.**

Let's enter to that directory and then copy some fasta files from the ```$COURSES/BIO326```, this is a share directory we (teachers) will use to upload data for you. In this class we are using the files from ```$COURSES/BIO326/BestPracticesOrion/BLASTExample``` path.

First take a look of the data

```bash
ls -l $COURSES/BIO326/BestPracticesOrion/BLASTExample
```

```console
total 2049
-rwxrwxr-x 1 auve bio326     838 Nov  7 23:16 amylase.Bgramini.fasta
-rwxrwxr-x 1 auve bio326 2085506 Nov  7 23:06 Bacteroides51.faa
```

*Tip: Having more than one terminal open does help to faster look into multiple directories*

As you can see there are multiple files here, lets copy the two fasta files **.faa and .fasta** into the $TMPDIR/$USER/work.dir.of.$SLURM_JOB_ID

```bash
cp /mnt/courses/BIO326/BestPracticesOrion/BLASTExample/*.fa* .
ls

```

```console
amylase.Bgramini.fasta  Bacteroides51.faa
```

*Remember that you can copy multiple files using regular expression (REGEX) in this case* * *.fa* * *means "everything that has .fa on it"*

No we can do some work on this files. Take a look of the **amylase.Bgramini.fasta** file 

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

It seems blastp is not installed as a default software in Orion.

### Conda envrionment:

Conda is an open source package management system and environment management system that runs on Windows, macOS, and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language. You can read more about conda [here](https://docs.conda.io/en/latest/).

For **BIO-326** we will use different Conda environment previously installed in Orion. However, as a Orion user you are allowed to install your own environment, please refere to the [Orion-Conda Environment](https://orion.nmbu.no/en/CondaEnvironment) for doing this.

1. Activate the Module miniconda3 to load all the conda basics

```bash
module load Miniconda3
```

Then let's configure the interactive session with the conda global environments

```bash
eval "$(conda shell.bash hook)"
```

What we are doing here is being sure Conda is loaded and then export all the conda configurations to our shell...**NB! Remember that the aim of this course is not to be a Linux expert so do not worry if this is a bit cryptic for you :-)** 

If everything was OK, you should now see the ```base``` conda environment loaded and the prompt shows this:

```console
(base) [bio326-2024-1@cn-14: work.dir.of.14221763]$
```

Now we can **activate** the conda environment:

```bash
conda activate $COURSES/BIO326/BestPracticesOrion/BLASTConda
```

Notice again how the prompt has changed:

```console
(/mnt/courses/BIO326/BestPracticesOrion/BLASTConda) [bio326-2024-1@cn-14: work.dir.of.14221763]$
```

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

And now lets run the BLAST,as we want to search for protein in a protein database the command we need to use is BLASTP:

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

It seems the amylase of *B. fragilis* has a match wiht the D0T87_RS12665 sequence of Bacteroides51. We can corroborate this by looking into the fasta file annotation header by doing something like this:

```bash
grep D0T87_RS12665 Bacteroides51.faa

```
We found the amylase!!!

### Copy results to the $SCRATCH, remove work.directory and exit the job.

**NB! Remember the $TMPDIR is a temprary directory so, we need to move the results back to our $SCRATCH partition. For this we can use the following commands:**

- Create a directory in the ```$SCRATCH``` to save all the outputs

```bash
mkdir $SCRATCH/MyInteractiveBlast.dir
```

```bash
rsync -aPlvh --no-owner --no-group *fasta.blastp.out $SCRATCH/MyInteractiveBlast.dir/
```

```console
sending incremental file list
amylase.Bgramini.fasta.blastp.out
             68 100%    0.00kB/s    0:00:00 (xfr#1, to-chk=0/1)

sent 185 bytes  received 35 bytes  440.00 bytes/sec
total size is 68  speedup is 0.31
```

Then let's sure this is copy back

```bash
ls $SCRATCH/MyInteractiveBlast.dir

```

```console
amylase.Bgramini.fasta.blastp.out
```

Finally, as the **$TMPDIR** is used for everyone a best practice is to delete all the temporary directories (i.e work.directory) from this location.

We can achive this by doing this:

* First go back to the main $TMPDI/$USER

```bash
cd $TMPDIR/$USER
rm -rf work.dir.of$SLURM_JOB_ID
```

Finally, we can logout of this node:

```bash
exit
```

You can see now we return to the main **login bio326-y-x** node.

## sbatch scripts.

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

To run your jobs on Orion, you should create a job script and submit it.

Let's go to our $SCRATCH, create a directory names SLURMTest and there a text file using vi or nano named myfirstsbatch.SLURM.sh

```bash
cd $SCRATCH
mkdir SLURMTest
nano myfirstsbatch.SLURM.sh
```

Basic commands on nano:

```^ = Control```
Save as:

```^O = Control + O```

Exit

```^x = Control + x```


```bash
#!/bin/bash


#####Everything after this will be the instructions to SLURM###
#################################################################
## Job name:
#SBATCH --job-name=MySbatchScript  #Name of the job
#
## Wall time limit:
#SBATCH --time=00:05:00  #Run for 5 minutes
#
##Partition
#SBATCH --partition=smallmem  #Partiion submitting the job

## Other parameters:
#SBATCH --cpus-per-task 2    #Number of cpus the job will use
#SBATCH --mem=1G             #Memory RAM
#SBATCH --nodes 1        #Number of computers
#SBATCH --mail-user=arturo.vera.ponce.de.leon@nmbu.no  #User email
#SBATCH --mail-type=begin      #notify by email when job starts
#SBATCH --mail-type=end        #notify by email when job ends
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

*You can copy this script from /mnt/courses/BIO326/BestPracticesOrion/myfirstsbatch.SLURM.sh by* ``` rsync -aPL $COURSES/BIO326/BestPracticesOrion/myfirstsbatch.SLURM.sh $SCRATCH/SLURMTest/```

Then submit the script by using the command ```sbatch```

```bash
sbatch myfirstsbatch.SLURM.sh

```

## Canceling Jobs 

Some times happens that we start a job but find some bugs in the script or simply we do not want to run for any reason. In this case there is a way to **cancel** jobs.
For this, we can use the **scancel** command following the **JOBID**

For example the following job 12315677:

```bash
squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
          12315677     orion MyFirstB bio326-2 PD       0:00      1 (Priority)
 ```

To cancel just type:

```bash
scancel 12315677
```

And then check for the status:

```bash

squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
```

If no slurm.out file is created and no job is showing by the squeue command, it means the job has been canceled.

### Monitoring the jobs by squeue

Users can check the status of the Job by the command ```squeue -u $USER``` 

```bash
 squeue -u $USER
```

```console
 JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          14221857  smallmem MySbatch bio326-2 PD       0:00      1 (Priority)
```
This means your job is pending to be submitted (PD)


```bash
squeue -u $USER

```

```console
          JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          14221857  smallmem MySbatch bio326-2  R       0:05      1 cn-6     
```

Now it is runing (R). 

When the job starts it produces an out and and err file ```slurm-JOBNAME-$JOB_ID.out  slurm-JOBNAME-$JOB_ID.err``` These are all the log files and possible errors.

And after 10 seconds it will finish and we should ended up with 3 files.

```console
slurm-MySbatchScript_14221857.err  slurm-MySbatchScript_14221857.out 10.txt
```
Let's check the output of the files:

```bash
more slurm-MySbatchScript_14221857.*
```

```console
::::::::::::::
slurm-MySbatchScript_14221857.err
::::::::::::::
::::::::::::::
slurm-MySbatchScript_14221857.out
::::::::::::::
Hello bio326-2024-1 this is my first JOB
I am running on the NODE cn-6
I am running with 2 cpus
Starting 14221857 at
Tue Feb 20 18:28:27 CET 2024
Ending 14221857 at
Tue Feb 20 18:28:37 CET 2024
```

This job also created a text (.txt) file names 10.txt, let's take a look:

```bash
more 10.txt
```

```console
I slept for 10 seconds
```

### Debunging errors during sbatch execution:

Let's run the following script that should launch a job to sleep for 20 seconds and then create a text file ```20.txt```

-Copy the script to the ```$SCRATCH```

```bash

rsync -aPL $COURSES/BIO326/BestPracticesOrion/mysecondsbatch.SLURM.sh $SCRATCH/SLURMTest/

```

And then submit the job:

```bash
sbatch mysecondsbatch.SLURM.sh
```

```console
Submitted batch job 14225128
```

Let's list the results:

```bash
ls -lrth
```

```console
total 3.0K
-rw-rw-r-- 1 bio326-2024-1 nobody 1.2K Feb 20 18:11 myfirstsbatch.SLURM.sh
-rw-rw-r-- 1 bio326-2024-1 nobody    0 Feb 22 11:05 slurm-MySbatchScript_14225084.err
-rw-rw-r-- 1 bio326-2024-1 nobody   23 Feb 22 11:05 10.txt
-rw-rw-r-- 1 bio326-2024-1 nobody  195 Feb 22 11:05 slurm-MySbatchScript_14225084.out
-rw-rw-r-- 1 bio326-2024-1 nobody 1016 Feb 22 11:22 mysecondsbatch.SLURM.sh
-rw-rw-r-- 1 bio326-2024-1 nobody    0 Feb 22 11:23 20.txt
-rw-rw-r-- 1 bio326-2024-1 nobody   74 Feb 22 11:23 slurm-MySbatchScript_14225128.err
-rw-rw-r-- 1 bio326-2024-1 nobody  194 Feb 22 11:23 slurm-MySbatchScript_14225128.out

```

We can notice the ```20.txt``` file is empty (value 0 in the 4th colum), we can check then the ```slurm-MySbatchScript_14225128.err``` to debug:

```bash
more slurm-MySbatchScript_14225128.err
```

```console
/var/tmp/slurmd/job14225128/slurm_script: line 32: ech: command not found
```

It looks like the error is in line 32, let's take a look of the slurm script:

```bash
less +32 -N mysecondsbatch.SLURM.sh
```

```console
32 sleep 20 && ech "I sleept for 20 seconds" > 20.txt
```

Changing the error in that line will correct the code. After changing, save the script, look for the line to be corrected and if it is correct resubmit as follow:

```sbatch 
sed -i.back '32s/ech/echo/' mysecondsbatch.SLURM.sh
less +32 -N mysecondsbatch.SLURM.sh
```
```console
32 sleep 20 && echo "I sleept for 20 seconds" > 20.txt
```
Now the command ```echo``` is well spelled let's resubmit:

```bash
sbatch mysecondsbatch.SLURM.sh
```

```console
Submitted batch job 14225840
```
Now the ```20.txt``` file is written and has info on it:

```bash
less 20.txt
```

```console
I sleept for 20 seconds
20.txt (END)
```


### Loading CONDA in a sbatch script:

Similar to using an ```interactive job``` users can load ```modules``` and ```conda``` in the sbatch script. Let's take as an example the followign simple sbatch template:

```bash
#!/bin/bash

#####################SLURM flags
## Job name:
#SBATCH --job-name=LoadingCONDA
## Wall time limit:
#SBATCH --time=00:10:00
## Other parameters:
#SBATCH --cpus-per-task 2  #Asking for 2 CPUS
#SBATCH --mem=4G     #4G of RAM
#SBATCH --nodes 1 #Only one node 
#SBATCH --partition=smallmem  #partition to run
#SBATCH --output=slurm-%x_%A.put  #slurm standar output file
#########################################

######Everything below this are the job instructions######

#Loading modules
module purge #This remove any module loaded 
module load Miniconda3 && eval "$(conda shell.bash hook)"  ###This line loads the conda module and configures the computer to run conda envrionments

##Launching BLAST as a conda environment and display 

##Activate the environment
conda activate $COURSES/BIO326/BestPracticesOrion/BLASTConda

###Displaying the help:
echo "This is the BLASTp basic help..."
blastp --help


```

*You can copy and launch this script from ```$COURSES/BIO326/BestPracticesOrion/``` as follow:

```bash
rsync -aPL $COURSES/BIO326/BestPracticesOrion/LoadingCONDA.sbatch.SLURM.sh $SCRATCH/SLURMTest/

```

```console
sending incremental file list
LoadingCONDA.sbatch.SLURM.sh
            910 100%    0.00kB/s    0:00:00 (xfr#1, to-chk=0/1)
```

Then launch the job:

```bash
sbatch LoadingCONDA.sbatch.SLURM.sh

```

```console
Submitted batch job 14225870
```

Let's check the output

```bash
more slurm-LoadingCONDA_14225870.out
```
```console
This is the BLASTp basic help...
USAGE
  blastp [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-ipglist filename]
    [-negative_ipglist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-qcov_hsp_perc float_value] [-max_hsps int_value]
    [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value] [-seg SEG_options]
```
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


## Bulletpoints

* Do not use the login node to run process (e.g. BLAST, SPADES, HMMER).
* Do not use the $HOME partition for lagre files storage.
* Use interactive jobs for testing and debugging.
* Use the $TMPDIR for faster computation.
* Monitoring your jobs by squeue.
* Delete intermediate results from the $TMPDIR.
* Use sbatch command to submit your "final" jobs scripts.

## Enjoy the Orion Cluster...
