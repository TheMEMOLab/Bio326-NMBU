# Basic exercises for best practices in Orion
v3, February, 2023 Author: Arturo Vera-Ponce de Leon

### What is this?
This document is intended to be a quick reference guide on the basic usage of the NMBU-Orion HPC cluster. For a complete reference please referer to the full documentation of [Orion](https://orion.nmbu.no/)

**Login into orion** 

To login into Orion cluster we need two things:
- Establish a [VPN](https://nmbuhjelp.nmbu.no/tas/public/ssp/content/detail/knowledgeitem?unid=4c4870bb-f0fc-4c9f-b0f5-0c036951436f) connection (if not using a NMBU network)
- Use a secure-shell command [ssh](https://en.wikipedia.org/wiki/SSH_(Secure_Shell)) via the command line

For login open a Command-line interfase (CLI) or Terminal  and type something like this. 

```bash
$ ssh bio326-2023-19@login.orion.nmbu.no
```
*Remember to change to your username bio326-y-x*

This will ask for your password. Type it

*Even you don't see anything the password is typed*

After the first attempt of login you will get something like:

```bash
ssh bio326-2023-19@login.orion.nmbu.no
bio326-2023-19@login.orion.nmbu.no's password:
You are required to change your password immediately (root enforced)
Last login: Tue Feb 21 17:39:16 2023 from 10.42.25.161
WARNING: Your password has expired.
You must change your password now and login again!
Changing password for user bio326-2023-19.
Changing password for bio326-2023-19.
(current) UNIX password:
```

Type the password again:

```bash
New password:
Retype new password:
passwd: all authentication tokens updated successfully.
Connection to login.orion.nmbu.no closed.
```

If this was done correctly, you will be logout and need to relogin:

```bash
ssh bio326-2023-19@login.orion.nmbu.no
bio326-2023-19@login.orion.nmbu.no's password:
Last login: Tue Feb 21 17:49:52 2023 from 10.42.25.161

Welcome to the NMBU Orion compute cluster environment.

You are logged in to a machine that can be used to access your home directory,
edit your scripts, manage your files, and submit jobs to the cluster environment.
Do not run any jobs on this machine, as they might be automatically terminated.

IMPORTANT:
  - Orion introduction: https://orion.nmbu.no/
  - Don’t use more than 150 threads/CPUs concurrently for more than one day. Please limit concurrently running array jobs to 25. Refer to Orion wiki how to do that. Orion can handle small-scale projects. Need more CPU hours? Please consider
    applying for national infrastructure resources: https://www.sigma2.no/
  - Please, PLEASE do compress your fastq, vcf and other non-compressed files
    using i.e. pigz.

NEWS:
  - 2022-07-28: Orion has been moved to a new storage system with a new data management policy. Please check Orion https://orion.nmbu.no/en/StorageChanges. We are still working out many details.
    Please email us if you miss anything or notice any issues.
  - Linux commands cp, rsync, scp, and mv are allowed to run for only 5 minutes in the login server. You must submit your jobs or use the interactive qlogin command to run these commands for longer than 5 minutes. You can also use filemanager.orion.nmbu.no to transfer files to your PC or arken.nmbu.no.

Apply for project storage https://tinyurl.com/OrionStorageApply

For any Orion related enquiry: orion-support@nmbu.no
We are on Teams: https://bit.ly/orion-teams


Your quota:
Volume          Type    ID                      Used    Quota   Limit   Grace
Home            USR     bio326-2023-19          16K     200G    300G    none

Have a nice evening!
-bash-4.2$
```

**Now you are logged into the Orion login-node.**

### Orion main configuration 

**The Orion HPC system currently consists of 1680 processor cores with more than 12 terabytes of RAM and 1 petabyte of storage accessible on a 10/1 Gbit network via NFS. The Orion HPC's computation nodes utilize various processors with Centos Linux 7.9 as an operating environment.**

Let's take a look into this figure: 

![Cluster](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/cluster.png)

**NEVER RUN A JOB IN THE LOGIN NODE!!! THE LOGIN NODE IS ONLY FOR LOOKING AND MANAGING FILES, INSTALL SOFTWARE AND WRITE SCRIPTS** 

How can I be sure of the number of CPUs and RAM of this "login" computer node and other nodes?

* CPUS: Use the command nproc

```Bash
-bash-4.2$ nproc
8
```

* RAM: We need to look for the "Total memory". All this info is allocated in the meminfo file at /proc directory. So we can use the grep command to look for this into the file.

```bash
-bash-4.2$ grep MemTotal /proc/meminfo |awk '{print $1,$2/1000000 " GB"}'
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
-bash-4.2$ sinfo -l -N
Tue Feb 21 17:53:54 2023
NODELIST   NODES   PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON
cn-1           1      orion*       mixed 144    4:18:2 309453        0      1 cpu_xeon none
cn-1           1     hugemem       mixed 144    4:18:2 309453        0      1 cpu_xeon none
cn-1           1       RSYNC       mixed 144    4:18:2 309453        0      1 cpu_xeon none
cn-1           1 interactive       mixed 144    4:18:2 309453        0      1 cpu_xeon none
cn-2           1      orion*       mixed 80     40:1:2 103186        0      1 cpu_xeon none
cn-2           1     hugemem       mixed 80     40:1:2 103186        0      1 cpu_xeon none
cn-2           1    smallmem       mixed 80     40:1:2 103186        0      1 cpu_xeon none
cn-2           1 interactive       mixed 80     40:1:2 103186        0      1 cpu_xeon none
cn-3           1      orion*       mixed 80     40:1:2 967353        0      1 cpu_xeon none
cn-3           1     hugemem       mixed 80     40:1:2 967353        0      1 cpu_xeon none
cn-3           1    smallmem       mixed 80     40:1:2 967353        0      1 cpu_xeon none
cn-3           1 interactive       mixed 80     40:1:2 967353        0      1 cpu_xeon none
cn-4           1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-4           1    smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-4           1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-5           1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-5           1    smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-5           1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-6           1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-6           1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-7           1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-7           1    smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-7           1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-8           1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-8           1    smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-8           1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-9           1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-9           1    smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-9           1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-10          1      orion*       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-10          1    smallmem       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-10          1 interactive       mixed 32     32:1:1 193230        0      1 cpu_xeon none
cn-11          1      orion*       mixed 64     2:16:2 257687        0      1 cpu_xeon none
cn-11          1    smallmem       mixed 64     2:16:2 257687        0      1 cpu_xeon none
cn-11          1 interactive       mixed 64     2:16:2 257687        0      1 cpu_xeon none
cn-12          1      orion*       mixed 32      2:8:2 257738        0      1 cpu_xeon none
cn-12          1    smallmem       mixed 32      2:8:2 257738        0      1 cpu_xeon none
cn-12          1 interactive       mixed 32      2:8:2 257738        0      1 cpu_xeon none
cn-13          1       RSYNC        idle 12      2:6:1  31917        0      1 cpu_xeon none
cn-14          1     hugemem       mixed 256    2:64:2 205153        0      1 cpu_amd, none
cn-14          1      orion*       mixed 256    2:64:2 205153        0      1 cpu_amd, none
cn-14          1 interactive       mixed 256    2:64:2 205153        0      1 cpu_amd, none
cn-15          1     hugemem   allocated 32      2:8:2 257738        0      1 cpu_xeon none
cn-15          1      orion*   allocated 32      2:8:2 257738        0      1 cpu_xeon none
cn-15          1 interactive   allocated 32      2:8:2 257738        0      1 cpu_xeon none
cn-16          1     hugemem       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-16          1      orion*       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-16          1 interactive       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-17          1     hugemem       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-17          1      orion*       mixed 256    2:64:2 103176        0      1 cpu_amd, none
cn-17          1 interactive       mixed 256    2:64:2 103176        0      1 cpu_amd, none
gn-0           1         gpu     drained 64     2:16:2 257710        0      1 cpu_amd, Job step not running
gn-1           1         gpu       mixed 64     2:16:2 257710        0      1 cpu_amd, none
gn-2           1         gpu       mixed 64     2:16:2 257710        0      1 cpu_amd, none
gn-3           1         gpu     drained 64     2:16:2 257710        0      1 cpu_amd, Kill task failed
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
[bio326-21-0@login ~]$ df -h
Filesystem                         Size  Used Avail Use% Mounted on
devtmpfs                            16G     0   16G   0% /dev
tmpfs                               16G     0   16G   0% /dev/shm
tmpfs                               16G  1.2G   15G   8% /run
tmpfs                               16G     0   16G   0% /sys/fs/cgroup
/dev/mapper/centos_login--0-root    28G   20G  7.6G  73% /
/dev/sda2                         1014M  287M  728M  29% /boot
/dev/sdb                           100G  2.5G   98G   3% /work
/dev/sda1                          200M   12M  189M   6% /boot/efi
tmpfs                              5.0G   29M  5.0G   1% /tmp
fs-1:/                             973M   11M  963M   2% /net/fs-1
fs-1:/projects01                    95T   95T  686G 100% /net/fs-1/projects01
fs-1:/home01                        76T   62T   15T  81% /net/fs-1/home01
cn-1:/mnt/SCRATCH                   77T   76T  1.3T  99% /net/cn-1/mnt/SCRATCH
tmpfs                              3.2G     0  3.2G   0% /run/user/1035
fs-1:/Geno                          29T   21T  8.1T  72% /net/fs-1/Geno
tmpfs                              3.2G     0  3.2G   0% /run/user/1018
tmpfs                              3.2G     0  3.2G   0% /run/user/10023
tmpfs                              3.2G     0  3.2G   0% /run/user/1004
cn-13:/mnt/BACKUP                   71T   43T   29T  61% /net/cn-13/mnt/BACKUP
fs-1:/results01                     95T   87T  8.3T  92% /net/fs-1/results01
fs-1:/Transpose                     38T   38T   18G 100% /net/fs-1/Transpose
tmpfs                              3.2G     0  3.2G   0% /run/user/10305
tmpfs                              3.2G     0  3.2G   0% /run/user/10274
tmpfs                              3.2G     0  3.2G   0% /run/user/50043
tmpfs                              3.2G     0  3.2G   0% /run/user/50045
tmpfs                              3.2G     0  3.2G   0% /run/user/50102
tmpfs                              3.2G     0  3.2G   0% /run/user/50067
tmpfs                              3.2G     0  3.2G   0% /run/user/10197
tmpfs                              3.2G     0  3.2G   0% /run/user/30055
tmpfs                              3.2G     0  3.2G   0% /run/user/50095
tmpfs                              3.2G     0  3.2G   0% /run/user/40023
tmpfs                              3.2G     0  3.2G   0% /run/user/1034
cvmfs2                             4.9G  2.2G  2.7G  45% /cvmfs/singularity.galaxyproject.org
cn-13:/mnt/SCRATCH2                148T  147T  691G 100% /net/cn-13/mnt/SCRATCH2
fs-1:/SandveLab                     38T   36T  1.7T  96% /net/fs-1/SandveLab
fs-1:/PEPomics01                    34T   29T  4.8T  86% /net/fs-1/PEPomics01
tmpfs                              3.2G     0  3.2G   0% /run/user/50090
tmpfs                              3.2G     0  3.2G   0% /run/user/30047
tmpfs                              3.2G     0  3.2G   0% /run/user/10191
tmpfs                              3.2G     0  3.2G   0% /run/user/50120
tmpfs                              3.2G     0  3.2G   0% /run/user/10209
tmpfs                              3.2G     0  3.2G   0% /run/user/10069
tmpfs                              3.2G     0  3.2G   0% /run/user/1080
tmpfs                              3.2G     0  3.2G   0% /run/user/50131
tmpfs                              3.2G     0  3.2G   0% /run/user/40015
tmpfs                              3.2G     0  3.2G   0% /run/user/10200
tmpfs                              3.2G     0  3.2G   0% /run/user/1010
tmpfs                              3.2G     0  3.2G   0% /run/user/10181
fs-1:/HumGut                        19T   15T  4.8T  76% /net/fs-1/HumGut
fs-1:/results03                     48T   44T  4.0T  92% /net/fs-1/results03
cn-1:/mnt/SALMON-SEQDATA            37T   37T  313G 100% /net/cn-1/mnt/SALMON-SEQDATA
fs-1:/Ngoc                         2.0T  1.2T  759G  61% /net/fs-1/Ngoc
fs-1:/TestFile                     973G  2.7G  971G   1% /net/fs-1/TestFile
cn-1:/mnt/labdata01                 81T   70T   11T  87% /net/cn-1/mnt/labdata01
fs-1:/IPVProjects01                 57T   54T  3.7T  94% /net/fs-1/IPVProjects01
cn-1:/mnt/homearchive               37T   24T   13T  66% /net/cn-1/mnt/homearchive
fs-1:/Foreco                        12T  1.9T  9.6T  17% /net/fs-1/Foreco
fs-1:/PreventADALL                 4.8T  2.4T  2.4T  50% /net/fs-1/PreventADALL
fs-1:/Home_turhamar                973G  691G  283G  71% /net/fs-1/Home_turhamar
fs-1:/Home_alme                    973G  396G  578G  41% /net/fs-1/Home_alme
fs-1:/Home_rush                    973G  167G  807G  18% /net/fs-1/Home_rush
//10.209.0.205/Completed_Projects  932G  583G  349G  63% /mnt/smb/GT2
//10.209.0.10/Completed_projects   932G  585G  348G  63% /mnt/smb/GT1
//10.209.0.203/Completed_Projects  932G  538G  395G  58% /mnt/smb/GT4
//10.209.0.204/Completed_Projects  932G  578G  355G  62% /mnt/smb/GT3
//10.209.0.202/Completed_Projects  932G  509G  424G  55% /mnt/smb/GT5
cvmfs2                             4.9G  2.2G  2.7G  45% /cvmfs/cvmfs-config.galaxyproject.org
tmpfs                              3.2G     0  3.2G   0% /run/user/1002
tmpfs                              3.2G     0  3.2G   0% /run/user/10028
tmpfs                              3.2G     0  3.2G   0% /run/user/10201
tmpfs                              3.2G     0  3.2G   0% /run/user/10289
tmpfs                              3.2G     0  3.2G   0% /run/user/1032
tmpfs                              3.2G     0  3.2G   0% /run/user/50157
tmpfs                              3.2G     0  3.2G   0% /run/user/1003
tmpfs                              3.2G     0  3.2G   0% /run/user/10205
tmpfs                              3.2G     0  3.2G   0% /run/user/50133
tmpfs                              3.2G     0  3.2G   0% /run/user/1027
tmpfs                              3.2G     0  3.2G   0% /run/user/4000
tmpfs                              3.2G     0  3.2G   0% /run/user/50094
```

As you can notice there are plenty of directories in Orion, but let's focus in the $HOME partition:

```bash
-bash-4.2$ df -h .
Filesystem      Size  Used Avail Use% Mounted on
fs-1:/home01     76T   62T   15T  81% /net/fs-1/home01
```
*syntaxis: d(isk)f(ree) -h(uman redable) .(the directory I am now) *

**All users have access to the $HOME, so please DO NOT USE THE $HOME FOR STORAGE LARGE FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and SLURM scripts**

### Where can I storage large files? 

There are two dierectories designed to this:
* $SCRATCH
* $PROJECT

As a student of this course, we are using the $SCRATCH folder to keep our raw sequencing files and final results. This directory, in contrast to the $HOME and $PROJECT, is not backed up. **Remember to make a copy of your important files into another (permanent) location!!!** All students have a directory in the $SCRATCH: /mnt/SCRATCH/bio326-22-x ; where x is the student number.

Let's move into that folder:

```bash
-bash-4.2$ cd $SCRATCH
-bash-4.2$ pwd
/mnt/SCRATCH/bio326-2023-19
```
*For moving here we use the variable ```$SCRATCH```. A variable is a character string to which the user assign a value. The value assigned could be a number, text, filename, device, or any other type of data. It is important you get familiarized with this and other Unix/Linux terminology. This [link](https://orion.nmbu.no/en/LinuxTutorial) can help you to navigate through these terms.

### Queues for different RAM usage:

Orion offers three primary job queues, or partitions: hugemem, smallmem and gpu. You are encouraged to choose the most appropriate partition.

|Partition 	|Job memory requirements| 	Available resources|
|-----------|-----------------------|----------------------|
|hugemem 	|> 100 GB RAM 	|2 x (1 TB RAM, 80 cores)|
|smallmem |	1-100 GB RAM 	|7 x (192 GB RAM, 32 cores)|
|gpu 	|1-100 GB RAM 	|4 x (256 GB RAM, 64 cores, 3 Quadro RTX 8000)|

To select the correct partition queue we can use the **SLURM** command: *--partition=\<partition>* and choosing the appropriate partition depending the amount of memory RAM the job will use to compute.
In general, there is a shorter job queue in smallmem than the hugemem queue. Also, be prepared to wait for your job to start if you need to use the hugemem queue.

## Running an interactive job to test programs and get used to working in the cluster

The easiest way to test software and look into huge files without messing the login node and other users, is by running an **interactive** job in Orion. This means you can "book" a compute node and type your commands directly in that node. Let's run an interactive job by the following commands:

```bash
-bash-4.2$ srun --partition=smallmem --cpus-per-task 2 --mem=2G --time=01:00:00 --pty bash -i
```

*Basic syntaxis the command:
 srun \<slurm-options> \<software-name/path>*
  
It might take a while to SLURM allocate the resources of this job. But as soon as it allocates the job a message like this will be displayed:

```bash
bash-4.2$
```
*NB! At the moment writing this documentation we found some issues with the ```$PS1``` variable displaying the prompt information to solve this:*

```bash
bash-4.2$ PS1="[\u@\h \W]$ "
[bio326-2023-19@cn-8 ~]$
```

You can notice that now the prompt has changed and shows the node (computer) we are running on. In this case the node "cn-x". Also if this is not displayed we can take advantage of the many [SLURM_environment_variables](https://slurm.schedmd.com/pdfs/summary.pdf). These are dynamic values that SLURM uses to control the computers. For example, if you would like to know what is the node I am working on and the number of CPUs requested for this job you can print the values of that by using different SLURM variables and the command "echo" follows by the name of the variable:

```bash
[bio326-2023-19@cn-8 ~]$ echo $SLURM_NODELIST
cn-8
[bio326-2023-19@cn-8 ~]$ echo $SLURM_CPUS_ON_NODE
2
```

Here we can run short parsing scripts, test software with a small datasets, etc. 

### Temporary working directory, faster and more efficient Jobs

Generaly any software can read (data) and write (results) from any partition of the cluster (i.e. $HOME, $SCRATCH, $PROJECT), however, I/O (reading and writing) from those locations uses a lot of network-trafic resources resulting in a high inefficenfy for heavy jobs (e.g mapping reads to large genomes/metagenomes or assembly genomes). Also if multiple users are running jobs in the same way the traffic in the network, even using the infiniband, makes the jobs super slow. 
To avoid this we can take advantage of the **$TMPDIR** partition. This is a physical hard-drive allocated in each of the computer nodes. We can migrate the data to there for faster I/O. Often, quite some efficiency can be gained by doing this.

Let's take a look, first we need to check if our **$USER** exists in that **$TMPDIR**

```
[bio326-2023-19@cn-8 ~]$ echo $TMPDIR/$USER
/home/work/bio326-2023-19
```

This means the user **bio326-y-u** has a directory in the **$TMPDIR** (/home/work). Move to that directory:

```bash
[bio326-2023-19@cn-8 ~]$ cd $TMPDIR/$USER
bash: cd: /home/work/bio326-2023-19: No such file or directory

```

It could happen that our $USER name is not in the $TMPDIR, if that is the case we can easily create this directory and go there:

```bash
[bio326-2023-19@cn-8 ~]$ mkdir $TMPDIR/$USER
[bio326-2023-19@cn-8 ~]$ cd $TMPDIR/$USER
[bio326-2023-19@cn-8 bio326-2023-19]$ pwd
/home/work/bio326-2023-19
```

Now we need to create an other directory **a work directory** to copy data for executing some commands. We can use another SLURM variable, let's say the JOBID to be consistent.


```bash
[bio326-2023-19@cn-8 bio326-2023-19]$ mkdir work.dir.of.$SLURM_JOB_ID
[bio326-2023-19@cn-8 bio326-2023-19]$ ls
work.dir.of.11818800
[bio326-2023-19@cn-8 work.dir.of.11818800]$ cd $TMPDIR/$USER/work.dir.of.$SLURM_JOB_ID
[bio326-2023-19@cn-8 work.dir.of.11818800]$ pwd
/home/work/bio326-2023-19/work.dir.of.11818800
```

By using the $SLURM_JOB_ID we can further identify what job we are running.

Let's enter to that directory and then copy some fasta files from the ```$COURSES/BIO326```, this is a share directory we (teachers) will use to upload data for you. In this class we are using the files from ```$COURSES/BIO326/BestPracticesOrion/BLASTExample``` path.

First take a look of the data

```bash
[bio326-2023-19@cn-8 bio326-2023-19]$ ls -l $COURSES/BIO326/BestPracticesOrion/BLASTExample
total 2049
-rwxrwxr-x 1 auve bio326     838 Nov  7 23:16 amylase.Bgramini.fasta
-rwxrwxr-x 1 auve bio326 2085506 Nov  7 23:06 Bacteroides51.faa
```

*Tip: Having more than one terminal open does help to faster look into multiple directories*

As you can see there are multiple files here, lets copy the two fasta files **.faa and .fasta** into the $TMPDIR/$USER/work.dir.of.$SLURM_JOB_ID

```bash
[bio326-21-0@cn-3 work.dir.of.12314866]$ cp /mnt/courses/BIO326/BestPracticesOrion/BLASTExample/*.fa* .
[bio326-21-0@cn-3 work.dir.of.12314866]$ ls
amylase.Bgramini.fasta  Bacteroides51.faa
```

*Remember that you can copy multiple files using regular expression (REGEX) in this case* * *.fa* * *means "everything that has .fa on it"*

No we can do some work on this files. Take a look of the **amylase.Bgramini.fasta** file 

```bash
[bio326-21-0@cn-3 work.dir.of.12314866]$ more amylase.Bgramini.fasta 
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
[bio326-2023-19@cn-8 work.dir.of.11818800]$ blastp
bash: blastp: command not found
```

It seems blastp is not installed as a default software in Orion.

### Conda envrionment:

Conda is an open source package management system and environment management system that runs on Windows, macOS, and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language. You can read more about conda [here](https://docs.conda.io/en/latest/).

For **BIO-326** we will use different Conda environment previously installed in Orion. However, as a Orion user you are allowed to install your own environment, please refere to the [Orion-Conda Environment](https://orion.nmbu.no/en/CondaEnvironment) for doing this.

Let's check if conda works in the node

```bash
[bio326-2023-19@cn-12 ~]$ conda --version
conda 4.14.0
```

Now we can **activate** the conda environment:

```bash
[bio326-2023-19@cn-12 work.dir.of.11818800]$ conda activate $COURSES/BIO326/BestPracticesOrion/BLASTConda

CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
To initialize your shell, run

    $ conda init <SHELL_NAME>

Currently supported shells are:
  - bash
  - fish
  - tcsh
  - xonsh
  - zsh
  - powershell

See 'conda init --help' for more information and options.

IMPORTANT: You may need to close and restart your shell after running 'conda init'.
```

If this happen we need to configure the Node to be able to work with conda. The following commands help to do that:

```bash
[bio326-2023-19@cn-12 work.dir.of.11818800]$ module load Miniconda3 && eval "$(conda shell.bash hook)"

(base) [bio326-2023-19@cn-12 work.dir.of.11818800]$

```
What we are doing here is being sure Conda is loaded and then export all the conda configurations to our shell...**NB! Remember the aim of this course is not to be a Linux expert so do not worry if this is a bit criptic for you :-)** 

Now conda is activate, we can see how our prompt has changed: ``` (base) [bio326-2023-19@cn-12 work.dir.of.11818800]$```

Then we can activate the environmet for run conda in this case this is at ```$COURSES/BIO326/BestPracticesOrion/BLASTConda```

```bash
(base) [bio326-2023-19@cn-12 work.dir.of.11818800]$ conda activate $COURSES/BIO326/BestPracticesOrion/BLASTConda
(/mnt/courses/BIO326/BestPracticesOrion/BLASTConda) [bio326-2023-19@cn-12 work.dir.of.11818800]$
```

We can then run a blast experiment:

* Index a database using ```makeblastdb``` and the molecule type prot:

```bash
makeblastdb -dbtype prot -in Bacteroides51.faa


Building a new DB, current time: 02/21/2023 19:54:33
New DB name:   /home/work/bio326-2023-19/work.dir.of.11818800/Bacteroides51.faa
New DB title:  Bacteroides51.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 4630 sequences in 0.182254 seconds.
```

Check the results by ls:

```bash

ls
amylase.Bgramini.fasta  Bacteroides51.faa.pdb  Bacteroides51.faa.pin  Bacteroides51.faa.pot  Bacteroides51.faa.ptf
Bacteroides51.faa       Bacteroides51.faa.phr  Bacteroides51.faa.pjs  Bacteroides51.faa.psq  Bacteroides51.faa.pto
```

And now lets run the BLAST,as we want to search for protein in a protein database the command we need to use is BLASTP:

```bash
blastp -query amylase.Bgramini.fasta -db Bacteroides51.faa -dbsize 1000000000 -max_target_seqs 1 -outfmt 6 -num_threads $SLURM_CPUS_ON_NODE -out amylase.Bgramini.fasta.blastp.out
Warning: [blastp] Examining 5 or more matches is recommended
```

Take a look into the results:

```bash
more amylase.Bgramini.fasta.blastp.out
WP_024997086.1  D0T87_RS12665   57.772  772     301     13      8       763     28      790     0.0     908
```

It seems the amylase of *B. fragilis* has a match wiht the D0T87_RS12665 sequence of Bacteroides51. We can corroborate this by looking into the fasta file annotation header by doing something like this:

```bash
grep D0T87_RS12665 Bacteroides51.faa
>D0T87_RS12665	alpha-amylase	WP_163175496.1
```
We found the amylase!!!

### Copy results to the $SCRATCH, remove work.directory and exit the job.

**NB! Remember the $TMPDIR is a temprary directory so, we need to move the results back to our $SCRATCH partition. For this we can use the following sintax:**

```bash
cp *fasta.blastp.out $SCRATCH
```

Then let's sure this is copy back

```bash
ls $SCRATCH
amylase.Bgramini.fasta.blastp.out
```

Finally, as the **$TMPDIR** is used for everyone a best practice is to delete all the temporary directories (i.e work.directory) from this location.

We can achive this by doing this:

* First go back to the main $TMPDI/$USER

```
[bio326-2-0@cn-3 work.dir.of.12314866]$ cd $TMPDIR/$USER
bio326-21-0@cn-3 bio326-21-0]$ ls
singularity  work.dir.of.12314866
```

Now we need to remove the work.dir.of 

```
[bio326-21-0@cn-3 bio326-21-0]$ rm -rf work.dir.of.*
[bio326-21-0@cn-3 bio326-21-0]$ ls
```

Finally, we can logout of this node:

```
exit
[bio326-2023-19@login ~]$
```

You can see now we return to the main **login bio326-y-x** node.

## Submit the same BLAST job but using a SLURM script.

Most of the time you do not use the interactive way for submiting jobs into the cluster. To submit jobs, you need to write all the instructions you want the computer to execute. This is what a script is.

SLURM uses a [bash](https://www.gnu.org/software/bash/) (computer language) base script to read the instructions. The first lines, are reserved words that SLURM needs to read inorder to launch the program:

```
-p --partition <partition-name>       --pty <software-name/path>
--mem <memory>                        --gres <general-resources>
-n --ntasks <number of tasks>         -t --time <days-hours:minutes>
-N --nodes <number-of-nodes>          -A --account <account>
-c --cpus-per-task <number-of-cpus>   -L --licenses <license>
-w --nodelist <list-of-node-names>    -J --job-name <jobname>
```

We can indicate these options by using the ```#SBATCH``` word following by any of these flag (e.g -c 10 ; means 10 CPUs).


```
#!/bin/bash

## Job name:
#SBATCH --job-name=Blast
#
## Wall time limit:
#SBATCH --time=00:00:00
#
## Other parameters:
#SBATCH --cpus-per-task 12
#SBATCH --mem=60G
#SBATCH --nodes 1
```

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
#SBATCH --output=slurm-%x_%A  #slurm standar output file
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


## Monitoring the jobs by squeue

A user can monitorate the status of the Job by the command ```squeue $USER``` 

```bash
squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          11818971  smallmem MyFirstB bio326-2  R       0:01      1 cn-12
```

This will show all the jobs, quite difficult to read. So instead we can indicate only to show our user jobs by adding the flag **-u** and the **$USER** variable:

When the job starts it produces an out file **slurm-JOBNAME-$JOB_ID.out**:

```bash
ls
jupyterhubLog  myfirstblast.SLURM.sh  slurm-MyFirstBlastp_11818971
```

We can check into this file:


```bash
more slurm-MyFirstBlastp_11818971
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

```
ls /mnt/SCRATCH/bio326-2023-19/MyFirstBlastp.dir 
amylase.Bgramini.fasta.blastp.out
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

## Bulletpoints

* Do not use the login node to run process (e.g. BLAST, SPADES, HMMER).
* Do not use the $HOME partition for lagre files storage.
* Use interactive jobs for testing and debugging.
* Use the $TMPDIR for faster computation.
* Monitoring your jobs by squeue.
* Delete intermediate results from the $TMPDIR.
* Use sbatch command to submit your "final" jobs scripts.

## Enjoy the Orion Cluster...
