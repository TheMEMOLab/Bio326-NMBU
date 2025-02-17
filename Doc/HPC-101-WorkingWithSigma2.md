# Wroking with SAGA Cluster Sigma2-NIRS

This workshop was doing by help of [Sigma2-NIRS](https://documentation.sigma2.no/index.html).

### What is this?

This document is intended to be a quick reference guide on the basic usage of the FRAM Sigma2 HPC cluster. For a complete reference please referer to the full documentation of [SAGA](https://documentation.sigma2.no/hpc_machines/saga.html).

## Login into SAGA

For login open a Command-line interfase (CLI) or Terminal  and type something like this. 

```bash
ssh auve@saga.sigma2.no
```
> [!Important]
> Rember to change your user name

This will ask for your password, the one you got from Sigma2 by text. Type it
*Even you don't see anything the password is typed*

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

**The supercomputer is named after the goddess in norse mythology associated with wisdom. Saga is also a term for the Icelandic epic prose literature. The supercomputer, placed at NTNU in Trondheim is designed to run both sequential and parallel workloads. It was made available to users right before the start of the 2019.2 period (October 2019).

Saga is provided by Hewlett Packard Enterprise and has a computational capacity of approximately 140 million CPU hours a year.**

Let's take a look into this figure: 

![Cluster](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/cluster.png)


**All users have access to the $HOME, so please DO NOT USE THE $HOME FOR STORAGE LARGE FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and SLURM scripts**

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
> **NEVER RUN A JOB IN THE LOGIN NODE!!! THE LOGIN NODE IS ONLY FOR LOOKING AND MANAGING FILES, INSTALL SOFTWARE AND WRITE SCRIPTS** 

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

### SAGA Resources and information through SLURM.

If we want to know the amount of CPU, RAM and other configuration in the cluster, we can use a set of tools (commands) SLURM provide to find available resources in SAGA.

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



