# Wroking with SAGA Cluster Sigma2-NIRS

This workshop was doing by help of [Sigma2-NIRS](https://documentation.sigma2.no/index.html).

### What is this?

This document is intended to be a quick reference guide on the basic usage of the FRAM Sigma2 HPC cluster. For a complete reference please referer to the full documentation of [SAGA](https://documentation.sigma2.no/hpc_machines/saga.html).

## Login into FRAM

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

>[!Warning]
> **NEVER RUN A JOB IN THE LOGIN NODE!!! THE LOGIN NODE IS ONLY FOR LOOKING AND MANAGING FILES, INSTALL SOFTWARE AND WRITE SCRIPTS** 

**All users have access to the $HOME, so please DO NOT USE THE $HOME FOR STORAGE LARGE FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and SLURM scripts**

### Where can I storage large files? 

>[!Important]
> During the BIN420 Course all data, scripts, results and so must be written in
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