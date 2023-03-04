# Software managment in Orion using CONDA

### This document provides a quick reference on how to use Conda environments to install software in Orion Cluster.

Conda is an open source package management system and environment management system that runs on Windows, macOS, and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language. Please refere to [CONDA](https://docs.conda.io/en/latest/) website for a complete documentation. 

### Using conda in Orion

**The following document uses $SCRATCH to create a "ToolBox" directory and then use mamaba to install a couple of bioinformatic software (minimap2 and samtools). However, this can be replicated in whaterver directory in the filesystem in Orion**

* Ask for a computer node to start working (remember "Don't get naked in the lobby"). The resources we are asking for are 2 CPUs (-c) 100 MB of RAM (--mem) 1 hr (-t) and using the smallmem partition (-p).

```console
qlogin -c 2 --mem=100M -t 01:00:00 -p smallmem
```

* Go to the ```$SCRATCH``` and create a directory named ToolBox:

```console
cd $SCRATCH && mkdir ToolBox
```

* Enter to the ToolBox directory:

```console
cd ToolBox/
pwd
/mnt/SCRATCH/bio326-2023-19/ToolBox
```

* We can use then [MAMBA](https://github.com/mamba-org/mamba) is a faster (c++) reimplementation of the conda package. For this we need to load the ```Miniconda3``` module and 
configure our ```shell``` (Terminal) to use conda using the following comands:

```bash

module load Miniconda3 && eval "$(conda shell.bash hook)"

``` 

After doing this our promt should change from ```bash-4.2$``` to ```(base) bash-4.2$``` This meand CONDA/MAMBA are ready to use. Then we can ```create``` an environment named for exmaples ```EUKVariantDetection``` and install on it minimap2, samtools  and sniffles. All these tools will be used for variant detection. All these tools are available at the BIOCONDA repository so we need to indicate this to mamba using ```-c bioconda``` flag:

```bash
mamba create -y --prefix ./EUKVariantDetection -c bioconda minimap2 samtools sniffles
```

If you type all correctly you will see something like this:

```
                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (0.27.0) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['minimap2', 'samtools', 'sniffles']

pkgs/r/linux-64                                               No change
pkgs/r/noarch                                                 No change
pkgs/main/noarch                                   819.6kB @   1.4MB/s  0.6s
bioconda/linux-64                                    4.6MB @   3.3MB/s  1.4s
pkgs/main/linux-64                                   5.2MB @   3.0MB/s  1.8s
bioconda/noarch                                      4.2MB @   1.9MB/s  2.3s
```

After the installation it will display something like:

```

To activate this environment, use

     $ mamba activate /net/fs-2/scale/OrionStore/Scratch/bio326-2023-19/ToolBox/EUKVariantDetection

To deactivate an active environment, use

     $ mamba deactivate
```

This means all the software were succesfully installed :-).

* Activate the environment and test if the software is installed by displaying the help/version:

```bash
conda activate /net/fs-2/scale/OrionStore/Scratch/bio326-2023-19/ToolBox/EUKVariantDetection
```
Here again our promt should change to:

```
(/net/fs-2/scale/OrionStore/Scratch/bio326-2023-19/ToolBox/EUKVariantDetection) bash-4.2$
```
* Ask for the minimap2, samtools and sniffles versions:

```bash
minimap2 -V
2.22-r1101
```
```bash
samtools --version
samtools 1.13
Using htslib 1.13
Copyright (C) 2021 Genome Research Ltd.
```
```niffles -V
Error:
         Couldn't find match for argument

Short usage:
       sniffles [options] -m <sorted.bam> -v <output.vcf>
Version: 1.0.12
```

* All the software is installed, now we can finish the ```qlogin``` session and use the conda environment in out SBATCH script.

```bash
exit
```

The following is the modified script that Marie Satiou wrote but using conda instead of singularity:

```bash
#!/bin/bash
#SBATCH --job-name=MappingVariants  # sensible name for the job
#SBATCH --mem=12G 
#SBATCH --ntasks=1   
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATHC -o slurm-%x-%A.out


##Activate conda environment

module load Miniconda3 && eval "$(conda shell.bash hook)"

### NB! Remember to use your own conda environment:

conda activate $SCRATCH/ToolBox/EUKVariantDetection 

GENOMEDIR='/mnt/courses/BIO326/EUK/bull_analysis/demo_data'

minimap2 \
    -t 8 \
    -a $GENOMEDIR/Bos_taurus.fa.gz \
    cleaned.bull.fastq.gz > bull.sam

# convert the sam file to bam format
samtools view -@ $SLURM_CPUS_ON_NODE -S -b bull.sam > bull0.bam

## sort the bam file
samtools sort -@ $SLURM_CPUS_ON_NODE bull0.bam -o bull.bam

# index the bam file
samtools index -@ $SLURM_CPUS_ON_NODE -M  bull.bam

# Variant Calling using Sniffles
sniffles -t $SLURM_CPUS_ON_NODE --input  bull.bam --vcf bull.vcf

##Finish


```

### You can use this script to perform your varian calling.

Enjoy :-)!





