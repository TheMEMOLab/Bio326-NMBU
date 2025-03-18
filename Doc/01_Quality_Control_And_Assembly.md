# Recovering MAGs from Metagenomic data

<img src="https://github.com/TheMEMOLab/Bin420-Bioinformatics-for-Functional-Meta-Omics/blob/main/img/Assemlby.webp" height="400">

## 1. Quality control of metagenomic raw reads.

In this part of the BIO326 course we will work with metagemonic data. Usually, metagenomes are very large experiments and contains multiple individuals, environments or time-points. As this is a small course we will use a single timepoint of a metagenomic experiment from the [SupaCow](https://www.nmbu.no/en/research/projects/supacow) project from the [MEMO group](https://www.nmbu.no/en/research/groups/memo-group-microbial-ecology-and-meta-omics) at NMBU.


Let's create a set of directories to work on the HPC:

```bash
cd /cluster/projects/nn9987k/$USER
mkdir metaG &&  mkdir -p metaG/data && mkdir -p metaG/results 
/cluster/projects/nn9987k/.share/conda_environments/BASICS/bin/tree
```

```console
.
└── metaG
    ├── data
    └── results

4 directories, 0 files
```

Create a softlink of the data:

```bash
ln -s /cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/D01T6_T.fq.gz /cluster/projects/nn9987k/$USER/metaG/data/
/cluster/projects/nn9987k/.share/conda_environments/BASICS/bin/tree
```

```console
.
└── metaG
    ├── data
    │   └── D01T6_T.fq.gz -> /cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/D01T6_T.fq.gz
    └── results

4 directories, 1 file
```


### Cleanning reads: Running Chopper and NanoPlot.

First, as we want to have some results let's work interactivelly with the computer.
>[!Warning]
> Remember do not get naked in the lobby. The following script will help us to ask for a computer:

```bash
/cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/srun.prarameters.Nonode.Account.sh 8 20G normal 120G nn9987k 02:00:00
```

Let's review a bit on sequencing data concepts:

The FastQ files and the Phred score:

![FQ](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/fastqC.png)
![ASCII](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/ASCII.png)


Now that we have all the data we can use a combination of tools to perform a QC and cleanning the reads.

### Running Chopper and NanoPlot.

Now that we have all the data we can use a combination of tools to perform a QC and cleanning the reads.

Use the [Chopper](https://github.com/wdecoster/chopper) tool to clean the reads and 
[Nanoplot](https://github.com/wdecoster/NanoPlot) to visualize the QC stats. 

- Load the toolbox (Conda environment)

```bash 
module load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/NANOPAKQC/
echo "I am working with this" $CONDA_PREFIX

```

- We are now in a room (computer node), with the toolbox, Now let's bring the luggage (data) to our closset in the room.

```bash
tree="/cluster/projects/nn9987k/.share/conda_environments/BASICS/bin/tree"
cd $LOCALSCRATCH
rsync -aPLhv /cluster/projects/nn9987k/$USER/metaG/data/D01T6_T.fq.gz .
$tree
```
You will end with someting like:

```console
.
└── D01T6_T.fq.gz

1 directory, 1 file
```

Now let's display the first 4 lines of the file:

```bash
zcat D01T6_T.fq.gz |head -4|less -S
```

Something like this should be displayed:

```console
@38d36892-c536-4918-9c33-d0c3e207ce07 runid=b750b83920e5186c0fefcc3ad3d1ff830aec4a68 sampleid=CP958b read=4675 ch=2199 start_time=2022-07-30T08:02:56Z model_version_id=dna_r10.4.1_e8.2_su>CTGGCCAACAGATTCATAGCACAAGCCATTCCAGCATAAGTGGAGAATGGTTTCTTTTTGGCATCAGATGGCGTCTTTTCTCTGTAGGATGTCTTTTGATGTCGCTTTCCCTATACGTCCTGAAAAGCCACGAAATTTCACCTATGTTTCACCTGATGTCTTTTTGGTTTCTTATAAATCACTGATT>+
HFHGIGLJJMJJHVIJPDEEDFKLKLML{{QKKOPJJIPKRML{SNKHJ{HGLFFGFFAAKHNIJIPKKJXPA@@@=:I=<6677++******-----10.-,,.0@=78C=DCBA?=1100577>A??>>FG{KR{P?==6434489;66878324110/17999:>NMDBBABM{N{YLFIP{J
(END) 
```
>[!Tip]
>To exit press the letter Q in your keyboard