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

Let's review a bit on sequencing data concepts:

The FastQ files and the Phred score:

![FQ](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/fastqC.png)



Now that we have all the data we can use a combination of tools to perform a QC and cleanning the reads.

