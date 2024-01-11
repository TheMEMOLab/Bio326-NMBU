# Using DRAM annotation to extract certain features.


Although DRAM provides a visualization in the html report. Sometimes users want to extract certain information (e.g carbohydrate metabolism) from the annotation and combine it with for example taxonomy annotation.

# Study case: Extracting and ploting the abundance of CAZy genes of microbial population after DRAM annotation.

For this study case we can use the MAGs produced by the study: [Nitrous oxide respiring bacteria in biogas
digestates for reduced agricultural emissions. ISME J 16(2):580–590.](https://doi.org/10.1038/s41396-021-01101-x)

Let's copy the DRAM annotations from this to our ```$SCRATCH``` area:

```console
mkdir $SCRATCH/prok/results/PARSINGDRAM
cd $SCRATCH/prok/results/PARSINGDRAM
cp -r $COURSES/BIO326/PROK/PARSINGDRAM/dram.genome_summaries.dir .
```
Let's check the data:

```console
cd dram.genome_summaries.dir/
ls
genome_stats.tsv  metabolism_summary.xlsx  product.html  product.tsv
```
Then copy the GTDBTK annotation and the CHECKM quality of the MAGs here: 

```console
cp /mnt/courses/BIO326/PROK/PARSINGDRAM/*.tsv .
ls
genome_stats.tsv  gtdbtk.ar53.summary.tsv  gtdbtk.bac120.summary.tsv  metabolism_summary.xlsx  product.html  product.tsv
```

Now we can use R to combine this and extract the information.

### The following R script does the trick:

https://github.com/TheMEMOLab/MetaGVisualToolBox/blob/main/scripts/CAZYheatmap.R

Let's work on Orion, to do this we should ask for a computing node:

```console
qlogin -c 4 -p smallmem
salloc: Pending job allocation 14111809
salloc: job 14111809 queued and waiting for resources
```

Then let's load the git Module and clone the repository:

```console
module load git
git clone https://github.com/TheMEMOLab/MetaGVisualToolBox.git
```

This script requires some libraries, we have previously installed this in a conda environment we can load this and run the script:

```console
module load Miniconda3 && conda activate $ORION/conda/CIGENE/R_env
Rscript MetaGVisualToolBox/scripts/CAZYheatmap.R metabolism_summary.xlsx quality_report.tsv gtdbtk.bac120.summary.tsv gtdbtk.ar53.summary.tsv CAZYHeatmap
```

After running we got the PDF ```CAZYHeatmap.pdf```

This PDF has the information of Taxonomy as columns and all the CAZy genes as rows.




