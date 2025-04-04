# Case-Study: Working with the BIO326 metagenomic samples.

## The effect of DNA extraction methods on ONT sequence quality.

DNA extraction in BIO326 for metagenomics 

![METAG](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/bio326metag.JPG)


In the wetlab samples were divided as follow:

![SAMPLES](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/DNAsamples.JPG)

These samples were sequenced using the PromethION and it produced something like this:

```bash
(base) prom@PC24B170:/data/20250319_Bio326_PROK/20250319_Bio326_PROK/20250319_1452_3C_PAY86999_d557822e$ ls fastq_pass/
barcode01  barcode03  barcode05  barcode07  barcode09  barcode11  barcode13  barcode15  barcode17  barcode19  barcode21  barcode23
barcode02  barcode04  barcode06  barcode08  barcode10  barcode12  barcode14  barcode16  barcode18  barcode20  barcode22  unclassified
```

As you see there are multiple folders with different barcodes and inside a lot of fastqfiles:

```bash
(base) prom@PC24B170:/data/20250319_Bio326_PROK/20250319_Bio326_PROK/20250319_1452_3C_PAY86999_d557822e/fastq_pass/barcode01$ ls |head -3
PAY86999_pass_barcode01_d557822e_20eb73bc_0.fastq.gz
PAY86999_pass_barcode01_d557822e_20eb73bc_100.fastq.gz
PAY86999_pass_barcode01_d557822e_20eb73bc_101.fastq.gz
```

We can use the information from the Lysis method, Group and Barcode and merge these fastq files into single fastq files for each lyisis method and Group using a loop as follow, having a text file like:

```bash
less barcodes.tsv

Vortex_SRE_1    barcode01
FastPrep_1      barcode02
Vortex_1        barcode03
Vortex_2        barcode04
FastPrep_2      barcode05
Vortex_SRE_2    barcode06
Vortex_3        barcode07
Vortex_SRE_3    barcode08
FastPrep_3      barcode09
Vortex_3        barcode10
FastPrep_4      barcode11
Vortex_4        barcode12
Vortex_SRE_4    barcode13

```

and use a loop:

```bash
 cat barcodes.tsv |while read -r line; do NAME=$(echo $line|awk '{print $1}'); BC=$(echo $line|awk '{print $2}'); echo $BC; zcat fastq_pass/$BC/*fastq.gz |pigz -p 12 >  rawdata/$NAME.fastq.gz; done
```
Then we will end with something like:

```bash
ls -1 rawdata/
FastPrep_1.fastq.gz
FastPrep_2.fastq.gz
FastPrep_3.fastq.gz
FastPrep_4.fastq.gz
Vortex_1.fastq.gz
Vortex_2.fastq.gz
Vortex_3.fastq.gz
Vortex_4.fastq.gz
Vortex_SRE_1.fastq.gz
Vortex_SRE_2.fastq.gz
Vortex_SRE_3.fastq.gz
Vortex_SRE_4.fastq.gz
```

Then using the usefull NanoPlot:

```bash
 parallel -j 12 "NanoPlot --fastq {} --N50 --loglength -o {}.Nanoplot.dir" ::: *.gz
```

And as we have been working we can collect the NanoStats.txt

```
FastPrep_1.NanoStats.txt
FastPrep_2.NanoStats.txt
FastPrep_3.NanoStats.txt
FastPrep_4.NanoStats.txt
Vortex_1.NanoStats.txt
Vortex_2.NanoStats.txt
Vortex_3.NanoStats.txt
Vortex_4.NanoStats.txt
Vortex_SRE_1.NanoStats.txt
Vortex_SRE_2.NanoStats.txt
Vortex_SRE_3.NanoStats.txt
Vortex_SRE_4.NanoStats.txt
```

Let's work with these files. Download the files from 

https://arken.nmbu.no/~auve/BIO326/bio326NanoStats.Prok.dir.zip



>[!Tip]
> If you have access to the terminal you can download the file by:
> ``` wget https://arken.nmbu.no/~auve/BIO326/bio326NanoStats.Prok.dir.zip ```

Decompress the Zip file and let's move to R and RStudio.
>[!Tip]
> If you have access to the terminal you unzip the file by:
> ``` unzip bio326NanoStats.Prok.dir.zip ```

```bash
Archive:  bio326NanoStats.Prok.dir.zip
   creating: bio326NanoStats.Prok.dir/
  inflating: bio326NanoStats.Prok.dir/FastPrep_2.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/FastPrep_3.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/FastPrep_1.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_SRE_1.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_3.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_4.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_SRE_4.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_1.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/FastPrep_4.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_2.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_SRE_2.NanoStats.txt
  inflating: bio326NanoStats.Prok.dir/Vortex_SRE_3.NanoStats.txt
```



### Using RStiudio to load these files:

```R
library(tidyverse)

# Define the directory containing your files
data_dir <- "."  # Change if needed

# Get the list of NanoStat output files
files <- list.files(data_dir, pattern = "\\.NanoStats\\.txt$", full.names = TRUE)

# Function to read and process each NanoStat file
read_nanostat <- function(file) {
  lines <- read_lines(file)  # Read file line by line
  lines <- lines[lines != ""]  # Remove empty lines
  
  # Extract key-value pairs using regex
  stats <- lines %>%
    str_trim() %>%  # Trim whitespace
    str_replace_all(",", "") %>%  # Remove thousands separator
    str_match("^(.*?):\\s+([0-9\\.]+)$") %>%  # Match key-value pairs
    as_tibble() %>%
    filter(!is.na(V2)) %>%
    select(Metric = V2, Value = V3) %>%
    mutate(Value = as.numeric(Value))
  
  sample_name <- gsub("\\.NanoStats\\.txt$", "", basename(file))  # Extract sample name
  stats <- stats %>% mutate(Sample = sample_name)
  
  return(stats)
}

# Read all files and combine into one dataframe
stats_data <- map_df(files, read_nanostat)
```

Pivot for easy ploting

```R
stats_wide <- stats_data %>%
  pivot_wider(names_from = Metric, values_from = Value)
head(stats_wide,2)
```

Convert Sample column to factor to maintain order

```R
stats_wide <- stats_wide %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample)))
```

Let's create a function to plot all the metrics:

```R
library(RColorBrewer)
# Define metrics of interest (must match labels in the files!)
metrics <- c(
  "Number of reads",
  "Total bases",
  "Median read length",
  "Mean read length",
  "STDEV read length",
  "Read length N50",
  "Mean read quality",
  "Median read quality"
)
# Define 8 colors from the "Dark2" palette into a factor
metric_colors <- setNames(brewer.pal(8, "Dark2"), metrics)

# Create sorted horizontal bar plots with metric-based colors
plot_list <- map(metrics, function(metric) {
  ggplot(stats_wide, aes(y = fct_reorder(Sample, !!sym(metric)), x = !!sym(metric))) +
    geom_bar(stat = "identity", fill = metric_colors[metric]) +  # Assign color per metric
    theme_minimal() +
    labs(title = metric, x = metric, y = "Sample") +
    theme(axis.text.y = element_text(size = 10))   # Adjust text size if needed
}) %>%
  set_names(gsub(" ","", metrics))
```

Check the results, by displaying the Mean read length and the N50

```R
plot_list$Numberofreads
plot_list$ReadlengthN50
```

We can use Patchwork to plot in the same canvas:

```r
library(patchwork)
plot_list$Numberofreads + plot_list$ReadlengthN50
```

As we can see the Bar plot is not the best plot to display and compare these results...

>[!Tip]
> Let's use a violin plot:

```R
ReadLengthVP <- stats_wide %>% 
  select(Sample,`Mean read length`)%>%
  mutate(Method = case_when(
    str_detect(Sample, "FastPrep") ~ "FastPrep",
    str_detect(Sample, "Vortex_SRE") ~ "Vortex_SRE",
    str_detect(Sample, "Vortex") ~ "Vortex"
  )) %>% 
  ggplot(aes(x = Method, y = `Mean read length`, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8)+ 
  ggtitle('Mean Read Length')

N50VP <- stats_wide %>% 
  select(Sample,`Read length N50`)%>%
  mutate(Method = case_when(
    str_detect(Sample, "FastPrep") ~ "FastPrep",
    str_detect(Sample, "Vortex_SRE") ~ "Vortex_SRE",
    str_detect(Sample, "Vortex") ~ "Vortex"
  )) %>% 
  ggplot(aes(x = Method, y = `Read length N50`, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8)+ 
  ggtitle('Read length N50')

ReadLengthVP + N50VP
```

We can use some statistics like ANOVA and t-test to check if there is differences in the groups. The library ggpubr allows us to calculate these kind of statistics and visualize on the ggplot2 object:

```R
library(ggpubr)

ReadLengthVP <- ReadLengthVP +  # Global test
  stat_compare_means(
    comparisons = list(
      c("FastPrep", "Vortex"),
      c("FastPrep", "Vortex_SRE"),
      c("Vortex", "Vortex_SRE")
    ),
    method = "t.test",  # Technically not perfect post-ANOVA but okay for vis
    label = "p.signif"
  ) 
N50VP <- N50VP +
  stat_compare_means(
    comparisons = list(
      c("FastPrep", "Vortex"),
      c("FastPrep", "Vortex_SRE"),
      c("Vortex", "Vortex_SRE")
    ),
    method = "t.test",  # Technically not perfect post-ANOVA but okay for vis
    label = "p.signif"
  ) 

ReadLengthVP + N50VP
```

And modify some Color and astetics of the plot:

```R
ReadLengthVP <- ReadLengthVP +
  theme_minimal() +
  scale_fill_brewer(palette =  "Dark2") +
  theme(legend.position = "none")

N50VP <- N50VP +
  theme_minimal() +
  scale_fill_brewer(palette =  "Dark2") +
  theme(legend.position = "none")
  
ReadLengthVP + N50VP
```

Saving the plot:

```R
MergedQCPlot <- ReadLengthVP + N50VP

ggsave(MergeQCPlot,file="MergedPlot.QualityScores.pdf")
```

## Taxonomy classification using Kraken2

What is Kraken2? According to [Wood and Salzberg, GEnome Biology 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46): "*Kraken is an ultrafast and highly accurate program for assigning taxonomic labels to metagenomic DNA sequences. Previous programs designed for this task have been relatively slow and computationally expensive, forcing researchers to use faster abundance estimation programs, which only classify small subsets of metagenomic data. Using exact alignment of k-mers, Kraken achieves classification accuracy comparable to the fastest BLAST program.*"

![K2ALG](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/k2algo.JPG)

We can run Kraken2 something like this:

```bash
kraken2 --db shared_databases/kraken2/PlusPF-8/ --threads 10 --output FastPrep_1.kraken2.nonames.out --report FastPrep_1.kraken2.report.tsv FastPrep_1.fastq.gz
```
It will produce a report like this:

```bash
head -11 FastPrep_1.fastq.kraken2.report.tsv
 89.13  77005   77005   U       0       unclassified
 10.87  9389    15      R       1       root
 10.84  9369    148     R1      131567    cellular organisms
 10.22  8827    1623    D       2           Bacteria
  5.01  4329    157     K       1783272       Bacillati
  3.47  2994    92      P       1239            Bacillota
  1.28  1102    70      C       186801            Clostridia
  0.78  673     6       O       3085636             Lachnospirales
  0.77  664     58      F       186803                Lachnospiraceae
  0.32  279     0       G       3031933                 Chordicoccus
  0.32  279     279     S       2709410                   Chordicoccus furentiruminis
```

We can parse this report and only keep the organisms with > 0.1 % of reads classified:

```bash
head  FastPrep_1.kraken2.species.tsv
  0.32  2709410 Chordicoccus_furentiruminis
  0.69  2487118 Intestinibaculum_porci
  0.66  907     Megasphaera_elsdenii
  0.12  1064535 Megasphaera_elsdenii_DSM_20460
  0.20  1280    Staphylococcus_aureus
  0.54  2903756 Streptomyces_sp_NBC_01167
  0.17  2913620 Prevotella_sp_E2-28
  0.26  2913614 Prevotella_communis
  0.41  839     Xylanibacter_ruminicola
  0.15  264731  Xylanibacter_ruminicola_23
```

## Comparing changes in the microbial composition in the different treatments (FastPrep, Vortex, Vortex + SRE)

All the Kraken2 results can be donwloaded here:

https://arken.nmbu.no/~auve/BIO326/kraken2Reports.zip


>[!Tip]
> If you have access to the terminal you can download the file by:
> ``` wget https://arken.nmbu.no/~auve/BIO326/kraken2Reports.zip ```

Decompress the Zip file and let's move to R and RStudio.
>[!Tip]
> If you have access to the terminal you unzip the file by:
> ``` unzip kraken2Reports.zip ```

Let's check this directory:

```bash
s kraken2Reports|head -4
FastPrep_1.fastq.kraken2.nonames.out
FastPrep_1.fastq.kraken2.report.tsv
FastPrep_1.kraken2.species.tsv
FastPrep_2.fastq.kraken2.nonames.out
```

We have the results, now let's load the kraken2.species.tsv filtered results into a table object in R:

```R
setwd("kraken2Reports")
# List files
files <- dir(pattern = "*.kraken2.species.tsv")
# Extract sample IDs from filenames
Names <- tibble(FileName = files) %>%
  mutate(sampleID = str_remove_all(FileName, ".kraken2..*"))

# Read in all tables
Tables <- map(files, ~ read_tsv(.,
                                col_names = c("Percentage", "TaxID", "SpeciesName"),
                                col_types = cols(.default = "c"))) %>%
  set_names(Names$sampleID)

```

Transform each table: keep only Species Name and Percentage of reads
```R
Tables2Abundance <- map(Tables, ~ select(., SpeciesName, Percentage) %>%
                          mutate(Percentage = as.numeric(Percentage)))
```

This is still an object of class list we can reduce into a dataframe:

```R
AbundanceTable <- reduce(Tables2Abundance, full_join, by = "SpeciesName") %>%
  rename_with(~ Names$sampleID, -SpeciesName) %>%
  mutate(across(-SpeciesName, ~ replace_na(.x, 0)))
```

Now we can use a hierarchical clustering to compare the samples. One of the most useful tool we can use is a heatmap.

Let's use the library pheatmap (prety heatmaps) to do this:

```R
library(pheatmap)
AbundanceTable %>% 
  column_to_rownames("SpeciesName") %>%
pheatmap()
```

This is not that prety...

We can normalize the abundances to log2 to really apreciate the changes and create a matrix for the heatmap

```R
MtrixForPH <- AbundanceTable %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1))) %>%
  column_to_rownames("SpeciesName") %>%
  as.matrix()
```

And then create the heatmap:

```R
pheatmap(MtrixForPH)
```

We want to compare the treatments so the clustering in the rows (species) is not that necesary:

```R
pheatmap(MtrixForPH, cluster_row=F, cellwidth = 16)
```

Finally let's have a more "cool" color than the default:

```R
library(viridis)
Color <- rev(inferno(500))
pheatmap(MtrixForPH, 
cluster_row=F,
cellwidth = 16,
color=Color)
```

We can even add a new object and define some colors for the Methods of DNA extraction

```r
ColAnnot <- AbundanceTable %>%
  select(-SpeciesName) %>%
  colnames() %>%
  enframe(value = "Sample") %>%
  mutate(Method=str_remove(Sample,"_\\d+.*")) %>%
  select(-name) %>% 
  column_to_rownames("Sample")


pheatmap(MtrixForPH, 
         cluster_row=F,
         cellwidth = 16,
         color=Color,
         annotation_col = ColAnnot)
```

Now we can use the same Color palete to the Violinplots by creating a list

```R
AnnotColor <- list(Method=AbundanceTable %>%
                     select(-SpeciesName) %>%
                     colnames() %>%
                     enframe(value = "Sample") %>%
                     mutate(Method=str_remove(Sample,"_\\d+.*")) %>%
                     select(-name) %>% 
                     column_to_rownames("Sample") %>%
                     distinct() %>%
                     mutate(Color=brewer.pal(3,"Dark2")) %>%
                     deframe())
```

And produce our Kraken2 heatmap:

```R
KrakenPH <- pheatmap(MtrixForPH, 
         cluster_row=F,
         cellwidth = 16,
         color=Color,
         annotation_col = ColAnnot,
         annotation_colors = AnnotColor)
```

We can then combine this plot wiht the Violin plots we previous produce:

First we need to convert the pheatmap object into a ggplot2 object. The library ggplotify is very helpful:

```R
library(ggplotify)
KrakenPH <- as.ggplot(KrakenPH)
```

And then combine everything using gridExtra:

```R
library(gridExtra)

QCK2HM <- grid.arrange(arrangeGrob(ReadLengthVP,N50VP),KrakenPH,ncol=2)
```

Save into a PDF for nice report:

```R
ggsave(QCK2HM,
       file="ViolinAndHeatmap.pdf",
       width=20,
       height=20
)
```
![PP](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/ViolinAndHeatmap.png)

# Working with Bio326 Metagenomes 

Last session we found the data generated in this course BIO326_2025 it is indeed usable. So now we can follow the following pipeline to recover MAGs and predict Taxonomy (who is there?) and Functional annotation (What are they doing?).

## 1. Cleanning the reads with Chopper:

As we have 12 fastq files:

```bash
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/FastPrep_1.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/FastPrep_2.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/FastPrep_3.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/FastPrep_4.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_1.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_2.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_3.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_4.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_SRE_1.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_SRE_2.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_SRE_3.fastq.gz
/cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/Vortex_SRE_4.fastq.gz
```
We should clean this for this, instead of running sbatch 12 times we can use a useful feature of the HPC that is parallelization by [Array jobs](https://documentation.sigma2.no/jobs/job_scripts/array_jobs.html).


<details>
<summary>The following template /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/1_chopper.SLURM.sh </summary>

```bash
#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=Chopper
#
## Wall time limit:
#SBATCH --time=02:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --mem=10G
#SBATCH --gres=localscratch:20G
#SBATCH --partition=normal,bigmem,hugemem
#SBATCH --output=slurm-%x_%A_%a.out
#########################################	



#Variables
RSYNC='rsync -aPLhv --no-perms --no-owner --no-group'
arraylist=$1
INDIR=$2
OUTDIR=$3

##Main script

#####Array list######

LIST=$arraylist
input=$(head -n $SLURM_ARRAY_TASK_ID $LIST | tail -n 1)

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2022.10


##Activate conda environments

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /cluster/projects/nn9987k/.share/conda_environments/NANOPAKQC/
echo "I am workung with this" $CONDA_PREFIX

###Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID_\_$SLURM_ARRAY_TASK_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

echo "copying files to" $LOCALSCRATCH

echo "Copy fq file"

time $RSYNC $INDIR/$input.*.gz ./$input.fq.gz

###

echo "Decompress ..."

time gzip -d $input.fq.gz

###
echo "Starting QC cleanning with chooper ..."
date +%d\ %b\ %T

time cat $input.fq | chopper \
--threads $SLURM_CPUS_ON_NODE \
-q 10 \
-l 1000  > $input.chopper.fq

echo "Compressing"

pigz -p $SLURM_CPUS_ON_NODE $input.chopper.fq



###
echo "Moving files to $OUTDIR"

###Creating a directory in $OUDIR for results

time $RSYNC $input.chopper.fq.gz $OUTDIR/

#######
echo "I've done"
date

```
</details>

To run this we need to add the Terminal some arguments: 
    - A list
    - Iniput dir
    -Output dir

And we can run as follow:

```bash
cd /cluster/projects/nn9987k/$USER/
ls -1 /cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/|sed 's/.fastq.gz//g' > SampleList.txt
sbatch -a 1-12 /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/1_chopper.SLURM.sh /cluster/projects/nn9987k/BIO326-2025/metaG/SampleList.txt /cluster/projects/nn9987k/BIO326-2025/metaG/rawdata/dataBIO326_2025/rawdata/  /cluster/projects/nn9987k/$USER/metaG/results/ChopperBio326_25 && mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/ChopperBio326_25
```

This will produce a directory With the following files:

```
FastPrep_1.chopper.fq.gz  FastPrep_3.chopper.fq.gz  Vortex_1.chopper.fq.gz  Vortex_3.chopper.fq.gz  Vortex_SRE_1.chopper.fq.gz  Vortex_SRE_3.chopper.fq.gz
FastPrep_2.chopper.fq.gz  FastPrep_4.chopper.fq.gz  Vortex_2.chopper.fq.gz  Vortex_4.chopper.fq.gz  Vortex_SRE_2.chopper.fq.gz  Vortex_SRE_4.chopper.fq.gz
```

## 2. Co-Assembly reads with MetaFlye

To extend the ONT reads we will use the [Flye]() assembler with the ```--meta``` flag:

<details>
<summary>This template /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/2_flye.SLURM.chr.sh</summary>

The arguments are:

    - Input Name
    - Iniput dir
    -Output dir



```bash
#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=MetaFly
#
## Wall time limit:
#SBATCH --time=04:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 14
#SBATCH --mem=20G
#SBATCH --gres=localscratch:250G
#SBATCH --partition=normal,bigmem,hugemem
#SBATCH --output=slurm-%x_%j.out
#########################################	

###Basic usage help for this script#######

print_usage() {
        echo "Usage: sbatch $0 input indir outputdir"
}

if [ $# -lt 3 ]
        then
                print_usage
                exit 1
        fi


###############Main SCRIPT####################

##Variables###

input=$1
INDIR=$2
outdir=$3
RSYNC='rsync -a Lhv --no-perms --no-owner --no-group'



##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2022.10

##Activate conda environments

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null

conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning

echo "I'm working with this CONDAENV"
echo $CONDA_PREFIX

###Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

echo "copying Reads to" $LOCALSCRATCH

zcat $INDIR/*.gz|pigz -p $SLURM_CPUS_ON_NODE > $input.fq.gz

####Assembly#######################

echo "Starting assembly by Flye...."
date +%d\ %b\ %T

time flye \
--nano-raw $input.*gz \
--meta \
--out-dir $input.flye.outdir \
-t $SLURM_CPUS_ON_NODE

echo "Final results are in: "$outdir

$RSYNC $input.flye.outdir $outdir/

####removing tmp dir. Remember to do this for not filling the HDD in the node!!!!###

echo "I've done at"
date


```

</details>

To run we can do something like:

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/2_flye.SLURM.chr.sh MetaAssBIO326_25 /cluster/projects/nn9987k/$USER/metaG/results/ChopperBio326_25 /cluster/projects/nn9987k/$USER/metaG/results/FlyAssemblyBIO326_25 && mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/FlyAssemblyBIO326_25
```
## 3. Polishing.

A basic model of how polishing works is that the polisher stacks all relevant reads on top of the genome and decides for each position whether the present nucleotide letter is the best representative for that position, or not. There are several sources of variation that make draft assemblies polishable. The main sources are multi-strain variation from closely related species as well as incorporation of sequencing errors during the sequencing process. Ideally, assemblers would be perfect, and we wouldn't have to perform polishing. But because of some noise or artefacts that are present in our data, we might make our genomes more truthful to their biological origin by performing these polishing steps.

Genome polishing is reminiscent of generating a consensus genome. Consensus genome creation is a term used in reference mapping. This is why you may incidentally see the term consensus being used in the tools that we're gonna run.

[medaka](https://github.com/nanoporetech/medaka) is a tool to create consensus sequences and variant calls from nanopore sequencing data. This task is performed using neural networks applied a pileup of individual sequencing reads against a reference sequence, mostly commonly either a draft assembly or a database reference sequence. It provides state-of-the-art results outperforming sequence-graph based methods and signal-based methods, whilst also being faster.

**Features**

    -Requires only basecalled data. (.fasta or .fastq)
    -Improved accuracy over graph-based methods (e.g. Racon).
    -50X faster than Nanopolish (and can run on GPUs).
    -Includes extras for implementing and training bespoke correction networks.
    -Works on Linux and MacOS.
    -Open source (Oxford Nanopore Technologies PLC. Public License Version 1.0)

### Running Medaka.

Medaka needs two parameters to run:

- Fasta file of the assembly.
- Fastq files used for the assembly.


<details>

<summary> The following SLURM script /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/3_Medaka.GPU.SLURM.sh runs MEDAKA </summary>

```bash
#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=MedakaPolishingGPU
#
## Wall time limit:
#SBATCH --time=04:00:00
###Account
#SBATCH --account=nn10039k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 14
#SBATCH --mem=50G
#SBATCH --gres=localscratch:200G
#SBATCH --partition=accel,a100
#SBATCH --output=slurm-%x_%j.out
#########################################

#Variables
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'
input=$1
READIR=$2
ASSDIR=$3
OUTDIR=$4

##Main script


##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
echo "SLURM partition assigned: $SLURM_JOB_PARTITION"

# Auto-detect the partition and load the correct module environment
if [[ "$SLURM_JOB_PARTITION" == "a100" ]]; then
    echo "Running on A100 partition - Swapping to Zen2Env"
    module --force swap StdEnv Zen2Env
else
    echo "Running on Accel (P100) partition - Using default Intel environment"
fi


##Activate conda environments
module load Anaconda3/2022.10
export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null

conda activate /cluster/projects/nn9987k/.share/conda_environments/MEDAKA
echo "I am workung with this" $CONDA_PREFIX

###Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

echo "copying files to" $LOCALSCRATCH

echo "Copy fq file"

time zcat $READIR/*.gz |pigz -p $SLURM_CPUS_ON_NODE > $input.fq.gz

echo "Copy assembly"

time $RSYNC $ASSDIR/assembly.fasta ./$input.assembly.fasta

##MEdaking

echo "Starting Medaka..."
date +%d\ %b\ %T

time medaka_consensus \
-i $input.fq.gz \
-d $input.assembly.fasta \
-o $input.medaka.dir \
-t $SLURM_CPUS_ON_NODE

echo "Cleaning and changing names..."

cd $input.medaka.dir
echo "I am on:"
pwd
mv consensus.fasta $input.toto

###Cleaning

ls -1|grep -v toto|\
while read -r line; 
    do
    rm -r $line;
done

mv $input.toto $input.medaka.consensus.fasta

##moving resutls

echo "Rsync results to $OUTDIR"
cd $LOCALSCRATCH
$RSYNC $input.medaka.dir $OUTDIR/

###
echo "I've done"
date

```

</details>

We can submit it by:

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/3_Medaka.GPU.SLURM.sh FlyAssemblyBIO326_25Polished /cluster/projects/nn9987k/$USER/metaG/results/ChopperBio326_25 /cluster/projects/nn9987k/$USER/metaG/results/FlyAssemblyBIO326_25/MetaAssBIO326_25.flye.outdir /cluster/projects/nn9987k/$USER/metaG/results/MedakaPolished && mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/MedakaPolished
```

## Comparing Assemblies before and after polishing:

> [!Important]
> As we will start working with files let's ask for an interactive session in SAGA:

```bash
/cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/srun.prarameters.Nonode.Account.sh 4 10G normal,bigmem,hugemem 20G nn9987k 02:00:00
```

Now we are logged into a computing node.

Display the content of the Flye assembly folder:

```bash
ls /cluster/projects/nn9987k/$USER/metaG/results/FlyAssemblyBIO326_25/MetaAssBIO326_25.flye.outdir
```

```
00-assembly  10-consensus  20-repeat  30-contigger  40-polishing  assembly.fasta  assembly_graph.gfa  assembly_graph.gv  assembly_info.txt  flye.log  params.json
``

And the ones in MEDAKA:

```bash

ls /cluster/projects/nn9987k/$USER/metaG/results/MedakaPolished/FlyAssemblyBIO326_25Polished.medaka.dir

```

```
FlyAssemblyBIO326_25Polished.medaka.consensus.fasta
```

As we can see these are fasta files, let's corroborate these are fasta files:

1) Let's assign this into a variables to easy manipulate files:

```bash
FLYE="/cluster/projects/nn9987k/$USER/metaG/results/FlyAssemblyBIO326_25/MetaAssBIO326_25.flye.outdir/assembly.fasta"
MEDAKA="/cluster/projects/nn9987k/$USER/metaG/results/MedakaPolished/FlyAssemblyBIO326_25Polished.medaka.dir/FlyAssemblyBIO326_25Polished.medaka.consensus.fasta"
```

To check we can use ```less```

```bash
less -S $MEDAKA
less -S $FLYE

```

We can use ```assembly-stats``` tool to check for the main stats in the assembly:

```bash
module load Miniconda3/23.10.0-1
eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning/
assembly-stats $FLYE $MEDAKA
```

```
stats for /cluster/projects/nn9987k/BIO326-2025/metaG/results/FlyAssemblyBIO326_25/MetaAssBIO326_25.flye.outdir/assembly.fasta
sum = 63963185, n = 2515, ave = 25432.68, largest = 1763124
N50 = 44610, n = 256
N60 = 31296, n = 431
N70 = 22114, n = 673
N80 = 15765, n = 1014
N90 = 10695, n = 1504
N100 = 508, n = 2515
N_count = 0
Gaps = 0
-------------------------------------------------------------------------------
stats for /cluster/projects/nn9987k/BIO326-2025/metaG/results/MedakaPolished/FlyAssemblyBIO326_25Polished.medaka.dir/FlyAssemblyBIO326_25Polished.medaka.consensus.fasta
sum = 63718690, n = 2515, ave = 25335.46, largest = 1763100
N50 = 44616, n = 255
N60 = 31117, n = 429
N70 = 22045, n = 672
N80 = 15701, n = 1012
N90 = 10639, n = 1503
N100 = 508, n = 2515
N_count = 0
Gaps = 0
```

These results are good but not very comparable ```assembly-stats``` can perform an output as a table using the ```-t``` flag:

```
assembly-stats -t $FLYE $MEDAKA

```

We can save this into a file:

```bash
assembly-stats -t $FLYE $MEDAKA > /cluster/projects/nn9987k/$USER/metaG/results/Flye.Medaka.stats.tsv
cat !$
```

```
cat /cluster/projects/nn9987k/$USER/metaG/results/Flye.Medaka.stats.tsv
filename        total_length    number  mean_length     longest shortest        N_count Gaps    N50     N50n    N70     N70n    N90     N90n
/cluster/projects/nn9987k/BIO326-2025/metaG/results/FlyAssemblyBIO326_25/MetaAssBIO326_25.flye.outdir/assembly.fasta    63963185        2515    25432.68        1763124 5080       0       44610   256     22114   673     10695   1504
/cluster/projects/nn9987k/BIO326-2025/metaG/results/MedakaPolished/FlyAssemblyBIO326_25Polished.medaka.dir/FlyAssemblyBIO326_25Polished.medaka.consensus.fasta  63718690  2515     25335.46        1763100 508     0       0       44616   255     22045   672     10639   1503
```

We can plot this result table:

```bash
conda activate /cluster/projects/nn9987k/.share/conda_environments/R_env/
cd /cluster/projects/nn9987k/$USER/metaG/results/
Rscript /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/assemblyStats.r  /cluster/projects/nn9987k/$USER/metaG/results/Flye.Medaka.stats.tsv
```

This producess the plot:

```
Plot saved as AssemblyStats.pdf
```

![AssemblyStats](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/AssStats.JPG)

**What can we say about these results?**

> [!Important]
> Remember to finish your interactive session by ```exit```