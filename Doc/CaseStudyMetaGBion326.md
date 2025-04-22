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
```

<details>

```

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

</details>

and use a loop:

```bash
 cat barcodes.tsv |while read -r line; do NAME=$(echo $line|awk '{print $1}'); BC=$(echo $line|awk '{print $2}'); echo $BC; zcat fastq_pass/$BC/*fastq.gz |pigz -p 12 >  rawdata/$NAME.fastq.gz; done
```
Then we will end with something like:

```bash
ls -1 rawdata/
```
<details>

```
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
</details>

Then using the usefull NanoPlot:

```bash
 parallel -j 12 "NanoPlot --fastq {} --N50 --loglength -o {}.Nanoplot.dir" ::: *.gz
```

And as we have been working we can collect the NanoStats.txt

<details>

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

</details>

Let's work with these files. Download the files from 

https://arken.nmbu.no/~auve/BIO326/bio326NanoStats.Prok.dir.zip



>[!Tip]
> If you have access to the terminal you can download the file by:
> ``` wget https://arken.nmbu.no/~auve/BIO326/bio326NanoStats.Prok.dir.zip ```

Decompress the Zip file and let's move to R and RStudio.
>[!Tip]
> If you have access to the terminal you unzip the file by:
> ``` unzip bio326NanoStats.Prok.dir.zip ```

<details>

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

</details>


- **Using RStiudio to load these files:**

<details>

<summary>This chunk of code will generate some useful plots for us:</summary>

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

# Pivot for easy ploting

stats_wide <- stats_data %>%
  pivot_wider(names_from = Metric, values_from = Value)
head(stats_wide,2)


# Convert Sample column to factor to maintain order

stats_wide <- stats_wide %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Let's create a function to plot all the metrics:
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

# Check the results, by displaying the Mean read length and the N50


plot_list$Numberofreads
plot_list$ReadlengthN50


# We can use Patchwork to plot in the same canvas:

library(patchwork)
plot_list$Numberofreads + plot_list$ReadlengthN50
```

</details>

As we can see the Bar plot is not the best plot to display and compare these results...

>[!Tip]
> Let's use a violin plot:

<details>

<summary> R code </summary>

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
</details>

We can use some statistics like ANOVA and t-test to check if there is differences in the groups. The library ggpubr allows us to calculate these kind of statistics and visualize on the ggplot2 object:
<details>

<summary>R code</summary>

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
</details>

And modify some Color and astetics of the plot:

<details>

<summary>R code</summary>

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
</details>

Saving the plot:

<details>
<summary>R code</summary>

```R
MergedQCPlot <- ReadLengthVP + N50VP

ggsave(MergeQCPlot,file="MergedPlot.QualityScores.pdf")

```
</details>

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
```

<details>

```
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

</details>

We can parse this report and only keep the organisms with > 0.1 % of reads classified:

```bash
head  FastPrep_1.kraken2.species.tsv
 
```
<details>

```
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

</details>


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
ls kraken2Reports|head -4

```

<details>

```
FastPrep_1.fastq.kraken2.nonames.out
FastPrep_1.fastq.kraken2.report.tsv
FastPrep_1.kraken2.species.tsv
FastPrep_2.fastq.kraken2.nonames.out

```

</details>


We have the results, now let's load the kraken2.species.tsv filtered results into a table object in R:

<details>
<summary>R code</summary>

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

</details>

Transform each table: keep only Species Name and Percentage of reads

<details>
<summary>R code</summary>

```R
Tables2Abundance <- map(Tables, ~ select(., SpeciesName, Percentage) %>%
                          mutate(Percentage = as.numeric(Percentage)))
```

</details>

This is still an object of class list we can reduce into a dataframe:

<details>
<summary>R code</summary>

```R
AbundanceTable <- reduce(Tables2Abundance, full_join, by = "SpeciesName") %>%
  rename_with(~ Names$sampleID, -SpeciesName) %>%
  mutate(across(-SpeciesName, ~ replace_na(.x, 0)))
```

</details>

Now we can use a hierarchical clustering to compare the samples. One of the most useful tool we can use is a heatmap.

Let's use the library pheatmap (prety heatmaps) to do this:

<details>
<summary>R code</summary>

```R
library(pheatmap)
AbundanceTable %>% 
  column_to_rownames("SpeciesName") %>%
pheatmap()
```

</details>

This is not that prety...

We can normalize the abundances to log2 to really apreciate the changes and create a matrix for the heatmap

<details>
<summary>R code </summary>

```R
MtrixForPH <- AbundanceTable %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1))) %>%
  column_to_rownames("SpeciesName") %>%
  as.matrix()
```

</details>

And then create the heatmap:

<details>
<summary>R code</summary>

```R
pheatmap(MtrixForPH)
```

</details>

We want to compare the treatments so the clustering in the rows (species) is not that necesary:

<details>
<summary>R code</summary>

```R
pheatmap(MtrixForPH, cluster_row=F, cellwidth = 16)
```

</details>

Finally let's have a more "cool" color than the default:

<details>
<summary>R code</summary>

```R
library(viridis)
Color <- rev(inferno(500))
pheatmap(MtrixForPH, 
cluster_row=F,
cellwidth = 16,
color=Color)
```

</details>

We can even add a new object and define some colors for the Methods of DNA extraction

<details>
<summary>R code</summary>

```R
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
</details>

Now we can use the same Color palete to the Violinplots by creating a list

<details>
<summary>R code</summary>

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

</details>

And produce our Kraken2 heatmap:

<details>
<summary>R code</summary>

```R
KrakenPH <- pheatmap(MtrixForPH, 
         cluster_row=F,
         cellwidth = 16,
         color=Color,
         annotation_col = ColAnnot,
         annotation_colors = AnnotColor)
```

</details>

We can then combine this plot wiht the Violin plots we previous produce:

First we need to convert the pheatmap object into a ggplot2 object. The library ggplotify is very helpful:

<details>
<summary>R code</summary>

```R
library(ggplotify)
KrakenPH <- as.ggplot(KrakenPH)
```

</details>

And then combine everything using gridExtra:

<details>
<summary>R code</summary>

```R
library(gridExtra)

QCK2HM <- grid.arrange(arrangeGrob(ReadLengthVP,N50VP),KrakenPH,ncol=2)
```

</details>

Save into a PDF for nice report:

<details>
<summary>R code</summary>
```R
ggsave(QCK2HM,
       file="ViolinAndHeatmap.pdf",
       width=20,
       height=20
)
```

</details>

![PP](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/ViolinAndHeatmap.png)

# The Bio326 Metagenomes: A tale of who is there and what are they doing?

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

<details>

```
FastPrep_1.chopper.fq.gz  FastPrep_3.chopper.fq.gz  Vortex_1.chopper.fq.gz  Vortex_3.chopper.fq.gz  Vortex_SRE_1.chopper.fq.gz  Vortex_SRE_3.chopper.fq.gz
FastPrep_2.chopper.fq.gz  FastPrep_4.chopper.fq.gz  Vortex_2.chopper.fq.gz  Vortex_4.chopper.fq.gz  Vortex_SRE_2.chopper.fq.gz  Vortex_SRE_4.chopper.fq.gz
```

</details>

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

### Comparing Assemblies before and after polishing:

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

## 4. Binning.

So far we created an assembly containing all contigs (continuous sequences) from each of the organisms in the microbial community that we sequenced from the cow rumen.

Presently, the assembly consists of thousands of contigs, each coming from a single species. By grouping together the contigs from each species present in our sample, we can create what is referred to as a MAG, short for Metagenome Assembled Genome.

Popular binning algorithms like the ones used in Metabat2 utilize contig depth as a tell tale to link the individual contigs together that come from the same species. This is done by mapping the original reads onto the assembly and then counting the read depth of each contig. The smart thing here is that contigs coming from the same species will have similar depth. Another vital contig statistic that binners use is the GC-content. Each species has its own intrinsic GC-content, and by grouping contigs further on GC-content -in this case by counting the tetranucleotide frequency- we might get good representatives for distinct species in our sample. If our bins live up to our requirements, we can refer to them as MAGs.

![Bining](https://github.com/TheMEMOLab/Bin420-Bioinformatics-for-Functional-Meta-Omics/blob/main/img/Binning.png)


>[!Note]
> In this course we will use 2 binning algorithms [Metabat2](https://bitbucket.org/berkeleylab/metabat/src) and  [Maxbin2](https://sourceforge.net/projects/maxbin2/)

Both algorithms relay in extracting the sequencing depth from the assemlby using a table like this:

<details>

```
contigName      contigLen       totalAvgDepth   MetaBiningBIO326_25Polished.assembly.sorted     MetaBiningBIO326_25Polished.assembly.sorted-var
contig_1784     16061   0       0       0
contig_2114     83571   3.45343 3.45343 2.98818
contig_2646     14442   2.83193 2.83193 0.380264
contig_1102     3538    2.63017 2.63017 0.516515
```

</details>

### Running Binning tools:

Where you need the length of each contig and the deepth sequenced in each experiment. To do this we will map all the reads to the MEDAKA polished assembly using ```minimap2``` and then we will use the ```jgi_summarize_bam_contig_depths``` script to get the depth file.

<details>

<summary>This SLURM script has all the instructions</summary>


```bash
#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=Binning
#
## Wall time limit:
#SBATCH --time=48:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=50G
#SBATCH --gres=localscratch:200G
#SBATCH --output=slurm-%x_%j.out
#########################################	

#Variables
RSYNC='rsync -aPLhv --no-perms --no-owner --no-group'
input=$1
READIR=$2
ASSDIR=$3
OUTDIR=$4

##Activate conda environments ## Arturo

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning/
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

##createing a directory for Metabat
mkdir $input.Binning.dir
cd $input.Binning.dir

echo "Copy fq file"

time $RSYNC $READIR/$input.chopper.fq.gz .

echo "Copy assembly"

time $RSYNC $ASSDIR/$input.medaka.dir/$input.medaka.consensus.fasta ./$input.assembly.fasta

##Align

echo "Start minimap2 "
date +%d\ %b\ %T

time minimap2 \
-ax map-ont \
-t $SLURM_CPUS_ON_NODE \
$input.assembly.fasta \
$input.chopper.fq.gz > $input.assembly.sam

#Get the bam and the files

echo "Samtools view.."
 
time samtools view \
-@ $SLURM_CPUS_ON_NODE \
-bS $input.assembly.sam > $input.assembly.bam

rm -r $input.assembly.sam
rm -r $input.chopper.fq.gz

echo "samtools sort"

time samtools sort \
-@ $SLURM_CPUS_ON_NODE $input.assembly.bam > $input.assembly.sorted.bam

rm -r $input.assembly.bam

###Binning

echo "Calculating depth..."

time jgi_summarize_bam_contig_depths \
        --outputDepth $input.depth.txt \
        $input.assembly.sorted.bam

echo "Modifying $input.depth.txt to a MaxBin2 format:"

cut -f1,3 $input.depth.txt | tail -n+2 > $input.depth_maxbin.txt  

##Using Parallel to run metaba2 and Maxbin2 at the same time with 8 CPUs each...

cpu=$(($SLURM_CPUS_ON_NODE/2))
echo "Running parallel with $cpu cpus"

time parallel -j2 ::: "metabat2 \
        -i $input.assembly.fasta \
        -a $input.depth.txt \
        -m 1500 \
        --seed 100 \
        -t $cpu \
        --unbinned \
        -o $input.Metabat2" \
        "run_MaxBin.pl \
        -contig $input.assembly.fasta \
        -out $input.MaxBin.out \
        -abund $input.depth_maxbin.txt \
        -thread $cpu"

#Changing suffix fa to fasta useful for further analysis

for i in *.fa;
        do
        a=$(basename $i .fa);
        mv $i $a.fasta;
done

#Remove the original assemlby

rm $input.assembly.fasta

###
echo "Moving Binnig files to $OUTDIR"

cd $LOCALSCRATCH

time $RSYNC $input.Binning.dir $OUTDIR/

#######
echo "I've done"
date

```

</details>

The script needs 4 arguments:
- input=$1  Name of the sample
- READIR=$2 Directory of the Choppered reads 
- ASSDIR=$3 Medaka directory
- OUTDIR=$4 Output directry

Let's run it:

```bash
 sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/4_Binning.SLURM.sh MetaBiningBIO326_25Polished /cluster/projects/nn9987k/$USER/metaG/results/ChopperBio326_25 /cluster/projects/nn9987k/$USER/metaG/results/MedakaPolished/FlyAssemblyBIO326_25Polished.medaka.dir  /cluster/projects/nn9987k/$USER/metaG/results/MetaBiningBIO326_25Polished && mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/MetaBiningBIO326_25Polished
```

After running for 20 min you will end up with a folder like:

```bash
cd /cluster/projects/nn9987k/auve/metaG/results/MetaBiningBIO326_25Polished/MetaBiningBIO326_25Polished.Binning.dir
ls
```
<details>

```

MetaBiningBIO326_25Polished.assembly.sorted.bam   MetaBiningBIO326_25Polished.MaxBin.out.019.fasta                  MetaBiningBIO326_25Polished.Metabat2.23.fasta
MetaBiningBIO326_25Polished.depth_maxbin.txt      MetaBiningBIO326_25Polished.MaxBin.out.log                        MetaBiningBIO326_25Polished.Metabat2.24.fasta
MetaBiningBIO326_25Polished.depth.txt             MetaBiningBIO326_25Polished.MaxBin.out.marker                     MetaBiningBIO326_25Polished.Metabat2.25.fasta
MetaBiningBIO326_25Polished.MaxBin.out.001.fasta  MetaBiningBIO326_25Polished.MaxBin.out.marker_of_each_bin.tar.gz  MetaBiningBIO326_25Polished.Metabat2.26.fasta
MetaBiningBIO326_25Polished.MaxBin.out.002.fasta  MetaBiningBIO326_25Polished.MaxBin.out.noclass                    MetaBiningBIO326_25Polished.Metabat2.27.fasta
MetaBiningBIO326_25Polished.MaxBin.out.003.fasta  MetaBiningBIO326_25Polished.MaxBin.out.summary                    MetaBiningBIO326_25Polished.Metabat2.28.fasta
MetaBiningBIO326_25Polished.MaxBin.out.004.fasta  MetaBiningBIO326_25Polished.MaxBin.out.tooshort                   MetaBiningBIO326_25Polished.Metabat2.29.fasta
MetaBiningBIO326_25Polished.MaxBin.out.005.fasta  MetaBiningBIO326_25Polished.Metabat2.10.fasta                     MetaBiningBIO326_25Polished.Metabat2.2.fasta
MetaBiningBIO326_25Polished.MaxBin.out.006.fasta  MetaBiningBIO326_25Polished.Metabat2.11.fasta                     MetaBiningBIO326_25Polished.Metabat2.30.fasta
MetaBiningBIO326_25Polished.MaxBin.out.007.fasta  MetaBiningBIO326_25Polished.Metabat2.12.fasta                     MetaBiningBIO326_25Polished.Metabat2.31.fasta
MetaBiningBIO326_25Polished.MaxBin.out.008.fasta  MetaBiningBIO326_25Polished.Metabat2.13.fasta                     MetaBiningBIO326_25Polished.Metabat2.3.fasta
MetaBiningBIO326_25Polished.MaxBin.out.009.fasta  MetaBiningBIO326_25Polished.Metabat2.14.fasta                     MetaBiningBIO326_25Polished.Metabat2.4.fasta
MetaBiningBIO326_25Polished.MaxBin.out.010.fasta  MetaBiningBIO326_25Polished.Metabat2.15.fasta                     MetaBiningBIO326_25Polished.Metabat2.5.fasta
MetaBiningBIO326_25Polished.MaxBin.out.011.fasta  MetaBiningBIO326_25Polished.Metabat2.16.fasta                     MetaBiningBIO326_25Polished.Metabat2.6.fasta
MetaBiningBIO326_25Polished.MaxBin.out.012.fasta  MetaBiningBIO326_25Polished.Metabat2.17.fasta                     MetaBiningBIO326_25Polished.Metabat2.7.fasta
MetaBiningBIO326_25Polished.MaxBin.out.013.fasta  MetaBiningBIO326_25Polished.Metabat2.18.fasta                     MetaBiningBIO326_25Polished.Metabat2.8.fasta
MetaBiningBIO326_25Polished.MaxBin.out.014.fasta  MetaBiningBIO326_25Polished.Metabat2.19.fasta                     MetaBiningBIO326_25Polished.Metabat2.9.fasta
MetaBiningBIO326_25Polished.MaxBin.out.015.fasta  MetaBiningBIO326_25Polished.Metabat2.1.fasta                      MetaBiningBIO326_25Polished.Metabat2.lowDepth.fasta
MetaBiningBIO326_25Polished.MaxBin.out.016.fasta  MetaBiningBIO326_25Polished.Metabat2.20.fasta                     MetaBiningBIO326_25Polished.Metabat2.tooShort.fasta
MetaBiningBIO326_25Polished.MaxBin.out.017.fasta  MetaBiningBIO326_25Polished.Metabat2.21.fasta                     MetaBiningBIO326_25Polished.Metabat2.unbinned.fasta

```

</details>

Let's count how many bins did we recover from MaxBin and howmany from Metabat2:

```bash
MB=$(ls -1|grep Metabat2|grep -v -E "lowDept|tooShort|unbin"|wc -l)
MX=$(ls -1|grep MaxBin|grep -v -E "log|marker|nocl|summ|too"|wc -l)

echo -e "Metabat2\t$MB\nMaxBin2\t$MX"

```
>[!Important]
>The number of Bins in each run (user) can change due to the binning algorithm and lack of seed setting options in MaxBin2.

## 5. Dereplication 

As we use 2 different binning strategies and we could be dealing with recoveing same organissm recovered in duplicated MAGs.

We will use then [dREP](https://drep.readthedocs.io/en/latest/index.html#) tool to dereplicate and get the best bin for each binning tool.

<details>

<summary> We can use this SLURM template: </summary>

```bash

#!/bin/bash

###################################
## Job name:
#SBATCH --job-name=dRep
#
## Wall time limit:
#SBATCH --time=24:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --cpus-per-task 16
#SBATCH --mem=120G
#SBATCH --gres=localscratch:150G
#SBATCH --partition=bigmem
#SBATCH --out slurm-%x_%j.out
#######################################


## Set up job environment:
#set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error

###Variables###
input=$1 #input string
indir=$2 #Input directory
ext=$3 #extension of fasta file e.g fasta
comp=$4 #completeness score 
con=$5 #Contamination score
outdir=$6 #output directory

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1

##Activate conda environments

export PS1=\$
eval "$(/cluster/software/Miniconda3/23.10.0-1/bin/conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/DEREPLICATION
echo "I am workung with this" $CONDA_PREFIX

####Do some work:########

echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "I am working with this enviroment loaded"
echo $CONDA_PREFIX
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

workdir=$(pwd)
echo "My working directory is" $workdir
echo "copying MAGs files ..."
time rsync -aL $indir/*.$ext .

##Check for complete contigs and unbined files and remove it from the list

unbin=$(ls -1 | grep "unbin")
contigs=$(ls -1 | grep "contigs.fasta")
lowDepth=$(ls -1 | grep "lowDepth.fasta")
tooshort=$(ls -1 | grep "tooShort.fasta")

# Array of file variables
files=("$contigs" "$unbin" "$lowDepth" "$tooshort")

# Loop through files and check existence
for file in "${files[@]}"; do
    if [[ -f "$file" ]]; then
        echo "File $file found, removing it..."
        rm "$file"
    fi
done

###

echo "Number of genomes to drep:"
ls -1|wc -l

mkdir $input.MAGs
mv *.$ext $input.MAGs/

#####dREP dereplicate pipeline################
echo "start dREP at"
date +%d\ %b\ %T

time dRep \
	dereplicate \
	$input.DREP.$comp.$con.out \
	-g $input.MAGs/*.$ext \
	-p $SLURM_CPUS_ON_NODE \
	-comp $comp \
	-con $con

##Copy files to the $SLURM_SUBMIT_DIR
cd $workdir

time rsync -aP $input.DREP.$comp.$con.out $outdir

echo "results are in: " $outdir/$input.DREP.$comp.$con.out

###
echo "I've done at"
date

```

</details>

The script requires  of parameters:

```
input=$1 #input string
indir=$2 #Input directory
ext=$3 #extension of fasta file e.g fasta
comp=$4 #completeness score 
con=$5 #Contamination score
outdir=$6 #output directory
```

To run:

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/5_drep.SLURM.sh DEREP_BIO326_25 /cluster/projects/nn9987k/$USER/metaG/results/MetaBiningBIO326_25Polished/MetaBiningBIO326_25Polished.Binning.dir fasta 70 5 /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION && mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION
```

After running we should end with a file structure like this:

```bash
tree -d -L 2 /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out
```

<details>

```

/cluster/projects/nn9987k/auve/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out
├── data
│   ├── checkM
│   ├── Clustering_files
│   ├── fastANI_files
│   ├── MASH_files
│   └── prodigal
├── data_tables
├── dereplicated_genomes
├── figures
└── log

11 directories

```


</details>


let's check how many MAGs were Dereplicated with > 70 % completeness and < 5 % contamination.

```bash

tree  /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/dereplicated_genomes

```

<details>

```
/cluster/projects/nn9987k/auve/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/dereplicated_genomes
├── MetaBiningBIO326_25Polished.MaxBin.out.001.fasta
├── MetaBiningBIO326_25Polished.Metabat2.11.fasta
├── MetaBiningBIO326_25Polished.Metabat2.14.fasta
├── MetaBiningBIO326_25Polished.Metabat2.27.fasta
└── MetaBiningBIO326_25Polished.Metabat2.8.fasta

0 directories, 5 files
```

</details>

We can go deeper and check how the DREPLICATION selected these Five genomes. First, let's check the CHECKM results:

```bash
ls  /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/data/checkM/checkM_outdir

```

<details>

```
bins  Chdb.tsv  checkm.log  lineage.ms  results.tsv  storage
```
</details>

we can take a look on the results by ```less results.tsv```

```bash
less  /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/data/checkM/checkM_outdir/results.tsv
```

This is a huge table, so let's just extract the fields we need the GenomeID (1), taxonomy (2), completeness (12) and contamination (13), and ask to retrieve only those that are > 70 % complete and < 5 % contaminated:

```bash
 cat /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/data/checkM/checkM_outdir/results.tsv |awk -F "\t" '{if($12 > 70  && $13 < 5) print $1,$2,$12,$13}'
```
<details>
```
MetaBiningBIO326_25Polished.MaxBin.out.001.fasta o__Clostridiales (UID1226) 98.57 0.00
MetaBiningBIO326_25Polished.MaxBin.out.002.fasta k__Bacteria (UID2372) 99.06 2.83
MetaBiningBIO326_25Polished.MaxBin.out.003.fasta g__Prevotella (UID2722) 97.81 4.26
MetaBiningBIO326_25Polished.MaxBin.out.006.fasta f__Lachnospiraceae (UID1255) 81.90 4.02
MetaBiningBIO326_25Polished.Metabat2.1.fasta o__Clostridiales (UID1226) 95.38 0.00
MetaBiningBIO326_25Polished.Metabat2.11.fasta k__Bacteria (UID2372) 99.06 0.94
MetaBiningBIO326_25Polished.Metabat2.14.fasta k__Bacteria (UID2372) 96.86 0.94
MetaBiningBIO326_25Polished.Metabat2.27.fasta g__Prevotella (UID2722) 95.32 0.76
MetaBiningBIO326_25Polished.Metabat2.8.fasta f__Lachnospiraceae (UID1255) 90.21 3.85
```
</details>


Here we have more than 5 (n=9) genomes, but dRep says only 5. Well check the ANI clustering analysis:

- CheckM save all the clustering analysis as dendrograms and NMDS plots in the ```figures``` folder as PDF. 

>[!Tip]
> We can open PDF directly in SAGA using the ```evince``` command. However, for this the X11 Display should be on in your computers. You can read more how to enable this in VS-code here [VS-CODE X Server](https://x410.dev/cookbook/enabling-ssh-x11-forwarding-in-visual-studio-code-for-remote-development/)

```bash
ls  /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/figures
```

<details>
```
Clustering_scatterplots.pdf  Primary_clustering_dendrogram.pdf     Secondary_clustering_MDS.pdf
Cluster_scoring.pdf          Secondary_clustering_dendrograms.pdf  Winning_genomes.pdf
```
</details>

Then using ```evince``` command to open the Primary_clustering_dendrogram.pdf file

```bash
evince /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/figures/Primary_clustering_dendrogram.pdf
```

![METAG](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/ANI.png)

- **How many dereplicated MAGs could we obtain if we change the parameters to 50 % Completeness 10 % Contamination?**

## 6. Taxonomy and functional annotation (Who are there and what can they do?)

### Running DRAM: Distilled and Refined Annotation of Metabolism

"[DRAM](https://github.com/WrightonLabCSU/DRAM) (Distilled and Refined Annotation of Metabolism) is a tool for annotating metagenomic assembled genomes and VirSorter identified viral contigs. DRAM annotates MAGs and viral contigs using KEGG (if provided by the user), UniRef90, PFAM, dbCAN, RefSeq viral, VOGDB and the MEROPS peptidase database as well as custom user databases..."



<img src="https://github.com/avera1988/NMBU-Bio-326/blob/main/images/DRAM.jpg" height="400">




<details>

<summary>This is the DRAM SLURM template:</summary>

```bash

#!/bin/bash
#########################################################################
#	SLURM scrip for running DRAM annotator on SAGA cluster using all MAGs
#########################################################################

###############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=DRAM
#
## Wall time limit:
#SBATCH --time=90:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=80G
#SBATCH --gres=localscratch:50G
#SBATCH --partition=bigmem
#SBATCH --out slurm-%x-%A.out

###########################################################

##########Variables

dir=$1 ##directory with fasta files
ext=$2 #Extension of the fasta files e.g. fasta
OUTDIR=$3 #Outputdirectory
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"

conda activate /cluster/projects/nn9987k/.share/conda_environments/DRAM/ 

####Do some work:########

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

#Copy the MAGs to the $LOCALSCRATCH for local computation

echo "copying MAGs to" $LOCALSCRATCH

time $RSYNC $dir/ ./MAGs


##################DRAM##############################

echo "DRAM for annotation at"
date +%d\ %b\ %T

time DRAM.py annotate \
-i MAGs'/*.'$ext \
-o dram.annotation.dir \
--min_contig_size 500 \
--threads $SLURM_CPUS_ON_NODE

echo "Distilling ..."
date +%d\ %b\ %T  

time DRAM.py distill \
-i dram.annotation.dir/annotations.tsv \
-o dram.genome_summaries.dir \
--trna_path dram.annotation.dir/trnas.tsv \
--rrna_path dram.annotation.dir/rrnas.tsv
echo "DRAM finished at"
date +%d\ %b\ %T

##Generate a master directory to keep all the results together 

mkdir DRAM.Results.dir
mv dram.annotation.dir DRAM.Results.dir
mv dram.genome_summaries.dir DRAM.Results.dir

###########Moving results to $SLURM_SUBMIT_DIR partition or anywhere the main script was submitted############

echo "moving results to" $OUTDIR/

cd $LOCALSCRATCH

time $RSYNC DRAM.Results.dir $OUTDIR/

echo "DRAM results are in: " $OUTDIR/DRAM.Results.dir


```

</details>

DRAM scipt requires the following arguments:

- dir=$1 ##directory with fasta files
- ext=$2 #Extension of the fasta files e.g. fasta
- OUTDIR=$3 #Outputdirectory

Use the following sbatch line to run:

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/6a_DRAM.SLURM.sh\
/cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/dereplicated_genomes \
fasta \
/cluster/projects/nn9987k/$USER/metaG/results/DRAM && \
mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/DRAM
```

>[!Warning]
> DRAM requires a lot of time (~4-6hrs) and resources (~80-500Gb RAM) to run, so plan accordingly for submitting the job.

After finish you will end up with something like:

```bash
tree /cluster/projects/nn9987k/$USER/metaG/results/DRAM
```

<details>

```
/cluster/projects/nn9987k/auve/metaG/results/DRAM
└── DRAM.Results.dir
    ├── dram.annotation.dir
    │   ├── annotate.log
    │   ├── annotations.tsv
    │   ├── genbank
    │   │   ├── MetaBiningBIO326_25Polished.MaxBin.out.001.gbk
    │   │   ├── MetaBiningBIO326_25Polished.Metabat2.21.gbk
    │   │   ├── MetaBiningBIO326_25Polished.Metabat2.27.gbk
    │   │   ├── MetaBiningBIO326_25Polished.Metabat2.32.gbk
    │   │   └── MetaBiningBIO326_25Polished.Metabat2.5.gbk
    │   ├── genes.faa
    │   ├── genes.fna
    │   ├── genes.gff
    │   ├── rrnas.tsv
    │   ├── scaffolds.fna
    │   └── trnas.tsv
    └── dram.genome_summaries.dir
        ├── distill.log
        ├── genome_stats.tsv
        ├── metabolism_summary.xlsx
        ├── product.html
        └── product.tsv

4 directories, 18 files
```

</details>

There are two directories: 
* dram.annotation.dir: It has all the "raw" annotations, gene sequences, protein preditions of the MAGs
* dram.genome_summaries.dir: It has the distilled part of the genomes with the sorted metabolic functions.

**dram.annotation.dir**

Here we can find a table with annotations (annotations.tsv) as well as 3 fasta files:

- genes.fna (All the predicted coding genes as nucleotides)
- genes.faa (All the predicted coding genes translated to proteins)
- scaffolds.fna (The total scaffolds/contigs of the bins)

**dram.genome_summaries.dir**

Although we can display the content of the *.tsv* files obtained by DRAM here in the terminal, the metabolism_summary.xlsx and product.html files are visually friendly, so it is recommended to export these to our personal computers and take a look. 


### CompareM2

An additional and brand new tool for annotation created here at the MEMO group is [CompareM2](https://comparem2.readthedocs.io/en/latest/).

According to the readthedocs:

**"🧬 CompareM2 is a genomes-to-report pipeline. It accepts prokaryotic (bacterial and archaeal) genomic assemblies and compares them in many different ways.

🦠 Being designed to analyze assemblies of both isolates and metagenomes (MAGs), it is useful for anyone working with microbial genomics. "

<details>

<summary>This is the SLURM template for CompareM2</summary>

```bash
#!/bin/bash
#########################################################################

###############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=COMPAREM2
#
## Wall time limit:
#SBATCH --time=90:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=80G
#SBATCH --gres=localscratch:200G
#SBATCH --partition=bigmem
#SBATCH --out slurm-%x-%A.out

###########################################################

##########Variables

dir=$1 ##directory with fasta files
OUTDIR=$2 #Outputdirectory
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'


##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
mmodule load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"

conda activate /cluster/projects/nn9987k/.share/conda_environments/COMPAREM2/ 

####Do some work:########

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

#Copy the MAGs to the $LOCALSCRATCH for local computation

echo "copying MAGs to" $LOCALSCRATCH

time $RSYNC $dir/ ./MAGs


##################COMPAREM2##############################
########Configuration of COMPAREM###############

export COMPAREM2_DATABASES="/cluster/projects/nn9987k/.share/db/COMPAREM2"
export APPTAINER_TMPDIR=$LOCALSCRATCH
export APPTAINER_CACHEDIR=$LOCALSCRATCH

echo "This configuration is running:"

echo "DB" $COMPAREM2_DATABASES
echo "APTTMPDID" $APPTAINER_TMPDIR
echo "APTCACHE" $APPTAINER_CACHEDIR

echo "Starting COMPAREM2"
date +%b\ %d\ %T

time comparem2 \
--cores $SLURM_CPUS_ON_NODE \
--use-singularity \
--singularity-prefix $(pwd) \
--config input_genomes="$(pwd)/MAGs/*.fasta" output_directory="CompareM.out.dir" \
--until assembly_stats checkm2 prokka gtdbtk

###########Moving results to $SLURM_SUBMIT_DIR partition or anywhere the main script was submitted############

echo "moving results to" $OUTDIR/

cd $LOCALSCRATCH

time $RSYNC CompareM.out.dir $OUTDIR/

echo "COMPAREM results are in: " $OUTDIR/CompareM.out.dir

```
</details>

This script requires the following arguments:

- dir=$1 ##directory with fasta files
- OUTDIR=$2 #Outputdirectory

**Running compareM2**

```bash
sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/6b_CompareM2.SLURM.sh /
/cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/dereplicated_genomes \
/cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2 && \
mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2
```

>[!Warning]
> compareM2 also requires a lot of resources and time, so again plan accordingly...

When this is done you should end-up with something like:

```bash
tree -d -L 2 /cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2
```

<details>

```
/cluster/projects/nn9987k/auve/metaG/results/COMPAREM2
└── CompareM.out.dir
    ├── assembly-stats
    ├── benchmarks
    ├── checkm2
    ├── gtdbtk
    ├── samples
    └── tables

7 directories

```

</details>


The main important result is the html report, that is here:

```
cd /cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2/ CompareM.out.dir/
ls
```

<details>

```
assembly-stats  benchmarks  checkm2  gtdbtk  metadata.tsv  report_14391991.html  samples  tables  versions.txt
```

</details>

>[!Note]
> Let's copy the report ot our local PC and take a look.

## 7. Visualization.

### 7.1 Phylogenomics:

> [!Tip]
> Taken from a book chapter published here at the MEMO group [Long-Read Metagenomics and CAZyme Discovery](https://link.springer.com/protocol/10.1007/978-1-0716-3151-5_19)

After assessment of quality (with CheckM2), taxonomic (with GTDB-Tk) and functional annotation (with DRAM and Prokka) of the MAGs, a visual representation of the phylogenetic clusterization (e.g., phylogenetic tree) and the annotations (e.g., heatmaps) are often desirable.

Different software can be used to generate phylogenetic trees of the MAGs, including PhyloPhlAn, Anv'o, and IQ-TREE. At the MEMO lab we have also created an easy Rscript that can combine the information from GTDBTk, PhyloPhlAn and Checkm2 to produce a phylogenetic tree with different layers. In the GitHub repo [MetaGVisualToolBox](https://github.com/TheMEMOLab/MetaGVisualToolBox) we can find this script.

The workflow, utilizes the PhyloPhlAn tool to produce a phylogenetic tree of the MAGs by searching and aligning 400 single-copy phylogenetic markers followed by the R package GGTree to merge and visualize the tree with quality metrics (completeness and contamination) and taxonomy annotations of each MAG.

We can run this into an interactive job:

```bash
/cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/srun.prarameters.Nonode.Account.sh 16 10G normal,bigmem,hugemem 20G nn9987k 04:00:00

```

Then PhyloPhlAn needs the amminoacid translated sequences from our MAGs. Let's use the ones from CompareM2:

```bash
mkdir -p $LOCALSCRATCH/ProteinPredictions && cd $LOCALSCRATCH
rsync -aLhv \
/cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2/CompareM.out.dir/samples/MetaBiningBIO326_25Polished.M*/prokka/*.faa \
ProteinPredictions/
```
>[!Note]
> If you were not able to produce some of the files you can copy from 
<details>

```
sending incremental file list
MetaBiningBIO326_25Polished.MaxBin.out.001.faa
MetaBiningBIO326_25Polished.Metabat2.21.faa
MetaBiningBIO326_25Polished.Metabat2.27.faa
MetaBiningBIO326_25Polished.Metabat2.32.faa
MetaBiningBIO326_25Polished.Metabat2.5.faa

```
</details>

Now we can call the PHYLOGENETICS tool box (Conda enviroment)

```bash
module load Miniconda3/23.10.0-1
eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/PHYLOGENETICS/
```

Then run PhyloPhlAn:

```bash
phylophlan \
-i ProteinPredictions \
-d /cluster/projects/nn9987k/.share/db/phylophlanDataBase/phylophlan \
-t a \
--diversity high \
-f /cluster/projects/nn9987k/.share/db/phylophlanDataBase/supermatrix_aa.cfg \
--verbose \
--nproc $SLURM_CPUS_ON_NODE

```

This will produce the directory ```ProteinPredictions_phylophlan```

```bash
ls ProteinPredictions_phylophlan/
```

<details>

```
ProteinPredictions_concatenated.aln  RAxML_bestTree.ProteinPredictions_refined.tre  RAxML_result.ProteinPredictions_refined.tre
ProteinPredictions_resolved.tre      RAxML_info.ProteinPredictions_refined.tre      tmp
ProteinPredictions.tre               RAxML_log.ProteinPredictions_refined.tre
```

</details>


We can use the the [MetaGVisualToolBox/GenoTaxoTree.R](https://github.com/TheMEMOLab/MetaGVisualToolBox/blob/main/scripts/GenoTaxoTree.R) script to then combine the phylogenetic tree produced by PhyloPhlAn with the Completeness and contamination information:

Let's gather the needed data from COMPAREM2

```bash
rsync -aLhv /cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2/CompareM.out.dir/checkm2/quality_report.tsv .
rsync -aLhv /cluster/projects/nn9987k/$USER/metaG/results/COMPAREM2/CompareM.out.dir/gtdbtk/gtdbtk.bac120.summary.tsv .
```

The GenoTaxoTree.R script needs:

```
Input: CheckM2 tabular results, GTDB-Tk classification tables of Bacteria and Archaea, PhyloPhlAn phylogenetic tree.
Output: A figure of the phylogenetic tree produced by PhyloPhlAn annotated with taxonomy and quality information from GTDB-Tk and CheckM2 by circular heatmaps
```
Like arguments:

```
CHECK <- args[1] #CheckM2 results table
BAC <- args[2]  #GTDBTk Bacterial classification
ARCH <- args[3] #GTDBTk Archeal classification
TREE <- args[4] #Phyloplhlan RAXMl tree
OUT <- args[5] #Output name
```

As we do no have archeas in our experiment lest create a mock file:

```bash
head -1  gtdbtk.bac120.summary.tsv > gtdbtk.ar53.summary.tsv
```

The R script is already installed in our PHYLOGENETICS enviroment so we can run it like:

```bash
Rscript $CONDA_PREFIX/bin/GenoTaxoTree.R quality_report.tsv gtdbtk.bac120.summary.tsv gtdbtk.ar53.summary.tsv ProteinPredictions_phylophlan/ProteinPredictions.tre BIO326Tree
```

Let's copy the PDF to our results:

```bash
 rsync -aPLhv BIO326Tree.pdf /cluster/projects/nn9987k/$USER/metaG/results/Visualization/
```

Also the PhyloPLHAN results:

```bash
rsync -aPLhv --exclude "tmp" ProteinPredictions_phylophlan /cluster/projects/nn9987k/$USER/metaG/results/Visualization/
```

We can now exit the interactive session.

Transfer the pdf to our local PC and explore:

```bash
scp auve@saga.sigma2.no:/cluster/projects/nn9987k/auve/metaG/results/Visualization/*.pdf .
```
>[!Important]
> This command must be run in your local PC not in SAGA!

![PHYL](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/Phylo.png)

### 7.2 distilR to summarize functional KEGG annotations

[distillR](https://github.com/anttonalberdi/distillR) is an R package for distilling functional annotations of bacterial genomes and metagenomes into meaningful quantitative metrics defined as Genome-Inferred Functional Traits (GIFT). The package relies on a curated database of ~500 metabolic pathways and gene clusters (collectivelly refered to as 'bundles') to calculate standardised genome-inferred functional traits using KEGG and Enzyme Commission (EC) identifiers. 

To do this we do not need a HPC we can run all the distillR commands locally in our PC. But we need some data from DRAM (```anntoations.tsv```) and GTDB-Tk (```gtdbtk.bac120.summary.tsv```). We can put this into a common place and rsync it together to our machine like:

```bash
ln -s /cluster/projects/nn9987k/auve/metaG/results/DRAM/DRAM.Results.dir/dram.annotation.dir/annotations.tsv /cluster/projects/nn9987k/auve/metaG/results/Visualization/
ln -s /cluster/projects/nn9987k/auve/metaG/results/COMPAREM2/CompareM.out.dir/gtdbtk/gtdbtk.bac120.summary.tsv /cluster/projects/nn9987k/auve/metaG/results/Visualization/
```

Then in our local PC we can do something like:

```bash
scp auve@saga.sigma2.no:/cluster/projects/nn9987k/auve/metaG/results/Visualization/*.tsv .
```
>[!Important]
> This command must be run in your local PC not in SAGA!



Then we can load these in R following the next code:

<details>

<summary>The following R code will use the gtdbtk and DRAM annotation to combine the taxonomy and functional annotation into a herachical clustering of GIFTs </summary>


```r

library(tidyverse)
library(distillR)
library(viridis)
library(Polychrome)
library(pheatmap)

#Read the annotation from DRAM
Annot <- read_tsv("annotations.tsv") %>%
  rename("GENE"=`...1`) %>%
  rename("Genome"=fasta) %>%
  select(GENE, Genome, ko_id, kegg_hit) %>%
  mutate(
    distill = ifelse(
      grepl("\\[EC:", kegg_hit), # Check if kegg_hit contains [EC:
      paste(ko_id, sub(".*(\\[EC:.*?])", "\\1", kegg_hit), sep = " "), # Combine ko_id and extracted EC number
      ko_id # If no [EC:], keep only ko_id
    )
  ) %>%
  select(Genome, distill) %>%
  drop_na()

#use distillR to obtain the GIFTS

GIFTs <- distill(Annot,GIFT_db,genomecol=1,annotcol=2)

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs,GIFT_db) %>%
  as.data.frame() %>%
  rownames_to_column("MAGs")

#Extract the EC number for some distill vaules use EC numbers and not KO

AnnotEC <- Annot %>%
  mutate(distill=str_replace(distill, ".*?\\s+", ""))


GIFTSEC <- distill(AnnotEC,GIFT_db,genomecol=1,annotcol=2)  

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elementsEC <- to.elements(GIFTSEC,GIFT_db)%>%
  as.data.frame() %>%
  rownames_to_column("MAGs")

GIFTS_ElementsTotal <- GIFTs_elements %>%
    bind_rows(GIFTs_elementsEC) %>%
  pivot_longer(-MAGs, names_to = "Sample", values_to = "Value") %>%
  group_by(MAGs, Sample) %>%
  filter(Value == max(Value)) %>%
  ungroup() %>%
  distinct(MAGs, Sample, .keep_all = TRUE) %>%  # in case max values tie
  pivot_wider(names_from = Sample, values_from = Value) %>%
  select(MAGs, sort(tidyselect::peek_vars()[-1]))
#Extract the functional trais from the database

Elements <- GIFT_db %>%
  select(Code_element,Domain, Function)

#This DF will help to create an annotad column
GIFTSannotCol <- Elements %>%
  select(Code_element,Function) %>%
  distinct() %>%
  column_to_rownames("Code_element")

#This DF will help to create an annotad column
GIFTSannotCol <- Elements %>%
  select(Code_element,Function) %>%
  distinct() %>%
  column_to_rownames("Code_element")

##Taxonomy from GTDBTK

Tax <- read_tsv("gtdbtk.bac120.summary.tsv") %>%
  dplyr::select(user_genome,classification) %>%
  rename("Genome"=user_genome)

#Modify the taxonomy

GenoTaxoInfo <- Tax %>%
  select(Genome,classification) %>%
  separate(classification,sep=";",
           into=c("D",
                  "P",
                  "C",
                  "O",
                  "Fa",
                  "G",
                  "S"),
           remove = FALSE) %>%
  mutate_at(.vars=vars(S),
            .funs = ~
              str_remove(.,
                         pattern = ".*[[:space:]]")) %>%
  mutate_at(.vars=vars(S),
            .funs = ~
              str_remove(.,
                         pattern = "(?<=sp)\\d.*")) %>%
  mutate_at(.vars = vars(D),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*d__")) %>%
  mutate_at(.vars = vars(P),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*p__")) %>%
  mutate_at(.vars = vars(C),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*c__")) %>%
  mutate_at(.vars = vars(O),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*o__")) %>%
  mutate_at(.vars = vars(Fa),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*f__")) %>%
  mutate_at(.vars = vars(G),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*g__")) %>%
  mutate_at(.vars = vars(S),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*s__")) %>%
  mutate_at(.vars = vars(D,P,C,O,Fa,G,S),
            .funs = ~
              str_remove_all(.,
                             pattern = "_..*")) %>%
  mutate_at(.vars = vars(D,P,C,O,Fa,G,S),
            .funs = ~
              str_replace(.,
                          pattern = "^$",
                          replacement = ";")) %>%
  unite("classification",D:P:C:O:Fa:G:S,sep = "_") %>%
  mutate_at(.vars = vars(classification),
            .funs = ~
              str_replace(.,
                          pattern = "_;.*",
                          replacement = ""))%>%
  separate(classification,sep = "_",
           remove = FALSE,
           into = c("Domain",
                    "Phylum",
                    "Class",
                    "Order",
                    "Family",
                    "Genus",
                    "Species")) %>%
  mutate(Genus=coalesce(Genus,Order)) %>%
  dplyr::arrange(Phylum)

#Save this table as RDS useful for other propuses

saveRDS(GenoTaxoInfo,file="GenoTaxoInfo.RDS")

#Tibble to get the tax for heatmap
TaxAnnot <- GenoTaxoInfo %>%
  select(Genome,Order) %>%
  column_to_rownames("Genome")


#List of nice colors  for rows (Taxa) and columns (Funcitions)


ColorAnnot <- list(Order=(TaxAnnot %>%
                            select(Order) %>%
                            distinct() %>%
                            mutate(Color=RColorBrewer::brewer.pal(length(TaxAnnot %>% 
                                                                           select(Order) %>% 
                                                                           distinct() %>% 
                                                                           deframe()),"Set1")) %>%
                            deframe()),
                   Function=(Elements %>%
                               select(Function) %>%
                               distinct() %>%
                               mutate(Color=palette36.colors(n=21)) %>%
                               deframe()))

Color <- rev(viridis(10))

#Generate the Heatmap

GIFTSHeatmap <- GIFTS_ElementsTotal %>% 
  column_to_rownames("MAGs") %>% pheatmap::pheatmap(,
                                   cluster_cols = F,
                                   color = Color,
                                   annotation_col = GIFTSannotCol,
                                   annotation_row = TaxAnnot,
                                   annotation_colors = ColorAnnot,
                                   show_colnames = F,
                                   show_rownames = F,
                                   cellwidth = 2,
                                   cellheight = 18, 
                                   border_color = F)

ggsave(GIFTSHeatmap,file="Gifts.DRAM.GTDBTk.pdf",width = 20,height = 20)

```

</details>


![GIFTS](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/GIFTS.png)

### 7.3 Coverage of MAGs in different MetaG samples:

As we build our MAGs from different set of fastq reads from 3 different technical samples (FastPrep, Vortex and Vortex-SRE), we test if there is difference in the abundance of MAGs accross these samples using the coverage of each genome and relative abundance in each sample using the tool CoverM.

CoverM aims to be a configurable, easy to use and fast DNA read coverage and relative abundance calculator focused on metagenomics applications.

CoverM calculates coverage of genomes/MAGs coverm genome or individual contigs coverm contig. Calculating coverage by read mapping, its input can either be BAM files sorted by reference, or raw reads and reference genomes in various formats.

The following script runs coverM

```bash
cd /cluster/projects/nn9987k/$USER
sbatch /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/7_CoverM.SLURM.sh /cluster/projects/nn9987k/$USER/metaG/results/ChopperBio326_25 /cluster/projects/nn9987k/$USER/metaG/results/DREPLICATION/DEREP_BIO326_25.DREP.70.5.out/dereplicated_genomes fasta /cluster/projects/nn9987k/$USER/metaG/results/COVERM && mkdir -p /cluster/projects/nn9987k/$USER/metaG/results/COVERM
```

This will produce a tsv table: 

```
/cluster/projects/nn9987k/$USER/metaG/results/COVERM/Abundance.coverM.tsv
```

We can copy to our computer and use R to see if there is difference in the abundance of MAGs

```bash
rsync -aLhv auve@saga.sigma2.no:/cluster/projects/nn9987k/auve/metaG/results/COVERM/*.tsv .
```
>[!Important]
> This command must be run in your local PC not in SAGA!

Let's load the data in R and produce a heatmap to see differences in relative abundance:


<details>

<summary>The following R code will summarize the differences in relative abundance of different MAGs in our MetaG samples </summary>

```R
library(tidyverse)
library(pheatmap)

CoverM <- read_tsv("Abundance.coverM.tsv")


#Lets remove the unmapped and keep pretty names:

CoverM <- CoverM %>%
  filter(Genome != "unmapped") %>%
  select(Genome, matches("Relative Abundance")) %>%
  rename_with(~ str_remove(., "\\.chopper.*$"), 
              .cols = -Genome)

#Load the Taxonomy

GenoTaxoInfo <- readRDS("GenoTaxoInfo.RDS")

#Tibble to get the tax for heatmap
TaxAnnot <- GenoTaxoInfo %>%
  select(Genome,Order) %>%
  column_to_rownames("Genome")

#Create an Annotation for the column:

ColAnnot <- CoverM %>%
  select(-Genome) %>%
  colnames() %>%
  enframe(value = "Sample") %>%
  mutate(Method=str_remove(Sample,"_\\d+.*")) %>%
  select(-name) %>% 
  column_to_rownames("Sample")

#List of nice colors  for rows (Taxa) and columns (Funcitions)


ColorAnnot <- list(Order=(TaxAnnot %>%
                            select(Order) %>%
                            distinct() %>%
                            mutate(Color=RColorBrewer::brewer.pal(length(TaxAnnot %>% 
                                                                           select(Order) %>% 
                                                                           distinct() %>% 
                                                                           deframe()),"Set1")) %>%
                            deframe()),
                   Method=(colnames(CoverM %>% 
                                      select(-Genome)) %>% 
                             enframe(value = "Sample") %>%
                             mutate(Method=str_remove(Sample,"_\\d+.*")) %>%
                             select(-name) %>% 
  column_to_rownames("Sample") %>%
    distinct() %>%
    mutate(Color=RColorBrewer::brewer.pal(3,"Dark2")) %>%
    deframe()))

#Create the heatmap:


Color <- rev(inferno(500))

CoverM %>%
  column_to_rownames("Genome") %>% 
  pheatmap::pheatmap(, 
                     cluster_row=F,
                     cellwidth = 10,
                     cellheight = 10,
                     color=Color,
                     annotation_col = ColAnnot,
                     annotation_row = TaxAnnot,
                     annotation_colors = ColorAnnot)

#Normalize the Abundance using log2

CoverMlog2 <- CoverM %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1))) %>%
  column_to_rownames("Genome") %>%
  as.matrix()

#Use this as a heatmap:

CoverMPHM <- CoverMlog2 %>% 
  pheatmap::pheatmap(, 
                     cluster_row=F,
                     cellwidth = 10,
                     cellheight = 10,
                     color=Color,
                     annotation_col = ColAnnot,
                     annotation_row = TaxAnnot,
                     annotation_colors = ColorAnnot,
                    show_colnames = F,border_color = F)

#Save the heatmap

ggsave(CoverMPHM,file="CoverM_heatmap.pdf")

# Pivot the table and set the treatments

CoverM_long <- CoverM %>%
  pivot_longer(-Genome, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Treatment = str_remove(Sample, "_\\d+$"))  #


# Run pairwise wilcox test for each genome

pairwise_wilcox <- CoverM_long %>%
  group_by(Genome) %>%
  summarise(
    test = list(pairwise.wilcox.test(Abundance, Treatment, p.adjust.method = "fdr")),
    .groups = "drop"
  )


# Extract and tidy each pairwise matrix

Wilcox_tidy_results <- pairwise_wilcox %>%
  mutate(p_matrix = map(test, ~ .x$p.value)) %>%
  select(Genome, p_matrix) %>%
  rowwise() %>%
  mutate(pairs = list(as_tibble(as.table(p_matrix), .name_repair = "unique"))) %>%
  unnest(pairs) %>%
  rename(Treatment1 = ...1, Treatment2 = ...2, p_adj = n) %>%
  filter(!is.na(p_adj))


```

</details>

![COVER](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/COVERM.png)