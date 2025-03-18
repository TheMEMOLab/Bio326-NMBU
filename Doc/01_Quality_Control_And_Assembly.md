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
/cluster/projects/nn9987k/BIO326-2025/HPC101/SLURM/srun.prarameters.Nonode.Account.sh 10 20G normal,bigmem,hugemem 120G nn9987k 02:00:00
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
@38d36892-c536-4918-9c33-d0c3e207ce07 runid=b750b83920e5186c0fefcc3ad3d1ff830aec4a68 
CTGGCCAACAGATTCATAGCACAAGCCATTCCAGCATAAGTGGAGAATGGTTTCTTTTTGGCATCAG
+
HFHGIGLJJMJJHVIJPDEEDFKLKLML{{QKKOPJJIPKRML{SNKHJ{HGLFFGFFAAKHNIJIPK
(END) 
```
>[!Tip]
>To exit press the letter Q in your keyboard

Now let's Run NanoPlot ond the rawReads and Chopper at the same time.

>[!Tip]
>Let's use multiple CPUs and GNU parallel to run multiple task at the same time. You can read more on parallel [here](https://bioinformaticsworkbook.org/Appendix/GNUparallel/GNU_parallel_examples.html#gsc.tab=0)

- We need to load a module Parallel

```bash
module load parallel/20240322-GCCcore-13.2.0
```

- Then run both NanoPlot and Chopper at the same time using half of the cpus for one task and other half for the second one:

```
parallel --jobs 2 ::: \
    "gunzip -c D01T6_T.fq.gz|chopper --threads 9 -q 10 -l 1000 > D01T6_T.chopper.fq" \
    "NanoPlot --fastq D01T6_T.fq.gz --N50 --loglength -o D01T6_T.Nanoplot.dir" \
    && NanoPlot --fastq D01T6_T.chopper.fq --N50 --loglength -o D01T6_T.chopper.Nanoplot.dir
```



>[!Warning]
>Chopper will run in 12 min but NanoPlot will take around 1 hr to run so, let's cancel the process by typing ```ctrl+c``` and use the pre-made results:

```bash

rsync -aPLhv /cluster/projects/nn9987k/BIO326-2025/metaG/results/D01T6_T.Nanoplot.dir $LOCALSCRATCH
rsync -aPLhv /cluster/projects/nn9987k/BIO326-2025/metaG/results/D01T6_T.chopper.Nanoplot.dir $LOCALSCRATCH
```

Let's list the results:

```
$tree
```

There are a bunch of files but we really need the stats to compare, for example

```bash
head -8 D01T6_T.Nanoplot.dir/NanoStats.txt
```

```
General summary:
Mean read length:                  5,242.1
Mean read quality:                    17.6
Median read length:                4,404.0
Median read quality:                  18.9
Number of reads:               2,332,038.0
Read length N50:                   5,894.0
STDEV read length:                 2,930.4

```

>[!Tip]
>Remember we can use R to do this:

```bash
cp D01T6_T.Nanoplot.dir/NanoStats.txt $LOCALSCRATCH/D01T6_T_Raw.txt && cp D01T6_T.chopper.Nanoplot.dir/NanoStats.txt $LOCALSCRATCH/D01T6_T_choper.txt

```
R is intalled as conda environment so let's call Miniconda3 and R_env

```bash
conda deactivate
conda activate /cluster/projects/nn9987k/.share/conda_environments/R_env/
```



<details>

And then we can run the following:

<summary> R script</summary>

```R
# Load required libraries
library(tidyverse)

# Function to read and process data from a file
read_summary <- function(file, label) {
  read_lines(file) %>%
    # Split each line into name and value
    str_split_fixed(":\\s+", 2) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    # Name the columns
    rename(Metric = V1, Value = V2) %>%
    # Convert Value column to numeric
    mutate(Value = as.numeric(gsub(",", "", Value)),
           File = label) # Add file label
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure two arguments are provided
if (length(args) != 2) {
  stop("Usage: Rscript plot_metrics.R <file_A> <file_B>")
}

# File paths from arguments
file_a <- args[1]
file_b <- args[2]

# Extract file labels (names without .txt)
label_a <- tools::file_path_sans_ext(basename(file_a))
label_b <- tools::file_path_sans_ext(basename(file_b))

# Read data from files with labels
data_a <- read_summary(file_a, label_a)
data_b <- read_summary(file_b, label_b)

# Combine data from both files
combined_data <- bind_rows(data_a, data_b)

# Filter and scale the required metrics
filtered_data <- combined_data %>%
  filter(Metric %in% c("Mean read length", "Mean read quality",
                       "Median read length", "Number of reads", "Total bases")) %>%
  # Scale metrics
  mutate(Value = case_when(
           Metric == "Mean read length" ~ Value / 1e3,
           Metric == "Median read length" ~ Value / 1e3,
           Metric == "Total bases" ~ Value / 1e9,
           Metric == "Number of reads" ~ Value / 1e6,
           TRUE ~ Value
         ),
         Metric = case_when(
           Metric == "Mean read length" ~ "Mean read length (Thousands)",
           Metric == "Median read length" ~ "Median read length (Thousands)",
           Metric == "Total bases" ~ "Total bases (Billions)",
           Metric == "Number of reads" ~ "Number of reads (Millions)",
           TRUE ~ Metric
         ))

# Check unique File values to debug issues
print(unique(filtered_data$File))

# Plot the data
plot <- ggplot(filtered_data, aes(x = Metric, y = Value, fill = File)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Selected Metrics",
       x = "Metric",
       y = "Value",
       fill = "File") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Use File names directly for colors
  scale_fill_manual(values = setNames(c("#1f78b4", "#33a02c"), c(label_a, label_b)))

# Save the plot to a file
output_file <- "comparison_NanoStats.pdf"
ggsave(output_file, plot, width = 10, height = 6)
cat("Plot saved to", output_file, "\n")


```

</details>

```bash
Rscript /cluster/projects/nn9987k/BIO326-2025/metaG/scripts/plotNanoStats.r $LOCALSCRATCH/D01T6_T_Raw.txt $LOCALSCRATCH/D01T6_T_choper.txt
```
This will produce a pdf file:

```
Plot saved to comparison_NanoStats.pdf
```
![NanoStats](https://github.com/TheMEMOLab/Bin420-Bioinformatics-for-Functional-Meta-Omics/blob/main/img/NanoStats.PNG)

Let's copy all the results to to our folder:

```bash
rsync -aPLhv comparison_NanoStats.pdf /cluster/projects/nn9987k/$USER/metaG/results/
rsync -aPLhv D01T6_T.Nanoplot.dir /cluster/projects/nn9987k/$USER/metaG/results/
rsync -aPLhv D01T6_T.chopper.Nanoplot.dir /cluster/projects/nn9987k/$USER/metaG/results/
```

And of course our Chopper results.
>[!Important]
>Chopper proudces a large fastq un-compressed file so before copy we need to compress it to save space. Remember always compress your files with tools like pigz:

```bash
pigz -p 10 D01T6_T.chopper.fq && rsync -aPLhv D01T6_T.chopper.fq.gz /cluster/projects/nn9987k/$USER/metaG/results/Chopper/
$tree /cluster/projects/nn9987k/$USER/metaG/results/Chopper/
```

We have done, now we need to exit our room and let the SLURM manager to clean,

```bash
exit
```

>[!Important]
>Always type ```exit``` command when you finish working.

