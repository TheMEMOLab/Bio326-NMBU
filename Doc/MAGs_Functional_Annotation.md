# Prokaryotic functional annotation

### So far we were able to obtain indiviual metagenome assembled genomes (MAGs) from the xantan enrichment gut environment. Let's review our workflow to see what is the next step: 
![workflow](https://github.com/avera1988/NMBU-Bio-326/blob/main/images/wrokflowmetagenome.png) 

After binning with [MetaBat](https://bitbucket.org/berkeleylab/metabat/src/master/) we obtained 9 bins. We used [CheckM](https://ecogenomics.github.io/CheckM/) to asses the quality of these MAGs and this was the result after run the "qa" pipeline:

```console
$ singularity exec /cvmfs/singularity.galaxyproject.org/c/h/checkm-genome:1.1.3--py_1 \
checkm \
qa \
$out_path/lineage.ms \
$out_path \
-o 2 > $out_path/ONT_qa_bins
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------
  Bin Id            Marker lineage         # genomes   # markers   # marker sets   Completeness   Contamination   Strain heterogeneity   Genome size (bp)   # ambiguous bases   # scaffolds   
# contigs   N50 (scaffolds)   N50 (contigs)   Mean scaffold length (bp)   Mean contig length (bp)   Longest scaffold (bp)   Longest contig (bp)     GC    GC std (scaffolds > 1kbp)   Coding d
ensity   Translation table   # predicted genes    0     1    2   3   4   5+  
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------
  ONT_bin.7   o__Bacteroidales (UID2657)      160         490           268           99.13            0.75               0.00               6455207                0                1        
    1           6455207          6455207               6455207                    6455207                  6455207                6455207         42.66              0.00                 90.7
0                11                 5055          3    485   2   0   0   0   
  ONT_bin.3   o__Clostridiales (UID1226)      155         278           158           98.67            1.27               0.00               6427786                0                1        
    1           6427786          6427786               6427786                    6427786                  6427786                6427786         50.27              0.00                 87.8
6                11                 5776          5    272   0   1   0   0   
  ONT_bin.1   o__Clostridiales (UID1212)      172         263           149           97.99            2.85              71.43               3054291                0                1        
    1           3054291          3054291               3054291                    3054291                  3054291                3054291         38.05              0.00                 89.6
7                11                 2667          4    252   7   0   0   0   
  ONT_bin.2    g__Bacteroides (UID2691)        33         839           309           81.93            0.46              62.50               5130773                0                33       
    33           211788           211788                155477                     155477                   405112                 405112         42.66              1.64                 89.1
9                11                 4231         155   676   8   0   0   0   
  ONT_bin.8   o__Bacteroidales (UID2621)      198         427           260           72.07            0.00               0.00               3472450                0                25       
    25           198197           198197                138898                     138898                   342879                 342879         45.20              1.71                 90.0
9                11                 2911         102   325   0   0   0   0   
  ONT_bin.5   o__Clostridiales (UID1226)      155         278           158           63.40            0.00               0.00               3278464                0                31       
    31           141330           141330                105756                     105756                   511058                 511058         49.37              1.66                 88.2
6                11                 3108          93   185   0   0   0   0   
  ONT_bin.9      k__Bacteria (UID203)         5449        104            58            8.62            0.00               0.00                382818                0                6        
    6            65540            65540                 63803                      63803                    130590                 130590         42.32              0.99                 91.0
7                11                 297           99    5    0   0   0   0   
  ONT_bin.6          root (UID1)              5656         56            24            0.00            0.00               0.00                753739                0                5        
    5            230855           230855                150747                     150747                   247841                 247841         45.27              0.86                 84.9
7                11                 773           56    0    0   0   0   0   
  ONT_bin.4          root (UID1)              5656         56            24            0.00            0.00               0.00                482281                0                3        
    3            166602           166602                160760                     160760                   166615                 166615         41.07              0.52                 87.0
3                11                 543           56    0    0   0   0   0   
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------

```

Now can start filtering the MAGs. A good criterion is to use the quality and contamination of the MAGs to sort them into *High*, *Medium* and *Low* quality MAGs. 
We can use the following table from [Bowers et al.,](https://www.nature.com/articles/nbt.3893) to classify the MAGs:
![tablemags](https://github.com/avera1988/NMBU-Bio-326/blob/main/images/mags.jpg)

We have all this information from the table generated by checkm, and we can use it to extract these quality score parameters. Now we have 9 MAGs and we can easily pick manually those MAGs ans sorted them for quality and contamination. Let's say >= 70 % completeness and =< 5 % contamination by looking into this table. But what if we had ended up with hundred or thousan of genomes ? Doing a manual selection would not be an easy option. So now let's code a little script to automatically select those *"Good quality MAGs"*:

- First, login to Orion and ask for an interactive node to work, we do not need many resources so 4G ram and 2 CPUs for 2 hrs will be enough:

```console
[bio326-21-0@login ~]$ srun --cpus-per-task 2 --mem=4G --time=02:00:00 --pty bash -i

Welcome to the NMBU Orion compute cluster environment.

You are logged in to a machine that can be used to access your home directory,
edit your scripts, manage your files, and submit jobs to the cluster environment.
Do not run any jobs on this machine, as they might be automatically terminated.

IMPORTANT:
  - Orion introduction: https://orion.nmbu.no/
  - Orion can handle small-scale projects. Need more CPU hours? Please consider
    applying for national infrastructure resources: https://www.sigma2.no/
  - Please, PLEASE do compress your fastq, vcf and other non-compressed files
    using i.e. pigz.

NEWS:
  - 2020-10-08: Orion has been re-built. We are still working out many details.
    Please email us if you miss anything, or notice any issues.

For any Orion related enquiry: orion-support@nmbu.no
PS: We are on Teams: https://bit.ly/orion-teams
[bio326-21-0@cn-16 ~]$
```
- Go the folder where are the final MAGs from MetaBat2 and the checkM results. In my case is in the ```$SCRATCH``` inside a folder named something like ```MetagenomicMAGs```, so go there and check you have the 9 MAGs and the checkM results folder: 

```console
[bio326-21-0@cn-16 ~]$ cd $SCRATCH/MetagenomicMAGS
[bio326-21-0@cn-16 MetagenomicMAGS]$ ls
checkM_results  ONT_bin.1.fa  ONT_bin.2.fa  ONT_bin.3.fa  ONT_bin.4.fa  ONT_bin.5.fa  ONT_bin.6.fa  ONT_bin.7.fa  ONT_bin.8.fa  ONT_bin.9.fa
```
*Remember to change your paths and directories to your user name*

- Now let's run again the ```checkm qa``` pipeline but indicating we want to print a *table separated values* file with the quality scores. How do we obtain that?, let's take a look into the checkm qa options. **We have installed a conda environment with checkM ```/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM```. So we need to first load that environment** 

```console
[bio326-21-0@cn-16 MetagenomicMAGS]$ source /mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/activate.conda.sh
[bio326-21-0@cn-16 MetagenomicMAGS]$ conda activate /net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM) [bio326-21-0@cn-16 MetagenomicMAGS]$ checkm qa --help
usage: checkm qa [-h] [-o {1,2,3,4,5,6,7,8,9}]
                 [--exclude_markers EXCLUDE_MARKERS] [--individual_markers]
                 [--skip_adj_correction] [--skip_pseudogene_correction]
                 [--aai_strain AAI_STRAIN] [-a ALIGNMENT_FILE]
                 [--ignore_thresholds] [-e E_VALUE] [-l LENGTH]
                 [-c COVERAGE_FILE] [-f FILE] [--tab_table] [-t THREADS] [-q]
                 [--tmpdir TMPDIR]
                 marker_file analyze_dir
```

- Assess bins for contamination and completeness.

**The options we need are: --tab_table and -f to save this into a file and not just printed to the standar input. We need the checkM results folder and the marker_files**

