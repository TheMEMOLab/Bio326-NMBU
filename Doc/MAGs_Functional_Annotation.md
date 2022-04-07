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

### Assess bins for contamination and completeness.

**The options we need are: --tab_table and -f to save this into a file and not just printed to the standar input. We need the checkM results folder and the marker_files**
```bash
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM) [bio326-21-0@cn-16 MetagenomicMAGS]$ checkm \
 qa \
 checkM_results/lineage.ms \
 -t $SLURM_CPUS_ON_NODE \
 checkM_results \
 -o 2 \
 --tab_table \
-f ONT_qa_bins.tsv
[2021-04-22 15:56:57] INFO: CheckM v1.1.3
[2021-04-22 15:56:57] INFO: checkm qa checkM_results/lineage.ms -t 2 checkM_results -o 2 --tab_table -f ONT_qa_bins.tsv
[2021-04-22 15:56:57] INFO: [CheckM - qa] Tabulating genome statistics.
[2021-04-22 15:56:57] INFO: Calculating AAI between multi-copy marker genes.
[2021-04-22 15:56:57] INFO: Reading HMM info from file.
[2021-04-22 15:56:57] INFO: Parsing HMM hits to marker genes:
    Finished parsing hits for 9 of 9 (100.00%) bins.
[2021-04-22 15:57:00] INFO: QA information written to: ONT_qa_bins.tsv
[2021-04-22 15:57:00] INFO: { Current stage: 0:00:02.594 || Total: 0:00:02.594 }

```

This will provide us with a *ONT_qa_bins.tsv* file where all the qa results are storage. Let's take a look:

```bash
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM) [bio326-21-0@cn-16 MetagenomicMAGS]$ more ONT_qa_bins.tsv 
Bin Id	Marker lineage	# genomes	# markers	# marker sets	Completeness	Contamination	Strain heterogeneity	Genome size (bp)	# ambiguous bases	# scaffolds	# contigs	N50 (scaffolds)	N50 (contigs)	Mean scaffold length (bp)	Mean contig length (bp)	Longest scaffold (bp)	Longest contig (bp)	GC	GC std (scaffolds > 1kbp)	Coding density	Translation table	# predicted genes	0	1	2	3	4	5+
ONT_bin.1	o__Clostridiales (UID1212)	172	263	149	97.99	2.85	71.43	3054291	0	1	1	3054291	3054291	3054291	3054291	3054291	3054291	38.0	0.00	89.67	11	2667	4	252	7	0	0	0
ONT_bin.2	g__Bacteroides (UID2691)	33	839	309	81.93	0.46	62.50	5130773	0	33	33	211788	211788	155477	155477	405112	405112	42.7	1.64	89.19	11	4231	155	676	8	0	0	0
ONT_bin.3	o__Clostridiales (UID1226)	155	278	158	98.67	1.27	0.00	6427786	0	1	1	6427786	6427786	6427786	6427786	6427786	6427786	50.3	0.00	87.86	11	5776	5	272	0	1	0	0
ONT_bin.4	root (UID1)	5656	56	24	0.00	0.00	0.00	482281	0	3	3	166602	166602	160760	160760	166615	166615	41.1	0.52	87.03	11	543	56	0	0	0	0	0
ONT_bin.5	o__Clostridiales (UID1226)	155	278	158	63.40	0.00	0.00	3278464	0	31	31	141330	141330	105756	105756	511058	511058	49.4	1.66	88.26	11	3108	93	185	0	0	0	0
ONT_bin.6	root (UID1)	5656	56	24	0.00	0.00	0.00	753739	0	5	5	230855	230855	150747	150747	247841	247841	45.3	0.86	84.97	11	773	56	0	0	0	0	0
ONT_bin.7	o__Bacteroidales (UID2657)	160	490	268	99.13	0.75	0.00	6455207	0	1	1	6455207	6455207	6455207	6455207	6455207	6455207	42.7	0.00	90.70	11	5055	3	485	2	0	0	0
ONT_bin.8	o__Bacteroidales (UID2621)	198	427	260	72.07	0.00	0.00	3472450	0	25	25	198197	198197	138898	138898	342879	342879	45.2	1.71	90.09	11	2911	102	325	0	0	0	0
ONT_bin.9	k__Bacteria (UID203)	5449	104	58	8.62	0.00	0.00	382818	0	6	6	65540	65540	63803	63803	130590	130590	42.3	0.99	91.0711	297	99	5	0	0	0	0

```

- So for filtering we need to select all genomes that have a *Completeness*(colum6) >= 70 and *Contamination	Strain* (colum 7). For this conditionals loops, we can use the AWK that is a programming language for data extraction and reporting tool. The goal of this course is not to learn AWK so just let's talk about the basics: 
* $ are references to colums (e.g. $6 meand colum 6)
* -F Command line option to specify input field delimiter (e.g. -F "\t" means the text is separated by tabs)
* awk '/pattern/ {action}' file↵Execute action for matched pattern 'pattern' on file 'file' (e.g awk -F "\t" '{if($6 >= 70 && $7 <= 5) print $1"\t"$6"\t"$7}' ONT_qa_bins.tsv means if colum 6 is greather than 70 and colum 7 is lower than 5 print: col.1 (ID), col.6(completeness) and col.7(Contamination), all separated by tabs ("\t")... 
* This [cheatsheet](https://www.shortcutfoo.com/app/dojos/awk/cheatsheet) is a useful resource if you are interested in learn a bit more. 

After this, let's filtering:

```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM) [bio326-21-0@cn-16 MetagenomicMAGS]$ awk -F "\t" '{if($6 >= 70 && $7 <= 5) print $1"\t"$6"\t"$7}' ONT_qa_bins.tsv
ONT_bin.1	97.99	2.85
ONT_bin.2	81.93	0.46
ONT_bin.3	98.67	1.27
ONT_bin.7	99.13	0.75
ONT_bin.8	72.07	0.00
```

After printing this we notice that five bins matches our quality filter condition. We then can move these MAGs to a new folder, let's named GoodQualityMAGs

```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM) [bio326-21-0@cn-16 MetagenomicMAGS]$ mkdir GoodQualityMAGs
```

And then whe need to move all these five MAGs to that folder. We can do it manually, but we are bioinformaticians, so let's use the computer to move this. For doing that we need to read the first colum (genome ID) and copy those genomes to the new folder GoodQualityMAGs. 

The best way to do this is by using a *while* loop:
* A basic example of while loop for reading "lines" in a file and print those lines (echo) to the standar ouput (screen) is: 
``` bash
cat file.txt | while read line; do
  echo $line
done
```

* In our specific case, the while loop needs to read (using the command *read*) the first colum after awk and copy (cp) that line (each line in the colum 1) to the folder (GoodQualityMAGs). **So far our awk prints 3 colums but for the loop we only need the column 1 ($1) so we need to modify the awk first and then apply the loop. Remember that the name of the file is the name in the column 1 plus the extension .fa so we also need to indicate this in the loop:**


```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/checkM) [bio326-21-0@cn-16 MetagenomicMAGS]$ awk -F "\t" '{if($6 >= 70 && $7 <= 5) print $1}' ONT_qa_bins.tsv| while read line; do cp $line.fa GoodQualityMAGs/;done
[bio326-21-0@cn-16 MetagenomicMAGS]$ cd GoodQualityMAGs/
[bio326-21-0@cn-16 GoodQualityMAGs]$ ls
ONT_bin.1.fa  ONT_bin.2.fa  ONT_bin.3.fa  ONT_bin.7.fa  ONT_bin.8.fa
```

6. **By applying this loop we were able to get all the genomes in the same folder at once. Try the loop if no don't panic, you can copy these genomes one by one using the normal cp command 😊**

We can exit now the interactive job:

```console
[bio326-21-0@cn-16 GoodQualityMAGs]$ exit
exit
```

Now that we have this we can start the annotation using DRAM.

## DRAM: Distilled and Refined Annotation of Metabolism

"[DRAM](https://github.com/shafferm/DRAM#dram) (Distilled and Refined Annotation of Metabolism) is a tool for annotating metagenomic assembled genomes and VirSorter identified viral contigs. DRAM annotates MAGs and viral contigs using KEGG (if provided by the user), UniRef90, PFAM, dbCAN, RefSeq viral, VOGDB and the MEROPS peptidase database as well as custom user databases..."
 ![dramaanot](https://github.com/avera1988/NMBU-Bio-326/blob/main/images/DRAM.jpg)


Until now we have the MAGs, the checkM results and last [session](https://github.com/liveha/NMBU-BIO326/blob/main/Binning_metaBAT.md) we learned how to use [GTDBTK](https://github.com/Ecogenomics/GTDBTk) for taxonomic classification. But we can run GTDBTK again using only the "GoodQuality" MAGs now. For this we need the script gtdbk.classifywf.SLURM.sh 

```bash
#!/bin/bash
#########################################################################
#	SLURM scrip for running gtdbtk on Orion cluster
#		Dependencies: gtdbk conda environment
#					  gtdbk.release95 db
#					  fasta files of genomes (MAGs)
#
#	It copies the database and the MAGs to a local disk on any Orion's node 
#	and runs locally. At the end it copies all the results to the PEP
#	or any other path in the cluster into a MAGs_gtdbk.dir folder.
#
#	to run: 
# sbatch gtdbk.sh path_to_MAGs_folder fasta_files_extension
# eg: sbatch gtdbk.sh /fs-1/PEPomics01/auve/Metasalmon/gtdbktest fna
#
# Author: Arturo Vera
# Jan 2021
#########################################################################

###############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=gtdbk_classifywf
#
## Wall time limit:
#SBATCH --time=08:00:00
#
## Other parameters:
#SBATCH --cpus-per-task 12
#SBATCH --mem=150G
#SBATCH	--nodes=1
#SBATCH --partition=hugemem
###Basic usage help for this script#######

print_usage() {
        echo "Usage: sbatch $0 path_to_MAGs fasta_files_extension"
        echo "eg: sbatch $0 /net/fs-1/PEPomics01/auve/Metasalmon/gtdbktest fna"
}

if [ $# -le 1 ]
        then
                print_usage
                exit 1
        fi


###############Main SCRIPT###################################

## Set up job environment:

module --quiet purge  # Reset the modules to the system default
module load Miniconda3

##Declaring variables

magsdir=$1
ext=$2

####Do some work:########

## For debuggin

echo "Hello " $USER 
echo "my submit directory is"
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

cd $TMPDIR

#Check if $USER exists in $TMPDIR

if [[ -d $USER ]]
then
        echo "$USER exists on $TMPDIR"
else
        mkdir $USER
fi


echo "copying files to $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID"

cd $USER
mkdir tmpDir_of.$SLURM_JOB_ID
cd tmpDir_of.$SLURM_JOB_ID

##Activate conda environments

export PS1=\$
source activate /mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/GTDBTK

#Copy the MAGs to the $TMPDIR

echo "copying MAGs to: " $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID
cd $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID
cp $magsdir/*$ext .   #the extension is the fasta file ext e.g. fna

#Create a MAGs dir and move there the MAGs
mkdir MAGs
mv *.$ext MAGs/

#####GTDBTK classify_wf pipeline################

time gtdbtk classify_wf \
--genome_dir MAGs \
--out_dir MAGs_gtdbk.dir \
-x $ext \
--cpus $SLURM_CPUS_ON_NODE \

cd $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID

############ Moving results to PEP partition or anywhere the main script was submitted ################

echo "moving results to" $SLURM_SUBMIT_DIR

cd $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID

time cp -r MAGs_gtdbk.dir $SLURM_SUBMIT_DIR

echo "gtdbk results are in" $SLURM_SUBMIT_DIR/MAGs_gtdbk.dir

####removing tmp dir. Remember to do this for not filling the HDD in the node############

cd $TMPDIR/$USER/
rm -r tmpDir_of.$SLURM_JOB_ID

echo "I've done at"
date
```
As always you can copy and paste this script from here or from ```/mnt/SCRATCH/bio326-21/MetaGenomeBinning/gtdbk.classifywf.SLURM.sh```. So if you have not produced the GTDBTK result copy the gtdbtk SLURM script and run it as follow:

```console
[bio326-21-0@login MetagenomicMAGS]$ sbatch gtdbk.classifywf.SLURM.sh /mnt/SCRATCH/bio326-21-0/MetagenomicMAGS/GoodQualityMAGs/GoodQualityMAGs fa
```

This will produce the ```MAGs_gtdbk.dir```. *GTDBTK takes >~ 45 min to analyze the MAGs, you either wait or can copy the MAGS_gtdbtk.dir to your folder by  ```cp -r /mnt/SCRATCH/bio326-21-0/MetagenomicMAGS/MAGs_gtdbk.dir .```*

Let's take a look on the Results:

```console
[bio326-21-0@login MetagenomicMAGS]$ cd MAGs_gtdbk.dir/
[bio326-21-0@login MAGs_gtdbk.dir]$ ls
align     gtdbtk.ar122.markers_summary.tsv  gtdbtk.bac120.filtered.tsv         gtdbtk.bac120.msa.fasta    gtdbtk.bac120.user_msa.fasta  gtdbtk.translation_table_summary.tsv  identify
classify  gtdbtk.bac120.classify.tree       gtdbtk.bac120.markers_summary.tsv  gtdbtk.bac120.summary.tsv  gtdbtk.log                    gtdbtk.warnings.log
```

As we can notice GTDBTK prouduces a lot of results. The file containing the final results is ```gtdbtk.bac120.summary.tsv```. Take a look on this:

```console
[bio326-21-0@login MAGs_gtdbk.dir]$ cat gtdbtk.bac120.summary.tsv
user_genome	classification	fastani_reference	fastani_reference_radius	fastani_taxonomy	fastani_ani	fastani_af	closest_placement_reference	closest_placement_radius	closest_placement_taxonomy	closest_placement_ani	closest_placement_af	pplacer_taxonomy	classification_method	note	other_related_references(genome_id,species_name,radius,ANI,AF)	msa_percent	translation_table	red_value	warnings
ONT_bin.1	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__CAG-41 sp900066215	GCF_003460745.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__CAG-41 sp900066215	98.76	0.81	GCF_003460745.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__CAG-41 sp900066215	98.76	0.81	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_001941225.1, s__CAG-41 sp001941225, 95.0, 78.23, 0.38	95.38	11	N/A	N/A
ONT_bin.2	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides thetaiotaomicron	GCF_000011065.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides thetaiotaomicron	97.92	0.8	GCF_000011065.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides thetaiotaomicron	97.92	0.8	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCF_900106755.1, s__Bacteroides faecis, 95.0, 89.07, 0.68; GCF_001688725.2, s__Bacteroides caecimuris, 95.0, 82.26, 0.39; GCF_002222615.2, s__Bacteroides caccae, 95.0, 82.05, 0.44; GCA_007197895.1, s__Bacteroides sp900066265, 95.75, 81.93, 0.36; GCF_000613385.1, s__Bacteroides acidifaciens, 95.0, 81.93, 0.39; GCF_003463205.1, s__Bacteroides sp003463205, 95.0, 81.75, 0.43; GCA_000210075.1, s__Bacteroides xylanisolvens, 95.0, 81.71, 0.43; GCA_007097645.1, s__Bacteroides sp007097645, 95.0, 81.61, 0.42; GCF_001314995.1, s__Bacteroides ovatus, 95.0, 81.52, 0.45; GCF_900130125.1, s__Bacteroides congonensis, 95.0, 81.51, 0.43; GCF_000156195.1, s__Bacteroides finegoldii, 95.75, 81.49, 0.36; GCA_900556625.1, s__Bacteroides sp900556625, 95.0, 81.18, 0.43; GCA_900547205.1, s__Bacteroides sp900547205, 95.0, 81.01, 0.38; GCF_900155865.1, s__Bacteroides bouchesdurhonensis, 95.0, 80.43, 0.32; GCF_003865075.1, s__Bacteroides sp003865075, 95.0, 80.31, 0.31; GCF_000614145.1, s__Bacteroides faecichinchillae, 95.0, 80.04, 0.34; GCF_000315485.1, s__Bacteroides oleiciplenus, 95.0, 79.79, 0.18; GCA_900557355.1, s__Bacteroides sp900557355, 95.0, 79.62, 0.13; GCF_003464595.1, s__Bacteroides intestinalis_A, 95.0, 79.57, 0.18; GCF_000172175.1, s__Bacteroides intestinalis, 95.0, 79.03, 0.18; GCF_000381365.1, s__Bacteroides salyersiae, 95.0, 78.99, 0.2; GCF_000154525.1, s__Bacteroides stercoris, 95.0, 78.98, 0.18; GCF_002849695.1, s__Bacteroides fragilis_A, 95.0, 78.88, 0.18; GCF_000613465.1, s__Bacteroides nordii, 95.0, 78.88, 0.19; GCF_000154205.1, s__Bacteroides uniformis, 95.0, 78.86, 0.16; GCF_003438615.1, s__Bacteroides sp003545565, 95.0, 78.85, 0.19; GCF_000614165.1, s__Bacteroides stercorirosoris, 95.0, 78.84, 0.18; GCF_000025985.1, s__Bacteroides fragilis, 95.0, 78.79, 0.19; GCF_004793475.1, s__Bacteroides sp002491635, 95.0, 78.73, 0.16; GCF_000158035.1, s__Bacteroides cellulosilyticus, 95.0, 78.71, 0.17; GCF_000155815.1, s__Bacteroides eggerthii, 95.0, 78.54, 0.17; GCF_000513195.1, s__Bacteroides timonensis, 95.0, 78.52, 0.17; GCF_900129655.1, s__Bacteroides clarus, 95.0, 78.51, 0.19; GCF_000614125.1, s__Bacteroides rodentium, 95.0, 78.49, 0.15; GCF_000517545.1, s__Bacteroides reticulotermitis, 95.0, 78.47, 0.22; GCF_900241005.1, s__Bacteroides cutis, 95.0, 78.47, 0.18; GCF_000186225.1, s__Bacteroides helcogenes, 95.0, 78.46, 0.16; GCA_000511775.1, s__Bacteroides pyogenes_A, 95.0, 78.43, 0.23; GCF_000374365.1, s__Bacteroides gallinarum, 95.0, 78.42, 0.13; GCF_000195635.1, s__Bacteroides fluxus, 95.0, 78.41, 0.17; GCF_000499785.1, s__Bacteroides neonati, 95.0, 78.41, 0.13; GCA_900555635.1, s__Bacteroides sp900555635, 95.0, 78.41, 0.15; GCA_900552405.1, s__Bacteroides sp900552405, 95.0, 78.38, 0.15; GCF_900108345.1, s__Bacteroides ndongoniae, 95.0, 78.35, 0.13; GCF_000428105.1, s__Bacteroides pyogenes, 95.0, 78.24, 0.23; GCA_900556215.1, s__Bacteroides sp900556215, 95.0, 78.12, 0.16; GCF_002998435.1, s__Bacteroides zoogleoformans, 95.0, 78.1, 0.14; GCF_900130135.1, s__Bacteroides togonis, 95.0, 78.05, 0.12; GCF_004342845.1, s__Bacteroides heparinolyticus, 95.0, 77.97, 0.14; GCF_900128475.1, s__Bacteroides massiliensis, 95.0, 77.94, 0.12; GCA_002471195.1, s__Bacteroides sp002471195, 95.0, 77.94, 0.08; GCA_002471185.1, s__Bacteroides sp002471185, 95.0, 77.92, 0.11; GCF_900128905.1, s__Bacteroides luti, 95.0, 77.47, 0.08; GCF_000428125.1, s__Bacteroides graminisolvens, 95.0, 77.46, 0.1; GCF_002160055.1, s__Bacteroides sp002160055, 95.0, 77.45, 0.08; GCA_002293435.1, s__Bacteroides sp002293435, 95.0, 77.41, 0.14; GCF_900104585.1, s__Bacteroides ihuae, 95.0, 77.36, 0.1; GCA_002307035.1, s__Bacteroides sp002307035, 95.0, 77.09, 0.1; GCA_900553815.1, s__Bacteroides sp900553815, 95.0, 76.96, 0.11	69.72	11	N/A	N/A
ONT_bin.3	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Enterocloster;s__Enterocloster sp000155435	GCF_000155435.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Enterocloster;s__Enterocloster sp000155435	98.63	0.83	GCF_000155435.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Enterocloster;s__Enterocloster sp000155435	98.63	0.83	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Enterocloster;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCF_003434055.1, s__Enterocloster aldenensis, 95.0, 90.97, 0.78; GCF_000233455.1, s__Enterocloster citroniae, 95.0, 81.66, 0.47; GCF_002234575.2, s__Enterocloster bolteae, 95.0, 80.93, 0.36; GCF_000424325.1, s__Enterocloster clostridioformis_A, 95.0, 80.76, 0.42; GCF_005845215.1, s__Enterocloster sp005845215, 95.0, 80.6, 0.37; GCF_900113155.1, s__Enterocloster clostridioformis, 95.0, 80.44, 0.39; GCA_001304855.1, s__Enterocloster sp001304855, 95.0, 78.85, 0.32; GCF_001517625.2, s__Enterocloster sp001517625, 95.0, 78.58, 0.22; GCF_003473545.1, s__Enterocloster sp000431375, 95.0, 78.48, 0.17; GCF_003024655.1, s__Enterocloster lavalensis, 95.0, 78.32, 0.21; GCA_900555045.1, s__Enterocloster sp900555045, 95.0, 78.3, 0.25; GCA_900543885.1, s__Enterocloster sp900543885, 95.0, 78.24, 0.15; GCF_000158075.1, s__Enterocloster asparagiformis, 95.0, 78.21, 0.22; GCA_900549235.1, s__Enterocloster sp900549235, 95.0, 78.2, 0.29; GCA_900551225.1, s__Enterocloster sp900551225, 95.0, 78.19, 0.25; GCA_900538485.1, s__Enterocloster sp900538485, 95.0, 78.15, 0.22; GCA_900540675.1, s__Enterocloster sp900540675, 95.0, 78.12, 0.3; GCA_900541315.1, s__Enterocloster sp900541315, 95.0, 78.08, 0.15; GCA_900555905.1, s__Enterocloster sp900555905, 95.0, 77.93, 0.21; GCA_900547035.1, s__Enterocloster sp900547035, 95.0, 77.3, 0.17	90.18	11	N/A	N/A
ONT_bin.7	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides intestinalis_A	GCF_003464595.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides intestinalis_A	97.03	0.79	GCF_003464595.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides intestinalis_A	97.03	0.79	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_900556215.1, s__Bacteroides sp900556215, 95.0, 94.74, 0.66; GCF_000172175.1, s__Bacteroides intestinalis, 95.0, 94.42, 0.75; GCA_900557355.1, s__Bacteroides sp900557355, 95.0, 92.47, 0.43; GCF_000158035.1, s__Bacteroides cellulosilyticus, 95.0, 89.86, 0.65; GCF_000513195.1, s__Bacteroides timonensis, 95.0, 88.61, 0.6; GCA_900552405.1, s__Bacteroides sp900552405, 95.0, 88.13, 0.6; GCF_000614165.1, s__Bacteroides stercorirosoris, 95.0, 83.88, 0.57; GCF_000315485.1, s__Bacteroides oleiciplenus, 95.0, 83.79, 0.64; GCF_003438615.1, s__Bacteroides sp003545565, 95.0, 81.16, 0.39; GCF_000195635.1, s__Bacteroides fluxus, 95.0, 80.62, 0.37; GCF_003865075.1, s__Bacteroides sp003865075, 95.0, 80.46, 0.17; GCA_007197895.1, s__Bacteroides sp900066265, 95.75, 80.04, 0.2; GCA_002293435.1, s__Bacteroides sp002293435, 95.0, 79.66, 0.44; GCF_000154205.1, s__Bacteroides uniformis, 95.0, 79.4, 0.33; GCF_000154525.1, s__Bacteroides stercoris, 95.0, 79.29, 0.34; GCF_000614125.1, s__Bacteroides rodentium, 95.0, 79.24, 0.32; GCF_900129655.1, s__Bacteroides clarus, 95.0, 79.23, 0.36; GCF_000156195.1, s__Bacteroides finegoldii, 95.75, 79.21, 0.2; GCF_000155815.1, s__Bacteroides eggerthii, 95.0, 79.2, 0.33; GCF_004793475.1, s__Bacteroides sp002491635, 95.0, 79.18, 0.32; GCF_900241005.1, s__Bacteroides cutis, 95.0, 79.13, 0.36; GCF_000381365.1, s__Bacteroides salyersiae, 95.0, 79.09, 0.19; GCF_000186225.1, s__Bacteroides helcogenes, 95.0, 78.86, 0.3; GCF_002222615.2, s__Bacteroides caccae, 95.0, 78.83, 0.22; GCF_000374365.1, s__Bacteroides gallinarum, 95.0, 78.82, 0.26; GCF_000011065.1, s__Bacteroides thetaiotaomicron, 95.0, 78.8, 0.19; GCF_900130135.1, s__Bacteroides togonis, 95.0, 78.79, 0.18; GCF_000613385.1, s__Bacteroides acidifaciens, 95.0, 78.76, 0.21; GCF_001688725.2, s__Bacteroides caecimuris, 95.0, 78.75, 0.21; GCF_002849695.1, s__Bacteroides fragilis_A, 95.0, 78.74, 0.19; GCF_900108345.1, s__Bacteroides ndongoniae, 95.0, 78.55, 0.19; GCA_000210075.1, s__Bacteroides xylanisolvens, 95.0, 78.53, 0.2; GCF_003463205.1, s__Bacteroides sp003463205, 95.0, 78.51, 0.17; GCF_000025985.1, s__Bacteroides fragilis, 95.0, 78.5, 0.19; GCF_900155865.1, s__Bacteroides bouchesdurhonensis, 95.0, 78.46, 0.16; GCA_007097645.1, s__Bacteroides sp007097645, 95.0, 78.45, 0.19; GCF_900130125.1, s__Bacteroides congonensis, 95.0, 78.44, 0.17; GCF_001314995.1, s__Bacteroides ovatus, 95.0, 78.42, 0.18; GCF_000613465.1, s__Bacteroides nordii, 95.0, 78.38, 0.2; GCA_900555635.1, s__Bacteroides sp900555635, 95.0, 78.33, 0.27; GCF_900106755.1, s__Bacteroides faecis, 95.0, 78.25, 0.2; GCF_000614145.1, s__Bacteroides faecichinchillae, 95.0, 78.22, 0.18; GCF_002998435.1, s__Bacteroides zoogleoformans, 95.0, 78.22, 0.26; GCF_004342845.1, s__Bacteroides heparinolyticus, 95.0, 78.2, 0.26; GCA_900547205.1, s__Bacteroides sp900547205, 95.0, 78.08, 0.19; GCA_000511775.1, s__Bacteroides pyogenes_A, 95.0, 77.92, 0.14; GCF_000499785.1, s__Bacteroides neonati, 95.0, 77.9, 0.13; GCF_002160055.1, s__Bacteroides sp002160055, 95.0, 77.8, 0.13; GCF_900128475.1, s__Bacteroides massiliensis, 95.0, 77.79, 0.17; GCF_000428105.1, s__Bacteroides pyogenes, 95.0, 77.74, 0.15; GCA_900556625.1, s__Bacteroides sp900556625, 95.0, 77.55, 0.18; GCA_002471185.1, s__Bacteroides sp002471185, 95.0, 77.53, 0.12; GCF_000517545.1, s__Bacteroides reticulotermitis, 95.0, 77.52, 0.13; GCA_002471195.1, s__Bacteroides sp002471195, 95.0, 77.46, 0.1; GCA_900553815.1, s__Bacteroides sp900553815, 95.0, 77.41, 0.16; GCF_900128905.1, s__Bacteroides luti, 95.0, 77.23, 0.1; GCF_900104585.1, s__Bacteroides ihuae, 95.0, 77.19, 0.11; GCF_000428125.1, s__Bacteroides graminisolvens, 95.0, 77.04, 0.12; GCA_002307035.1, s__Bacteroides sp002307035, 95.0, 76.99, 0.12	97.28	11	N/A	N/A
ONT_bin.8	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis	GCF_000012845.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis	98.42	0.87	GCF_000012845.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis	98.42	0.87	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCF_004793765.1, s__Parabacteroides distasonis_A, 95.0, 93.66, 0.62; GCF_900186615.1, s__Parabacteroides bouchesdurhonensis, 95.0, 80.12, 0.2; GCF_000154105.1, s__Parabacteroides merdae, 95.0, 79.7, 0.22; GCF_000969835.1, s__Parabacteroides goldsteinii, 95.0, 79.61, 0.24; GCF_000156495.1, s__Parabacteroides johnsonii, 95.0, 79.51, 0.23; GCF_003479145.1, s__Parabacteroides sp003479145, 95.0, 79.4, 0.23; GCF_000969825.1, s__Parabacteroides gordonii, 95.0, 79.38, 0.24; GCA_900549585.1, s__Parabacteroides sp900549585, 95.0, 79.28, 0.09; GCF_900128505.1, s__Parabacteroides timonensis, 95.0, 79.24, 0.24; GCF_003480915.1, s__Parabacteroides sp003480915, 95.0, 79.23, 0.23; GCF_003363715.1, s__Parabacteroides acidifaciens, 95.0, 79.21, 0.21; GCF_900155425.1, s__Parabacteroides sp900155425, 95.0, 79.09, 0.22; GCF_900108035.1, s__Parabacteroides chinchillae, 95.0, 78.83, 0.18; GCF_003473295.1, s__Parabacteroides sp003473295, 95.0, 78.63, 0.18; GCA_900540715.1, s__Parabacteroides sp900540715, 95.0, 78.32, 0.21; GCA_900541965.1, s__Parabacteroides sp900541965, 95.0, 78.18, 0.17; GCF_002159645.1, s__Parabacteroides sp002159645, 95.0, 77.94, 0.08; GCA_000436495.1, s__Parabacteroides sp000436495, 95.0, 77.69, 0.1; GCA_004562445.1, s__Parabacteroides sp004562445, 95.0, 77.68, 0.1; GCA_900552465.1, s__Parabacteroides sp900552465, 95.0, 77.55, 0.2; GCA_900552415.1, s__Parabacteroides sp900552415, 95.0, 77.39, 0.17; GCA_900548175.1, s__Parabacteroides sp900548175, 95.0, 77.24, 0.09; GCA_900547435.1, s__Parabacteroides sp900547435, 95.0, 76.84, 0.09	77.78	11	N/A	N/A

```

This is a **HUGE** file and is quite confusing so just let's take a closer look into the header:

```console
[bio326-21-0@login MAGs_gtdbk.dir]$ head -1 gtdbtk.bac120.summary.tsv
user_genome	classification	fastani_reference	fastani_reference_radius	fastani_taxonomy	fastani_ani	fastani_af	closest_placement_reference	closest_placement_radius	closest_placement_taxonomy	closest_placement_ani	closest_placement_af	pplacer_taxonomy	classification_method	note	other_related_references(genome_id,species_name,radius,ANI,AF)	msa_percent	translation_table	red_value	warnings
```

The columns describing the genome id, taxonomy classification and ANI are 1, 2 and 6, we can pars this to have a more tidy useful table by "cutting" the columns 1, 2 and 6 using the *cut* command:

```console
[bio326-21-0@login MAGs_gtdbk.dir]$ cat gtdbtk.bac120.summary.tsv |cut -f 1,2,6
user_genome	classification	fastani_ani
ONT_bin.1	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__CAG-41 sp900066215	98.76
ONT_bin.2	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides thetaiotaomicron	97.92
ONT_bin.3	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Enterocloster;s__Enterocloster sp000155435	98.63
ONT_bin.7	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides intestinalis_A	97.03
ONT_bin.8	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis	98.42
```
*Most of the MAGs have the closest reference species in the database with ANI values > 97 %, so we can say that these MAGs are the same bacterial species as those previously described in the database. And all of them seem to be human gut-associated bacteria...*

### We have now all the imputs for DRAM. Let's see what does DRAM need to run. 

1. DRAM is installed as a conda environment we need to activate this first:

```bash
[bio326-21-0@login MAGs_gtdbk.dir]$ module load Miniconda3 
[bio326-21-0@login MAGs_gtdbk.dir]$ source activate /net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login MAGs_gtdbk.dir]$
```

2. Print the help of DRAM

```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login MAGs_gtdbk.dir]$ DRAM.py --help
usage: DRAM.py [-h] {annotate,annotate_genes,distill,strainer,neighborhoods} ...

positional arguments:
  {annotate,annotate_genes,distill,strainer,neighborhoods}
    annotate            Annotate genomes/contigs/bins/MAGs
    annotate_genes      Annotate already called genes, limited functionality compared to annotate
    distill             Summarize metabolic content of annotated genomes
    strainer            Strain annotations down to genes of interest
    neighborhoods       Find neighborhoods around genes of interest

optional arguments:
  -h, --help            show this help message and exit
  ```

This is a multi-modular software: 

![DRAMTOTAL](https://github.com/avera1988/NMBU-Bio-326/blob/main/images/DRAMtotal.jpg)

We need to run the annotation and then DRAM is capable to extract and classify Metabolic trends (e.g carbohydrate active enzymes (CAZy) coding genes) and produce visual plots with the distill command.

Let's display the help of these two commands:

```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login MAGs_gtdbk.dir]$ DRAM.py annotate --help
usage: DRAM.py annotate [-h] -i INPUT_FASTA [-o OUTPUT_DIR] [--min_contig_size MIN_CONTIG_SIZE] [--prodigal_mode {train,meta,single}]
                        [--trans_table {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}] [--bit_score_threshold BIT_SCORE_THRESHOLD]
                        [--rbh_bit_score_threshold RBH_BIT_SCORE_THRESHOLD] [--custom_db_name CUSTOM_DB_NAME] [--custom_fasta_loc CUSTOM_FASTA_LOC] [--gtdb_taxonomy GTDB_TAXONOMY]
                        [--checkm_quality CHECKM_QUALITY] [--use_uniref] [--low_mem_mode] [--skip_trnascan] [--keep_tmp_dir] [--threads THREADS] [--verbose]
 

(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login MAGs_gtdbk.dir]$ DRAM.py distill --help
usage: DRAM.py distill [-h] [-i INPUT_FILE] [-o OUTPUT_DIR] [--rrna_path RRNA_PATH] [--trna_path TRNA_PATH] [--groupby_column GROUPBY_COLUMN] [--custom_distillate CUSTOM_DISTILLATE]
                       [--distillate_gene_names] [--genomes_per_product GENOMES_PER_PRODUCT]
```

Basically the software only needs the fasta files, the checkM result and the Taxonomy. Then Distill needs the annotations from DRAM annotate to produce the summary tables and plots. The following is an script for running both the annotation and the distill part:

```bash
#!/bin/bash
#########################################################################
#	SLURM scrip for running DRAM annotator on Orion cluster
#		Dependencies: 
#	DRAM conda environment
#	gtdbtk.bac120.summary.tsv For taxonomy annotation
#	results.tsv	CheckM results for quality annotation
#
#	This script copies all the MAGs to a local disk on any Orion's node 
#	and runs locally on the node. At the end it copies all the results to the PEP
#	or any other path in the cluster into a DRAM.Results.dir folder.
#
#	to run: 
# sbatch  dram.GTDB.CM.SLURM.sh inputdir extension_of_fasta_files trans_table gtdbtk.tsv checkm.tsv
#eg: sbatch  dram.GTDB.CM.SLURM.sh MAGS fa 11 gtdbtk.bac120.summary.tsv ONT.tsv
#
# Author: Arturo Vera
# April 2021
#########################################################################

###############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=DRAM
#
## Wall time limit:
#SBATCH --time=24:00:00
#
## Other parameters:
#SBATCH --cpus-per-task 12
#SBATCH --mem=20G
#SBATCH --partition=smallmem

###########################################################

###Basic usage help for this script#######

print_usage() {
        echo "Usage: sbatch $0 inputdir extension_of_fasta_files trans_table gtdbtk.tsv checkm.tsv"
        echo "eg: sbatch $0 MAGS fa 11 gtdbtk.bac120.summary.tsv ONT.tsv"
}

if [ $# -le 1 ]
        then
                print_usage
                exit 1
        fi


## Set up job environment:

module --quiet purge  # Reset the modules to the system default
module load Miniconda3

##Activate conda environments

export PS1=\$
source activate /net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM

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

##Variables

input=$1 #Directory with genomes
ext=$2 #extention of the fasta files eg. .fa
trans_table=$3 #Translation table used by prodigal
gtdbtk=$4 #gtdbtk results table .tsv
checkm=$5 #checkm results table .tsv


## Copying data to local node for faster computation

cd $TMPDIR

#Check if $USER exists in $TMPDIR

if [[ -d $USER ]]
	then
        	echo "$USER exists on $TMPDIR"
	else
        	mkdir $USER
fi


echo "copying files to" $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID

cd $USER
mkdir tmpDir_of.$SLURM_JOB_ID
cd tmpDir_of.$SLURM_JOB_ID

#Copy the MAGs to the $TMPDIR

echo "copying MAGs to" $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID

cp -r $SLURM_SUBMIT_DIR/$input .
cp -r $SLURM_SUBMIT_DIR/$gtdbtk .
cp -r $SLURM_SUBMIT_DIR/$checkm .

echo "This are the files:"
ls -1

##################DRAM##############################

echo "DRAM started at"
date +%d\ %b\ %T

time DRAM.py annotate \
-i $input'/*.'$ext \
--trans_table $trans_table \
--gtdb_taxonomy $gtdbtk   \
--checkm_quality $checkm \
-o dram.annotation.$input.dir \
--threads $SLURM_CPUS_ON_NODE

echo "Distilling..."

time DRAM.py distill \
-i dram.annotation.$input.dir/annotations.tsv \
-o dram.genome_summaries.$input.dir \
--trna_path dram.annotation.$input.dir/trnas.tsv \
--rrna_path dram.annotation.$input.dir/rrnas.tsv

echo "DRAM finished at"
date +%d\ %b\ %T

mkdir DRAM.Results.$input.dir
mv dram.annotation.$input.dir DRAM.Results.$input.dir
mv dram.genome_summaries.$input.dir DRAM.Results.$input.dir

###########Moving results to PEP partition or anywhere the main script was submitted############

echo "moving results to" $SLURM_SUBMIT_DIR/$input

cd $TMPDIR/$USER/tmpDir_of.$SLURM_JOB_ID

time cp -r *.dir $SLURM_SUBMIT_DIR/$input

echo "DRAM results are in: " $SLURM_SUBMIT_DIR/$input/DRAM.Results.dir.$input.dir

####removing tmp dir. Remember to do this for not filling the HDD in the node!!!!###

cd $TMPDIR/$USER/
rm -r tmpDir_of.$SLURM_JOB_ID

echo "I've done at"
date

```
**You can copy this script to your folder by: ```cp /mnt/SCRATCH/bio326-21-0/MetagenomicMAGS/dram.GTDB.CM.SLURM.sh .```**

Let's run DRAM:

```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login MetagenomicMAGS]$ sbatch dram.GTDB.CM.SLURM.sh GoodQualityMAGs fa 11 gtdbtk.bac120.summary.tsv ONT_qa_bins.tsv
```


*DRAM took > 3 hrs to annotate all the genomes but you can obtain a copy of the results by ```cp -r /mnt/SCRATCH/bio326-21-0/MetagenomicMAGS/DRAM.Results.GoodQualityMAGs.dir .```*

Let's have a look of the DRAM results: 
```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login MetagenomicMAGS]$ cd DRAM.Results.GoodQualityMAGs.dir/
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login DRAM.Results.GoodQualityMAGs.dir]$ ls
dram.annotation.GoodQualityMAGs.dir  dram.genome_summaries.GoodQualityMAGs.dir
```
There are two directories: 
* dram.annotation.GoodQualityMAGs.dir: It has all the "raw" annotations, gene sequeces, protein preditions of the MAG's
* dram.genome_summaries.GoodQualityMAGs.dir: It has the destilled part of the genomes with the sorted metabolic functions


Take a look into the summaries directory:

```console
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login DRAM.Results.GoodQualityMAGs.dir]$ cd dram.genome_summaries.GoodQualityMAGs.dir/
(/net/cn-1/mnt/SCRATCH/bio326-21/GenomeAssembly/condaenvironments/DRAM) [bio326-21-0@login dram.genome_summaries.GoodQualityMAGs.dir]$ ls
genome_stats.tsv  metabolism_summary.xlsx  product.html  product.tsv
```
The files have different information:

* genome_stats.tsv: Basic annotaion stats of the genomes, as # of contigs/scaffolds, taxonomy, RNAgenes etc.
* metabolism_summary.xlsx: An excel file with all the Metabolic summary in each genome.
* product.html: Interactive heatmaps of the metabolic summaries
* product.tsv: Tables to reproduce the heatmaps of above

Let's check the summary of genomic stats:

```console
[bio326-21-0@login dram.genome_summaries.GoodQualityMAGs.dir]$ more genome_stats.tsv 
genome	number of scaffolds	taxonomy	completeness score	contamination score	5S rRNA	16S rRNA	23S rRNA	tRNA count	assembly quality
ONT_bin.1	1	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__CAG-41 sp900066215	97.99	2.85	4 present	4 present	4 present	57	med
ONT_bin.2	33	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides thetaiotaomicron	81.93	0.46	4 present	3 pres
ent	3 present	60	med
ONT_bin.3	1	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Enterocloster;s__Enterocloster sp000155435	98.67	1.27	5 present	5 pres
ent	5 present	69	med
ONT_bin.7	1	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides intestinalis_A	99.13	0.75	4 present	4 pres
ent	4 present	68	med
ONT_bin.8	25	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis	72.07	0.0	5 present	7 pres
ent	6 present	68	med
```
But if we would like to know how many genes were annotated per each genome, this info is not here. For doing that we need to extract them from the raw annotations. This file is located in ```dram.annotation.GoodQualityMAGs.dir``` so let's bring that file here and look into the file"

```console
[bio326-21-0@login dram.genome_summaries.GoodQualityMAGs.dir]$ cp ../dram.annotation.GoodQualityMAGs.dir/annotations.tsv .
[bio326-21-0@login dram.genome_summaries.GoodQualityMAGs.dir]$ head -2 annotations.tsv 
	fasta	scaffold	gene_position	start_position	end_position	strandedness	rank	kegg_id	kegg_hit	peptidase_id	peptidase_family	peptidase_hit	peptidase_RBH	peptidase_identity	peptidase_bitScore	peptidase_eVal	pfam_hits	cazy_hits	vogdb_description	vogdb_categories	heme_regulatory_motif_count	bin_taxonomy	bin_completeness	bin_contamination
ONT_bin.1_tig00000006_1	ONT_bin.1	tig00000006	1	3	218	-1	E													0d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Monoglobales_A;f__UBA1381;g__CAG-41;s__CAG-41 sp900066215	97.99	2.85
```

So this file, has all the genes that were predicted and annotated by DRAM. The first column has all the genes per each genome so we need to count those in order to know the number of genes per each genome. How to obtain that info, is not trivial and we need to do a little scripting for that. The following command line is an example of how we can use some simple ```bash``` commands we have already used (more, cut, awk, while, grep ...) to extract this information from the annotations file:

```console
[bio326-21-0@login dram.annotation.GoodQualityMAGs.dir]$ more annotations.tsv |cut -f 1|cut -d _ -f 1,2|sort|uniq|awk '{if($1 ~ /^[[:alnum:]]/) print $1}'|while read line; do echo $line;grep -c $line annotations.tsv; done
ONT_bin.1
2692
ONT_bin.2
4253
ONT_bin.3
5856
ONT_bin.7
5076
ONT_bin.8
2901
````

**Again, the main aim of this course is not to learn how to code, so do not panic, just look into the commands and try to make sense of them ...**

Although we can display the content of the *.tsv* files obtainded by DRAM here in the terminal, the  metabolism_summary.xlsx and product.html files are visually friendly, so it is recommendable to export these data to our personal computers and take a look. **A guide on how to copy files from Orion to our personal computers can be found in the [BacterialGenomeAssemblyMiniON](https://github.com/avera1988/NMBU-Bio-326/blob/main/Doc/BacterialGenomeAssemblyMiniON.md) document.**

Once in your computer, you can open the product.html, to explore the metabolic potential of your MAGs.

![dramhtml](https://github.com/avera1988/NMBU-Bio-326/blob/main/images/dramhtml.png)

is there any special metabolic pathway would you like to look at? 

Now is time for the funny part that is parsing the information and to interpret the biological meaning encoding in these MAGs ...

### Enjoy DRAM and have fun looking through your annotated MAGs


