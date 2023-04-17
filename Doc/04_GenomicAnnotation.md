
### So far have obtained indiviual metagenome assembled genomes (MAGs) from the rummen metagenome. Let's review our workflow to see what is the next step: 

![workflow](https://github.com/avera1988/NMBU-Bio-326/blob/main/images/wrokflowmetagenome.png) 

After binning with [MetaBat2](https://bitbucket.org/berkeleylab/metabat/src/master/) we obtained 9 bins. We used [CheckM](https://ecogenomics.github.io/CheckM/) to asses the quality of these MAGs.

Metabat2 produce a set of bins as the following example:

```bash
ls -1|grep .fa|sort -V
bin.1.fa
bin.2.fa
bin.3.fa
bin.4.fa
bin.5.fa
bin.6.fa
bin.7.fa
bin.8.fa
bin.9.fa
bin.10.fa
bin.11.fa
bin.12.fa
bin.13.fa
bin.14.fa
bin.15.fa
bin.16.fa
bin.17.fa
bin.18.fa
bin.19.fa
bin.20.fa
bin.21.fa
bin.22.fa
bin.23.fa
bin.24.fa
bin.25.fa
bin.26.fa
```
We then can use the report from the ```assemblycomparator2``` to check the quality report by CheckM, lets use 

