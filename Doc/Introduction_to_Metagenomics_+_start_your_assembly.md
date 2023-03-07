


## Filter raw reads

```

filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    output/samtools/mapped.sorted.bam.fastq.gz \
    | gzip \
    > output/filtlong/output.fastq.gz

```
