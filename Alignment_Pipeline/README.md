Aligns FASTQ sequence data using a simple bash script.

Download reference genome : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

Assign reference genome through -r
Assign input/output folders through -i -o
Assign threads through -t

Simply call in terminal using:

```
python Alignment_SAMBAM.py \
  -r ~/reference/hg38/hg38.fa \
  -i /mnt/f/scRNA_PreSAMBAM \
  -o /mnt/f/scRNA_PreSAMBAM/aligned_output \
  -t 4
```

