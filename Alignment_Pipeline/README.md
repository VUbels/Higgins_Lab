This script automates the alignment of scRNA-seq fastq files using kb count. It processes multiple subfolders in an input directory and saves the results to a specified output directory. Optionally, the script can build reference files using kb ref. Run Alignment_DropletRemoval if filter is not used and/or droplet based sequencing technique is used. 

Usage

```
bash
./alignment_script.sh -i <input_dir> -o <output_dir> -x <index_file> -g <t2g_file> [-t threads] [-m memory] [-r] [-h]
```

Required Parameters:

    -i: Input directory containing subfolders with fastq files.
    -o: Output directory for aligned results.
    -x: Path to index file (optional if using -r).
    -g: Path to transcript-to-gene mapping file (optional if using -r).

Optional Parameters:

    -t: Number of threads (default: 8).
    -m: Memory in GB (default: 8).
    -r: Build reference files using kb ref.
    -h: Display help message.

List of supported single-cell technologies

```
short name       description
----------       -----------
10xv1            10x version 1 chemistry
10xv2            10x version 2 chemistry
10xv3            10x version 3 chemistry
Bulk             Bulk RNA-seq
SmartSeq2        Smart-seq2
BDWTA            BD Rhapsody WTA
CELSeq           CEL-Seq
CELSeq2          CEL-Seq version 2
DropSeq          DropSeq
inDropsv1        inDrops version 1 chemistry
inDropsv2        inDrops version 2 chemistry
inDropsv3        inDrops version 3 chemistry
SCRBSeq          SCRB-Seq
SmartSeq3        Smart-seq3
SPLiT-seq        SPLiT-seq
SureCell         SureCell for ddSEQ
Visium           10x Visium Spatial Transcriptomics
```

```
Output example

/mnt/f/scRNA_Aligned/
│
├── TEST/
│   ├── counts_unfiltered/
│   │   ├── cells_x_genes.mtx
│   │   ├── genes.names.txt
│   │   ├── cells_x_genes.genes.names.txt
│   │   ├── barcodes.txt
│   ├── counts_filtered/
│   │   ├── cells_x_genes.mtx
│   │   ├── genes.names.txt
│   │   ├── cells_x_genes.genes.names.txt
│   │   ├── barcodes.txt
│   ├── counts_filtered_droplet_removed/
│   │   ├── filtered_cells_x_genes.mtx
│   │   ├── filtered_genes.names.txt
│   │   ├── filtered_cells_x_genes.genes.names.txt
│   │   ├── filtered_barcodes.txt
│   │   ├── pre_post_filter_dimensions.txt
│   │   ├── rank_total_UMI_plot.png


```
