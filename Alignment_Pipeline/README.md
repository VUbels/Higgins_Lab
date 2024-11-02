This script automates the alignment of scRNA-seq fastq files using kb count. It processes multiple subfolders in an input directory and saves the results to a specified output directory. Optionally, the script can build reference files using kb ref. Run Alignment_DropletRemoval if filter is not used and/or droplet based sequencing techniquen is used. 

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

```



```
