#!/bin/bash

# Default values for optional parameters
threads=8
memory=8
build_ref=false

# Function to display usage message
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -x <index_file> -g <t2g_file> [-t threads] [-m memory] [-r] [-h]"
    echo "  -i: Input directory containing subfolders with fastq files (e.g., /mnt/f/scRNA_PreSAMBAM)"
    echo "  -o: Output directory for aligned results (e.g., /mnt/f/scRNA_Alligned)"
    echo "  -x: Path to the index file (e.g., index.idx)"
    echo "  -g: Path to the transcript-to-gene (t2g) mapping file (e.g., t2g.txt)"
    echo "  -t: Number of threads to use (default: 8)"
    echo "  -m: Memory in GB to allocate (default: 8)"
    echo "  -r: Build the reference files using 'kb ref' before running 'kb count'"
    echo "  -h: Display this help message"
    exit 1
}

# Parse command-line arguments
while getopts ":i:o:x:g:t:m:rh" opt; do
    case ${opt} in
        i ) # Input directory
            input_dir=$OPTARG
            ;;
        o ) # Output directory
            output_dir=$OPTARG
            ;;
        x ) # Index file
            index=$OPTARG
            ;;
        g ) # Transcript-to-gene file
            t2g=$OPTARG
            ;;
        t ) # Threads
            threads=$OPTARG
            ;;
        m ) # Memory
            memory=$OPTARG
            ;;
        r ) # Build reference files
            build_ref=true
            ;;
        h ) # Help
            usage
            ;;
        \? ) # Invalid option
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        : ) # Missing argument
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Ensure mandatory arguments are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    echo "Error: -i and -o options are required."
    usage
fi

# If the user opts to build reference files, perform this step
if $build_ref; then
    echo "Building reference files using kb ref..."
    
    # Fetch necessary files from gget and build reference files
    ref_files=$(gget ref --ftp -w dna,gtf homo_sapiens)
    fasta_file=$(echo "$ref_files" | grep -m1 ".fa")   # Extract the FASTA file path
    gtf_file=$(echo "$ref_files" | grep -m1 ".gtf")    # Extract the GTF file path

    echo "Fetched reference files: $fasta_file and $gtf_file"

    # Create the index and t2g files using kb ref
    kb ref -i "$output_dir/index.idx" -g "$output_dir/t2g.txt" -f1 "$fasta_file" "$gtf_file" \
        -t "$threads" -m "$memory"

    echo "Reference files built: index.idx and t2g.txt"

    # Update the index and t2g variables to point to the newly created files
    index="$output_dir/index.idx"
    t2g="$output_dir/t2g.txt"
fi

# Ensure the index and t2g files exist before proceeding
if [ ! -f "$index" ] || [ ! -f "$t2g" ]; then
    echo "Error: Index file (-x) and t2g file (-g) must be provided or built using -r."
    exit 1
fi

# Loop through each subfolder in the input directory
for subfolder in "$input_dir"/*; do
    if [ -d "$subfolder" ]; then  # Check if it's a directory
        # Construct the input file paths
        f1_file="$subfolder/*_f1.fastq.gz"
        r2_file="$subfolder/*_r2.fastq.gz"

        # Create the output directory for the current subfolder
        subfolder_name=$(basename "$subfolder")
        mkdir -p "$output_dir/$subfolder_name"

        # Run the kb count command
        echo "Running kb count for $subfolder_name with $threads threads and $memory GB memory..."
        kb count -i "$index" -g "$t2g" -x 10xv3 -t "$threads" -m "$memory" -o "$output_dir/$subfolder_name" $f1_file $r2_file --filter

        echo "Processed: $subfolder_name"
    fi
done

echo "All samples processed!"
