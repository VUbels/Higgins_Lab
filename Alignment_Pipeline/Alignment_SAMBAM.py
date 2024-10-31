#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
from pathlib import Path
import logging
import glob

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def check_dependencies():
    """Check if required tools are installed."""
    for tool in ['bwa', 'samtools']:
        try:
            subprocess.run([tool, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error(f"{tool} not found. Please install it using conda: conda install -c bioconda {tool}")
            sys.exit(1)

def index_reference(reference_path):
    """Index reference genome if not already indexed."""
    if not os.path.exists(f"{reference_path}.bwt"):
        logging.info("Indexing reference genome...")
        try:
            subprocess.run(['bwa', 'index', reference_path], check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error indexing reference genome: {e}")
            sys.exit(1)

def find_fastq_pairs(base_dir):
    """Find all FASTQ pairs in nested directories."""
    pairs = []
    # Walk through all subdirectories
    for subdir in os.listdir(base_dir):
        full_path = os.path.join(base_dir, subdir)
        if os.path.isdir(full_path):
            # Look for f1 and r2 files in this directory
            f1_files = glob.glob(os.path.join(full_path, "*_f1.fastq*"))
            for f1 in f1_files:
                # Construct expected r2 filename
                r2 = f1.replace("_f1.fastq", "_r2.fastq")
                if os.path.exists(r2):
                    pairs.append((f1, r2, subdir))
    return pairs

def process_fastq_pair(f1_path, r2_path, sample_name, reference_path, output_dir, threads=4):
    """Process a pair of FASTQ files to generate sorted BAM."""
    try:
        # Create sample-specific output directory
        sample_output_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Define output files
        sam_output = os.path.join(sample_output_dir, f"{sample_name}.sam")
        bam_output = os.path.join(sample_output_dir, f"{sample_name}.sorted.bam")
        
        logging.info(f"Processing sample: {sample_name}")
        
        # Generate SAM file
        logging.info(f"Generating SAM file for {sample_name}")
        bwa_cmd = ['bwa', 'mem', '-t', str(threads), reference_path, f1_path, r2_path]
        with open(sam_output, 'w') as sam_file:
            subprocess.run(bwa_cmd, stdout=sam_file, check=True)
        
        # Convert SAM to sorted BAM
        logging.info(f"Converting SAM to sorted BAM for {sample_name}")
        subprocess.run([
            'samtools', 'sort',
            '-@', str(threads),
            '-o', bam_output,
            sam_output
        ], check=True)
        
        # Index BAM
        logging.info(f"Indexing BAM file for {sample_name}")
        subprocess.run(['samtools', 'index', bam_output], check=True)
        
        # Generate stats
        stats_file = os.path.join(sample_output_dir, f"{sample_name}.flagstat.txt")
        logging.info(f"Generating statistics for {sample_name}")
        with open(stats_file, 'w') as f:
            subprocess.run(['samtools', 'flagstat', bam_output], stdout=f, check=True)
        
        # Optionally remove SAM file to save space
        # os.remove(sam_output)  # Uncomment if you want to remove SAM files
        
        logging.info(f"Completed processing {sample_name}")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Error processing {sample_name}: {e}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error processing {sample_name}: {e}")
        return False
    return True

def main():
    parser = argparse.ArgumentParser(description='Process FASTQ files to SAM/BAM files')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('-i', '--input_dir', required=True, help='Base directory containing FASTQ files in subdirectories')
    parser.add_argument('-o', '--output_dir', default='processed_alignments', help='Output directory for SAM/BAM files')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to use')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging()
    
    # Check dependencies
    check_dependencies()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Index reference genome if needed
    index_reference(args.reference)
    
    # Find all FASTQ pairs
    fastq_pairs = find_fastq_pairs(args.input_dir)
    
    if not fastq_pairs:
        logging.error(f"No matching FASTQ pairs found in {args.input_dir}")
        sys.exit(1)
    
    # Process each pair
    for f1, r2, sample_name in fastq_pairs:
        logging.info(f"Found pair: {f1} and {r2} for sample {sample_name}")
        success = process_fastq_pair(f1, r2, sample_name, args.reference, args.output_dir, args.threads)
        if not success:
            logging.warning(f"Failed to process {sample_name}")
    
    logging.info("All samples processed!")

if __name__ == '__main__':
    main()
