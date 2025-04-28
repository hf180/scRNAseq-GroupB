#!/bin/bash
#SBATCH --job-name=fastq_pr
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to run prinseq to filter fastq files

# Load perl for prinseq
module load perl

# Set input and output directories, and path to prinseq install location
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_out"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_prinseq_cleaned"
pathseq="/home/h/hp329/Software/prinseq-lite-0.20.4/prinseq-lite.pl"
mkdir -p "$output"

# Find all forward and reverse files, process in 6 jobs at once
find "$input" -name "*_1.fastq" | parallel -j 6 '
	forward={}
	file_prefix=$(basename "$forward" "_1.fastq")

    # Define the matching reverse paired read
	reverse="'$input'/${file_prefix}_2.fastq"

    # Run prinseq to filter reads
	perl "'${pathseq}'" -fastq "${forward}" -fastq2 "${reverse}" \
	-min_len 30 -trim_left 10 -trim_qual_right 25 -lc_method entropy -lc_threshold 65 \
	-out_good "'${output}'/${file_prefix}_cleaned" -out_bad null -out_format 3
'
