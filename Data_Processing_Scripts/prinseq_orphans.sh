#!/bin/bash
#SBATCH --job-name=lonePairs
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59

# Set input and output directories, and load perl for running prinseq
module load perl

input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_diced"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_orphanage"
pathseq="/home/h/hp329/Software/prinseq-lite-0.20.4/prinseq-lite.pl"

# Create output directory if it doesn't exist yet
mkdir -p "$output"

# Export variables to use in parallel jobs
export input output pathseq

# Find all forward fastq files and process 6 jobs in parallel
find "$input" -name "*_1.fastq" | parallel -j 6 '
	forward={}
    # Extract unique file prefix for each cell
	file_prefix=$(basename "$forward" "_1.fastq")
	# Find the matching paired end read
    reverse="$input/${file_prefix}_2.fastq"
    # Run prinseq to filter reads
	perl "${pathseq}" -fastq "${forward}" -fastq2 "${reverse}" \
	-min_len 30 \
	-out_good "${output}/${file_prefix}_adopted" -out_bad null -out_format 3
'
