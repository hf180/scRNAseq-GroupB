#!/bin/bash
#SBATCH --job-name=fastqc_run
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=22:00:00
# A script to run fastqc on the cleaned fastq files

# Setting input and output directories
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_prinseq_cleaned"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastQC_logs"

mkdir -p "$output"

module load fastqc

export output

# FastQC is run on each file
find "$input" -name "*.fastq" | parallel -j 6 fastqc {} -o "$output"
