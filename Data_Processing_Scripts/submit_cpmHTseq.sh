#!/bin/bash
#SBATCH --job-name=salmon_count
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to call a python script ro get CPM gene counts and add correct gene ids

# Load and activate conda environment
source /home/h/hp329/anaconda3/etc/profile.d/conda.sh

conda activate SRP_env

input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_GeneCounts_STAR"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/geneCounts_HTseqSTAR"
scripto="/scratch/alice/h/hp329/Steered_Project/Scripts/cpmHTSeq.py"
identifier="/scratch/alice/h/hp329/Steered_Project/Scripts/gene_mapping.csv"

# Call cpmHTSeq.py on every STAR output filex
mkdir -p "$output"

for file in "$input"/*.txt; do
        filename=$(basename "$file" .txt)
        output_file="$output/${filename}_CPM.csv"

        python3 "$scripto" "$file" "$output_file" --identifier_file "$identifier"
done
