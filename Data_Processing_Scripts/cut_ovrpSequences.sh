#!/bin/bash
#SBATCH --job-name=rmv_sequences
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to run Cutadept to remove overrepresented sequences from paired end read files

# Input and output directories
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_prinseq_cleaned"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_OVRP_diced"

# Text file containing overrepresented sequences, extracted from MultiQC report
OVRPseq="/scratch/alice/h/hp329/Steered_Project/Scripts/OVRP.txt"

mkdir -p "$output"

module load cutadapt

# Export variables needed for parallel jobs
export OVRPseq output

# Define the function to be run in parallel
remove_OVRP() {
	forward="$1"
	reverse="$2"

    # Identify the unique file prefix and assign the forward and reverse reads to variables
	file_prefix=$(basename "$forward" | sed 's/_cleaned_1.fastq//')
	outforward="$output/${file_prefix}_diced_1.fastq"
	outreverse="$output/${file_prefix}_diced_2.fastq"

    # Run Cutadapt to remove overrepresented sequences in text file
	cutadapt -e 0.15 -m 30 \
	-b file:"$OVRPseq" -B file:"$OVRPseq" \
	-o "$outforward" -p "$outreverse" \
	"$forward" "$reverse"
}

# Export the function 
export -f remove_OVRP

# Run remove_OVRP in 6 parallel jobs
parallel -j 6 --link remove_OVRP ::: "$input"/*_cleaned_1.fastq ::: "$input"/*_cleaned_2.fastq
