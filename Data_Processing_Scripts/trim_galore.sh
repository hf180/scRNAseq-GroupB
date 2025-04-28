#!/bin/bash
#SBATCH --job-name=rmv_adaptors
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to run trim galore and remove nextera adaptors

# Specify the input and output directories
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_orphanage"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_PREPPED"

# Assign file path for Trim Galore installation
trimseq="/home/h/hp329/Software/TrimGalore-0.6.10/trim_galore"

# Create output file
mkdir -p "$output"

export input output trimseq

# Define function to run Trim Galore in parallel
trim_galley() {
    # Extract unique file prefixes and pair with the corresponding file
	forward="$1"
	file_prefix=$(basename "${forward}" "_1.fastq")
	reverse="${input}/${file_prefix}_2.fastq"

    # Run Trim Galore with the built in option to remove Nextera sequencing adaptors
	"${trimseq}" --paired --nextera --stringency 1 --output_dir "$output" "$forward" "$reverse"
}
# Export the functio to use in parallel
export -f trim_galley

# Call the function on each forward read fastq file into 6 parallel jobs
find "$input" -name "*_1.fastq" | parallel -j 6 trim_galley
