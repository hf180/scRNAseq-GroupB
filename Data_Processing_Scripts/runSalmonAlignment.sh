#!/bin/bash
#SBATCH --job-name=salmon_run
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to run salmon quantification on the fastq reads

# Specify the input and output file paths, as well as the salmon index
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_PREPPED"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/salmon_QUANT_output"
salmonIndex="/scratch/alice/h/hp329/Steered_Project/Scripts/salmonIndex/salmonIndex_refseq"

mkdir -p "$output"

module load salmon

# Loop over all forward files and extract the corresponding reverse
for forward in "$input"/*_1.fastq; do
	file_prefix=$(basename "$forward" "_1.fastq")
        reverse="${input}/${file_prefix}_2.fastq"
        outDir="${output}/${file_prefix}"

    # Run Salmon in quantification mode
	salmon quant -i "$salmonIndex" \
		-l A -1 "$forward" -2 "$reverse" \ # Set forward and reverse reads, and automatically detect library type
		-p 16 --validateMappings -o "$outDir" # Use selective alignment mapping and set thread usage to 16
done

