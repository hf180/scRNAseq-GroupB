#!/bin/bash
#SBATCH --job-name=STAR_alignment
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to run STAR on fastq data

# Specify input and output directories as well as genome index file path
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_PREPPED"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/fastq_STAR"
starryindex="/scratch/alice/h/hp329/Steered_Project/Scripts/starIndex"

mkdir -p "$output"

module load star

# Loop over all forward reads
for forward in "$input"/*_1.fastq; do
        # Extract corresponding reverse read
        file_prefix=$(basename "${forward}" "_1.fastq")
        reverse="${input}/${file_prefix}_2.fastq"

        # Run STAR with 16 threads
        STAR --runThreadN 16 \  
        --genomeDir "$starryindex" \
        --readFilesIn "$forward" "$reverse" \
        --outFileNamePrefix "$output/${file_prefix}_" \
        --outSAMtype BAM SortedByCoordinate \   # Output as BAM files
        --outFilterType BySJout \   # Filter alignments by splicing junctions
        --outFilterMultimapNmax 20 \    # Set a maximum of 20 multiple alignments for each read
        --alignSJoverhangMin 8 \    # Set minimum overhang across novel splicing junctions 
        --alignSJDBoverhangMin 1 \  # Set minimum overhang across annotated splicing junctions
        --outFilterMismatchNmax 999 \   # Set maximum number of mismatches
        --outFilterMismatchNoverLmax 0.04 \ # Set maximum number mismatches as a proportion of read length
        --alignIntronMin 20 \   # Set minimum size of introns
        --alignIntronMax 1000000 \  # Maximum size of introns
        --alignMatesGapMax 1000000 \    # Set maximum allowed gap between paired reads
        --outSAMstrandField intronMotif \   # Include stranded info in output BAM
        --outSAMattributes Standard # Include standard set of attributes in output BAM
done
