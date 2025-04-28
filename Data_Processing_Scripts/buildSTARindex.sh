#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=03:00:00
# A script to make a star index

# Specify the reference genome directory
input="/scratch/alice/h/hp329/Steered_Project/hg19_genome"

# Specify the fasta files with the genome, and the annotation file
genomeFasta="$input/GRCh37.p13.genome.fa"
annotation="$input/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
index="starIndex"
mkdir -p "$index"

module load star

# STAR is run in generate genome index mode with the assembly from original analysis 
STAR --runMode genomeGenerate --genomeDir "$index" --genomeFastaFiles "$genomeFasta" --sjdbGTFfile "$annotation" --sjdbOverhang 100 --runThreadN 16
