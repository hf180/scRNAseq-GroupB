#!/bin/bash
#SBATCH --job-name=salmonIndex
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=02:00:00
# A script to make a salmon index

# Specifying file paths to genome fasta and RNA fasta files
fishIndex="/scratch/alice/h/hp329/Steered_Project/hg38/ncbi_dataset/data/GCF_000001405.40/"

genome="/scratch/alice/h/hp329/Steered_Project/hg38/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
transcripts="/scratch/alice/h/hp329/Steered_Project/hg38/ncbi_dataset/data/GCF_000001405.40/rna.fna"

module load salmon
# Load salmon, create output directory and navigate into it
mkdir -p salmonIndex
cd salmonIndex

# Generate a list of decoy sequences from the genome fasta
grep "^>" "$genome" | cut -d " " -f 1 | sed 's/>//' > decoys.txt 

# Concatenate genome and transcripts into a single fasta file for indexing
cat "$transcripts" "$genome" > salmon_decoy.fa 

# Build the salmon index in decoy aware mode
salmon index -t salmon_decoy.fa -d decoys.txt \
	-i salmonIndex_refseq \
	--threads 8 --keepDuplicates
