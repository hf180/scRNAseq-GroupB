#!/bin/bash
#SBATCH --job-name=salmon_count
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hp329@student.le.ac.uk
#SBATCH --export=NONE
#SBATCH --time=23:59:59
# A script to call an R script to quantify transcript abundance and gene counts in salmon quant files

source /home/h/hp329/anaconda3/etc/profile.d/conda.sh

conda activate SRP_env

Rscript /scratch/alice/h/hp329/Steered_Project/Scripts/quantifyTranscripts.R
