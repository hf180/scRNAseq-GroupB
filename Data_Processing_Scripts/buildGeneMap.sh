#!/bin/bash
# A script to extract ENSG IDs and gene names from a GTF file into CSV format
input="/scratch/alice/h/hp329/Steered_Project/hg19_genome/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
# Specify path to the genome gtf annotation file

# Create CSV and name headers, then extract gene ID and gene name pairs
(echo "ENSG_ID,GENE_NAME"
	awk '$3 == "gene" { # Only pull data from the third column with genes
		match($0, /gene_id "([^"]+)\./, gid); # Expression to take gene ID up to the .
		match($0, /gene_name "([^"]+)"/, gname); # Expression to do the same with gene names
		if (gid[1] && gname[1]) print gid[1] "," gname[1]; # If both IDs are found print them into a CSV file
}' "$input" | sort -u) > gene_mapping.csv # Filter out duplicate lines and output to CSV


