#!/usr/bin/env python3

# Script to convert HTSeq-count output into CPM and log10-CPM, as in original analysis

import pandas as pd
import numpy as np
import argparse

def process_CPM(input, output, identifier_file=None):
    # Read the HTSeq-count file into a dataframe	
    df = pd.read_csv(input, sep='\t', header=None, names=["gene_id", "count"])
    # Remove rows where gene_id starts with underscores
	df = df[~df["gene_id"].str.startswith('__')]

    # Ensure count column is numeric data
	df["count"]=pd.to_numeric(df["count"], errors="coerce").fillna(0)
	total_reads = df["count"].sum()
    # If no reads, set CPM and log10_CPM to zero
	if total_reads == 0:
		df["CPM"]=0
		df["log10_CPM"]=0
	else:
        # Else perform normalisation as in original analysis
		df["CPM"]=df["count"]/total_reads * 1e6
		df["log10_CPM"]=np.log10(df["CPM"]+1)

    # Rename columns by mapping ENSG IDs to gene names
	if identifier_file:
		identifier = pd.read_csv(identifier_file, dtype=str)
        # Remove version numbers from gene IDs
		df["gene_id"] = df["gene_id"].astype(str).str.replace(r"\.\d+$","",regex=True)
		if "ENSG_ID" in identifier.columns and "GENE_NAME" in identifier.columns:
            # If both columns are valid merge the dataframes for proper gene name and count matrix
			df = df.merge(identifier, left_on="gene_id", right_on="ENSG_ID", how="left")
            # Assemble gene count matrix with recalculated values
			df = df[["gene_id","GENE_NAME", "count","CPM","log10_CPM"]]
		else:
			print("Incorrect column labels")
	df.to_csv(output, index=False)

# Set up command line argument parsing with argparse
if __name__=="__main__":
	parser= argparse.ArgumentParser(description="Convert htseq counts and change columns to proper gene names")
	parser.add_argument("input", help="Input the txt htseqcount file")
	parser.add_argument("output", help="Output as a csv file")
	parser.add_argument("--identifier_file", help="Get correct gene id into output csv files",default=None)

	args=parser.parse_args()
	process_CPM(args.input, args.output, args.identifier_file)
