#!/usr/bin/env python3
# Script to convert tximport gene count output into log10-TPM values

import pandas as pd
import numpy as np
import argparse
import re

    
def process_CPM(input, output):
    # Read input CSV file with gene counts
	df=pd.read_csv(input, index_col=0)
    # Clean row names
	df.index=df.index.str.replace("^gene-", "", regex=True)
	df.columns=[column.replace('"',"") for column in df.columns]

    # Convert all data to numeric
	df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
	total_reads = df.sum(axis=0)

    # Transform counts, and re-create new CSV files
	logCPM = np.log10(CPM+1)
	logCPM.insert(0, "Gene_ID", logCPM.index)
	logCPM.reset_index(drop=True, inplace=True)
	logCPM.to_csv(output,index=False)

# Set up argparse to take command input
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Convert tximport gene counts to log10 CPM")
	parser.add_argument("input", help="Input CSV gene count files")
	parser.add_argument("output", help="Output the converted csv file")
	args = parser.parse_args()
	process_CPM(args.input, args.output)
