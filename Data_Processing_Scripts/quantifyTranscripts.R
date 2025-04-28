#!/usr/bin/env Rscript

library(tximport)
library(readr)

# Specify directories and 
input="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/salmon_QUANT_output"
output="/scratch/alice/h/hp329/Steered_Project/raw_mappedReads/salmon_quantified_GeneCounts"
txIndex="/scratch/alice/h/hp329/Steered_Project/Scripts/tx2gene.tsv"

# Create output directory if it does not exist
if (!dir.exists(output)) {dir.create(output, recursive =TRUE)}

txTable <- read_tsv(txIndex, col_names=c("TXNAME", "GENEID"))

# Find all Salmon quantification files under input directory
salmon_quants <- list.files(path=input, pattern="quant.sf", recursive=TRUE, full.names= TRUE)

# Extract sample names based on folder names
sampleVector <- basename(dirname(salmon_quants))
names(salmon_quants) <- sampleVector

# Use lists to record missing transcripts
missingTranscripts <- data.frame(Sample=character(), MissingCount=integer(), stringsAsFactors=FALSE)
allMissing <- list()

# Loop through each quantification file
for (quant in seq_along(salmon_quants)) {
	sample_quince <- names(salmon_quants)[quant]
	quant_file <- salmon_quants[quant]

	quant_data <- read_tsv(quant_file)
    # Identify transcripts present in quant.sf but missing from tx2gene table
	missing <- setdiff(quant_data$Name, txTable$TXNAME)
    # Record number of missing transcripts for current sample
	missingTranscripts <- rbind(missingTranscripts, data.frame(Sample=sample_quince, MissingCount=length(missing)))
    # Save list of missing transcript IDs for each sample
	allMissing[[sample_quince]] <- missing	

    # Run tximport to summarise transcript-level counts to gene-level counts
	txImp <- tximport(files=quant_file, type="salmon", tx2gene=txTable, countsFromAbundance="lengthScaledTPM")
	
	outputQuince <- file.path(output, paste0(sample_quince, "_geneCounts.csv"))
	write.csv(txImp$counts, file=outputQuince)
    # Write gene-level counts to CSV
}
# Write a summary table of number of missing transcripts per sample
summary_file <- file.path(output, "missingTranscripts_summary.csv")
write.csv(missingTranscripts, file=summary_file, row.names=FALSE)

# Convert all missing transcripts into a single list
missing_flat <- unlist(allMissing, use.names = FALSE)
# Calculate frequency of each missing transcript across all samples. 
missing_freq <- sort(table(missing_flat), decreasing = TRUE)
freq_file <- file.path(output, "missing_transcript_frequencies.csv")
write.csv(as.data.frame(missing_freq), file = freq_file, row.names = FALSE)
# Error handling step to identify dropped transcripts. This comes out to 0 in final run.
