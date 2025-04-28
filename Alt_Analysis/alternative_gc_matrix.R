#Obtaining the gene count matrix using mapped read data generated using an alternative approach

library(tidyr)
library(dplyr)
library(readr)
library(biomaRt)
library(stringr)

#changing working directory
setwd('/home/hf180/Desktop/scRNAseq')

#obtaining files from github
system("git clone https://github.com/hf180/scRNAseq-GroupB.git")

#set paths
mapping_file <- '/home/hf180/Desktop/scRNAseq/SraRunTable.csv'
input_dir <- '/home/hf180/Desktop/scRNAseq/scRNAseq-GroupB/finished_GeneMatrices/alternative_GeneCounts'
output_dir <- "/home/hf180/Desktop/scRNAseq/scRNAseq-GroupB/finished_GeneMatrices/alternative_GeneCounts/cleaned_counts"

#create output dir
if (!dir.exists(output_dir)) dir.create(output_dir)

#load mapping 
map <- read.csv(mapping_file, stringsAsFactors = FALSE)

#rename columns for easier reference
map <- map %>%
  rename(
    SRR_ID = Run,
    GSE_ID = GEO_Accession..exp.  # Adjust this if your GSEs are in another column
  )

#list all input files
all_files <- list.files(input_dir, pattern = "_alternative_geneCounts.csv$", full.names = TRUE)

#process each file
for (file in all_files) {
  #extract SRR ID from filename
  filename <- basename(file)
  srr <- sub("_alternative_geneCounts.csv$", "", filename)
  
  #look up GSE ID
  gse <- map$GSE_ID[match(srr, map$SRR_ID)]
  #output filename
  output_file <- file.path(output_dir, paste0(gse, "_", srr, ".csv"))
  
  #read, clean, and write
  counts <- read_csv(file, col_names = FALSE)
  colnames(counts) <- c("Gene", "Count")
  
  counts_clean <- counts %>%
    group_by(Gene) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  write_csv(counts_clean, output_file)
}

#get protein-coding gene list from Ensembl
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")

# Check if connection was successful
if (!is.null(mart)) {
  # Retrieve protein-coding genes
  coding_genes <- getBM(
    attributes = c("hgnc_symbol"),
    filters = "biotype",
    values = "protein_coding",
    mart = mart
  ) %>% pull(hgnc_symbol) %>% unique()
  
  # Print the first few coding genes
  print(head(coding_genes))
} else {
  message("Failed to connect to Ensembl.")
}

#find all cleaned files
files <- list.files("/home/hf180/Desktop/scRNAseq/scRNAseq-GroupB/finished_GeneMatrices/alternative_GeneCounts/cleaned_counts", pattern = "*.csv", full.names = TRUE)

#named the files
names(files) <- stringr::str_split(files, pattern = "/", simplify = TRUE)[ ,7] %>% stringr::str_replace(".csv","")

#create an empty dataframe
results<- data.frame()

#loop through each file, read it and store it in the list
for (file in files){
  x <- read.csv(file, sep = ",", header = TRUE, stringsAsFactors = FALSE) %>% 
    filter(!grepl('rRNA', Gene, ignore.case = TRUE)) %>% 
    filter(Gene %in% coding_genes)
  #if the file has two columns, rename them
  if (ncol(x) == 2) {
    colnames(x) <- c("Gene", "Expression")
  }
  #add the sample name as a new column
  sample_name <- stringr::str_replace(file, "/home/hf180/Desktop/scRNAseq/scRNAseq-GroupB/finished_GeneMatrices/alternative_GeneCounts/cleaned_counts", "") %>%
    stringr::str_replace(".csv", "")
  x$sample <- sample_name
  #combine the dataframe
  results <- bind_rows(results, x)
}

#reshape to have genes as rows and samples as columns
#handle duplicates
results_wide <- results %>%
  group_by(Gene, sample) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = sample, values_from = Expression)

#save the cleaned dataframe to CSV
write.csv(results_wide, "/home/hf180/Desktop/scRNAseq/scRNAseq-GroupB/finished_GeneMatrices/alternative_GeneCounts/alt_results.csv", row.names = FALSE)

