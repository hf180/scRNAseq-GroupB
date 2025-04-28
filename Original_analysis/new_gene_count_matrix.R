#this is code used for the scRNA-seq project. 
#this code converts all gene matrices to a big gene count matrix, where rows are genes and cells are columns.
#as this matrices come with their SRR code they were converted to their allocated GSM code, duplicated were dropped and genes were filtered to protein coding ones. 
#this code can be used for both original and alternitive pipeline

#load these libraries
library(tidyr)
library(dplyr)
library(readr)
library(biomaRt)
library(stringr)

#set paths
mapping_file <- '/home/ex9/SraRunTable.csv'
input_dir <- "/home/ex9/Steered_project/original_geneCounts"
output_dir <- "/home/ex9/Steered_project/original_geneCounts/cleaned_counts"

#create output dir
if (!dir.exists(output_dir)) dir.create(output_dir)

#load mapping 
map <- read.csv(mapping_file, stringsAsFactors = FALSE)

#rename columns for easier reference
map <- map %>%
  rename(
    SRR_ID = Run,
    GSE_ID = GEO_Accession..exp.)

#list all input files
all_files <- list.files(input_dir, pattern = "_original_geneCounts\\.csv$", full.names = TRUE)

#process each file
for (file in all_files) {
  #extract SRR ID from filename
  filename <- basename(file)
  srr <- sub("_original_geneCounts\\.csv$", "", filename)
  
  #look up GSE ID
  gse <- map$GSE_ID[match(srr, map$SRR_ID)]
  #output filename
  output_file <- file.path(output_dir, paste0(gse, "_", srr, ".csv"))
  
  #read, clean, and write
  counts <- read_csv(file, col_names = FALSE, show_col_types = T)
  colnames(counts) <- c("Gene", "Count")
  
  counts_clean <- counts %>%
    group_by(Gene) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  write_csv(counts_clean, output_file)
}

#get protein-coding gene list from Ensembl
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
coding_genes <- getBM(
  attributes = c("hgnc_symbol"),
  filters = "biotype",
  values = "protein_coding",
  mart = mart
) %>% pull(hgnc_symbol) %>% unique()

#find all cleaned files
files <- list.files("/home/ex9/Steered_project/original_geneCounts/cleaned_counts", pattern = "*.csv", full.names = TRUE)

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
  sample_name <- stringr::str_replace(file, "/home/ex9/Steered_project/original_geneCounts/cleaned_counts", "") %>%
    stringr::str_replace(".csv", "")
  x$sample <- sample_name
  #combine the dataframe
  results <- bind_rows(results, x)
}

#if you want to reshape to have genes as rows and samples as columns
results_wide <- results %>%
  pivot_wider(names_from = sample, values_from = Expression)

#save the dataframe to a CSV file
write.csv(results_wide, "/home/ex9/Steered_project/original_geneCounts/cleaned_counts/original_results.csv", row.names = FALSE)
