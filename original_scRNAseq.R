#This code attempts to replicate the scRNA-seq analysis of Darmanis et al., 
#This code attempts to stick as much as possible to the what the authors did and try to do it with the tools that were available in 2015. 
#SCDE, flexmix were the only packaging needing downgrading the rest are up to date
#This project started on the 12th of March 2025.

#**************************************

#load all relevant packages
library(readr)
library(lattice)
library(flexmix)
library(scde)
library(parallel)
library(Rtsne)
library(mclust)
library(igraph)
library(plotly)
library(FactoMineR)
library(ggplot2)
library(SingleCellExperiment)
library(GEOquery)
library(pheatmap)
library(Matrix.utils)
library(magrittr)
library(purrr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(tidyr)
library(scales)

#read the dataframe 
df <- read.csv("/home/ex9/Steered_project/original_geneCounts/cleaned_counts/original_results.csv", sep =",", header = TRUE, row.names = 1)

#reading counts matrix into R
mat <- as.matrix(df)
mat[is.na(mat)] <- 0 #checking and removing any null values
#turning values to integers
mat<-apply(mat,2,function(x) {storage.mode(x) <- 'integer'; x})
#check the number of rows and columns
dim(mat)

#load data as a single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = mat))

#to get the count matrix again
count_matrix <- counts(sce)

#simple transformation of the data
assay(sce, "logcounts") <- log2(counts(sce) + 1)

#further more matrix stats
colMeans(counts(sce))
#adding the means as a new column 
colData(sce)$mean_counts <- colMeans(counts(sce))
#adding another column with total counts
colData(sce)$total_data <- colSums(counts(sce))
#adding another column but this one is cell per million
assay(sce, "cpm") <- sweep(counts(sce),2,sce$total_data/1e6,'/')
#check that columns are 1e6 
colSums(cpm(sce))

#calculate the mean counts per gene
gene_means <- rowMeans(counts(sce))
#turning them to TRUE/FALSE
gene_means > 0.01
#counts of at least 1 
counts(sce) > 0
#total number of detected genes per cell
total_detected_per_cell <- colSums(counts(sce) > 0)
#print the first 50 values
total_detected_per_cell[1:50]

#load metadata attempt
gse <- getGEO("GSE67835", GSEMatrix = TRUE)
#check structure
length(gse)#would be 2 as there are 2 platforms
names(gse) #print IDs of the files
#access each file's data
gse15520 <- gse[[1]]
gse18573 <- gse[[2]]
#extracting the metadata of each platform
metadata_15520 <- pData(phenoData(gse15520))
metadata_18573 <- pData(phenoData(gse18573))
#checking for any overlapping samples if any 
shared_samples <- intersect(rownames(metadata_15520), rownames(metadata_18573))
#no overlaps so the metadatas will be combined together
metadata <- rbind(metadata_15520, metadata_18573)
#extract critical variables
metadata$cell_type <- gsub(".*: ", "", metadata$characteristics_ch1.1)
metadata$age <- gsub(".*: ", "", metadata$characteristics_ch1.2)
head(metadata[, c("title", "cell_type", "age", "geo_accession")])

#add metadata to the SingleCellExperiment
sce$sample_id <- gsub("_.*", "", colnames(sce))
sce$cell_type <- metadata[match(sce$sample_id, metadata$geo_accession), "cell_type"]
#add metadata to colData
colData(sce) <- cbind(colData(sce), metadata[match(sce$sample_id, metadata$geo_accession), c("cell_type", "age")])

#making some graphs 
cell_info <- as.data.frame(colData(sce))
head(cell_info)
#histograms
ggplot(data = cell_info, aes(x = total_data/1e6)) +
  geom_histogram(fill = "black", bins = 8) +
  theme(axis.title.x = element_text(angle = 0))

#clean the original matrix
cd <- clean.counts(mat)
dim(cd)

#building an error model 
#err <- knn.error.models(ss_counts, groups = NULL, k = ncol(cd)/2, n.cores=1, verbose = 1)
error_model <- scde.error.models(counts = cd,groups = NULL, n.cores = 1, threshold.segmentation = TRUE,
                                 save.crossfit.plots = FALSE, save.model.plots = FALSE,
                                 verbose = 1)
save(error_model, file = "fitted_errmod_original.Rdata")

#filter out cells that dont show positive correlation
valid_cells <- error_model$corr.a > 0
table(valid_cells)
error_model <- error_model[valid_cells,]

#loading the error model in the sce class
common_cells <- colnames(sce)[colnames(sce) %in% rownames(error_model)]
error_model <- error_model[common_cells,]
sce <- sce[,common_cells]
colData(sce) <- cbind(colData(sce), error_model)

#calculating pairwise distances
e_prior <- scde.expression.prior(models = error_model, counts = cd)
jp <- scde.posteriors(models = error_model, counts = cd, e_prior,return.individual.posterior.modes = TRUE,
                      n.cores =1)
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
pred_fail <- scde.failure.probability(models = error_model, magnitudes = jp$jp.modes)
pred_sfail <- scde.failure.probability(models = error_model, counts = cd)

#weight matrix
matw <- 1-sqrt(pred_sfail*sqrt(pred_sfail*pred_fail))
#magnitude matrix(using individual posterior modes here)
mat <- log10(exp(jp$modes)+1)
#weighted distance
cell.names <- colnames(cd)
names(cell.names)<-cell.names

require(boot)
mode_fail_dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1){
  unlist(lapply(cell.names,function(nam2){
    corr(cbind(mat[,nam1],mat[,nam2]),w = sqrt(sqrt(matw[,nam1]*matw[,nam2])))
  }))
},mc.cores = 1)),upper = FALSE)

#save the distances as a file
save(mode_fail_dist, file = "mode_fail_dist_original.RData")

#perform hierarchial clustering 
clust<- hclust(mode_fail_dist,method = "ward.D2")
plot(as.dendrogram(clust),cex=0.7)

#running Rtsne package and performing reduction of dimensions
tsne_map<- Rtsne(mode_fail_dist,dims = 3, perplexity=30)
#extracting clusters centroids
tsne_embedding <- tsne_map$Y

#using mcLust and clustering the cells
seqBIC <- mclustBIC(tsne_embedding, G= 1:40)
plot(seqBIC)

#fit a GMM to the t-SNE embedding
gmm_result <- Mclust(tsne_embedding, G = 1:40)

#view the optimal number of clusters selected by BIC
print(gmm_result$G) 

#view the cluster assignments
cluster_assignments <- gmm_result$classification
print(cluster_assignments)

#plotting box plots of clustering uncertainty 
#combined cluster assignments and uncertainty into a dataframe
uncertainty_df <- data.frame( Cluster = as.factor(gmm_result$classification), 
                              Uncertainty = gmm_result$uncertainty)

#visualisisng the boxplots
ggplot(uncertainty_df, aes(x = Cluster, y = Uncertainty, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Clustering Uncertainty",
    x = "Cluster",
    y = "Uncertainty"
  ) +
  scale_fill_brewer(palette = "Set2")

#plot the t-SNE embedding with cluster assignments
tsne_df <- data.frame(
  X = tsne_map$Y[, 1],
  Y = tsne_map$Y[, 2],
  Z = tsne_map$Y[, 3]
)
head(tsne_df)

colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE")
hover_text <- paste(
  "Digit:",
  "Dimension 1:", round(tsne_df$X, 3),
  "Dimension 2:", round(tsne_df$Y, 3),
  "Dimension 3:", round(tsne_df$Z, 3)
)

plot_ly(
  data = tsne_df,
  x = ~X,
  y = ~Y,
  z = ~Z,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 6),
  text = hover_text,
  hoverinfo = "text",
  color = ~cluster_assignments,
  colors = colors
) %>% layout(
  title = "t-SNE 3-Dimensional Digit Visualization",
  scene = list(
    xaxis = list(title = "t-SNE Dimension 1"),
    yaxis = list(title = "t-SNE Dimension 2"),
    zaxis = list(title = "t-SNE Dimension 3")
  )
)

#add tSNE results in SCE
reducedDim(sce,'tsne_3D') <- tsne_map$Y
#add GMM cluster assignments to colData
colData(sce)$cluster <- as.factor(cluster_assignments)

#doing a differential expression of genes based on clusters
#results need to be stored in a list 
cluster_levels <- levels(colData(sce)$cluster)
cluster_DEG <- list()

#make a loop to go through each cluster and perform a one cluster vs all
#loop  through each cell type
for (cluster in cluster_levels){
  #creating group labels
  cell_groups <- factor(ifelse(colData(sce)$cluster == cluster, cluster,'rest of cells'))
  names(cell_groups) <- colnames(cd)
  #build prior
  e_prior1 <- scde.expression.prior(models = error_model, counts = cd)
  #run DE
  g_diff1 <- scde.expression.difference(error_model,cd,e_prior1, groups = cell_groups, n.randomizations = 100, n.cores =1, verbose =1)
  #differental expression analysis
  p_values <- 2*pnorm(abs(g_diff1$Z),lower.tail=F) # 2-tailed p-value
  #adjusting p values to control for false discovery rate (FDR)
  p_values_adj <- 2*pnorm(abs(g_diff1$cZ),lower.tail=F)
  sign_genes <- which(p_values_adj<0.05)
  ord <- order(p_values_adj[sign_genes]) #order by p values
  d <- cbind(g_diff1[sign_genes, 1:3], p_values_adj[sign_genes])[ord,]
  colnames(d) <- c("Lower bound","log2 fold change","Upper bound","p-value")
  #store results
  cluster_DEG[[cluster]]<-d
}

#adding it into the sce object
metadata(sce)$Clusters_DEG <- cluster_DEG
#extracting the top 20 genes from each cluster
top20 <- lapply(metadata(sce)$Clusters_DEG, function(top_df){
  top_df <- top_df[order(top_df[,'p-value'],-abs(top_df[,'log2 fold change'])), ]
  return(head(top_df, 20))
})

#combine all the results into one dataset
top20_table <- do.call(rbind, lapply(names(top20), function(cluster_name){
  cluster_df <- as.data.frame(top20[[cluster_name]])
  cluster_df$cluster <- cluster_name
  cluster_df$gene <- rownames(cluster_df)
  return(cluster_df)
}))

#reordering the columns
column_order <- c('cluster', 'gene', 'log2 fold change', 'p-value', 'Upper bound', 'Lower bound')
top20_table <- top20_table[,column_order]
#sort by cluster, p-value and FDR
top20_table <- top20_table[order(top20_table$cluster, top20_table$'p-value', top20_table$'log2 fold change'), ]
#resetting row names
rownames(top20_table) <- NULL
#saving the results
write.csv(top20_table, "top20_DEGs_per_cluster_original.csv", row.names = FALSE)

#next step involves manually annotating each cluster
#these genes were extracted from PANGOLADB
brain_markers <- list(
  Astrocytes = c("GFAP", "AQP4", "SLC1A2", "ALDH1L1", "GJA1", "S100B",
                 "SLC1A3", "FABP7", "GLUL", "SOX9", "CD44", "SPARCL1"),
  Microglia = c("CX3CR1", "AIF1", "P2RY12", "TMEM119", "CSF1R", "C3",
                "ITGAM", "TREM2", "LY86", "CD68", "CD14", "SPI1"),
  Oligodendrocytes = c("MBP", "PLP1", "MOG", "CLDN11", "MOBP", "CNP",
                       "MAG", "RTN4", "SLC44A1", "OPALIN", "QKI", "TF"),
  Neurons = c("SYT1", "SNAP25", "RBFOX3", "GAD1", "SLC17A7", "NEUROD6",
              "MAP2", "TUBB3", "GRIN1", "CAMK2A", "SYN1", "ELAVL4"),
  Endothelial_cells = c("CLDN5", "FLT1", "PECAM1", "CDH5", "VTN", "ESAM",
                        "VWF", "TIE1", "KDR", "ICAM2", "ENG", "EPHB4"),
  OPCs = c("PDGFRA", "CSPG4", "VCAN", "SOX10", "OLIG2", "GPR17",
           "NKX2-2", "PTPRZ1", "DLL3", "ST8SIA1", "S100A10", "ANXA2"),
  Neuronal_Progenitors = c("ASCL1", "DCX", "SOX2", "NEUROD1", "PROX1", "PAX6",
                           "HES5", "VIM", "TBR2", "EMX2", "FABP7", "NES")
) %>% unlist()

#trim whitespace in both
top20_table$gene <- trimws(top20_table$gene)
brain_markers <- lapply(brain_markers, trimws)

#create a reverse lookup table
marker_ref <- stack(brain_markers) %>% 
  setNames(c("gene", "cell_type"))

#annotate DEGs
annotated_deg <- top20_table %>%
  left_join(marker_ref, by = "gene") %>%
  arrange(cluster)

#view annotated DEGs
head(annotated_deg)

#another step that the authors took was to do supervised clustering
#to do this they mention usage of the mouse data set from Zhang et.al., however their data is unavailable
#from their paper all the top 50 genes were extracted instead
library(stringr)

mouse_markers <- list(
  Astrocytes = str_to_title(c(
    "I-GI", "AAP4", "ITIN3", "BINGT10", "IGAR7", "PLOS4", "GM3", "SIC14A1",
    "PFE1G1", "PILAG3", "CES", "PAER6", "ATL4111", "CM", "COD680", "FINO1",
    "SIC30A10", "SIC6A11", "FPI3", "SIE4A4", "GIGBZ2", "PIP1C1A", "GN1F1",
    "ENIPB2", "EGFR", "ALV64131", "OX1", "NWD1", "ALP13A4", "KUM3", "PHD",
    "SORC2", "TRE", "SOE9", "ABOL2", "FZF10", "LNG1", "MIC1", "CHU1", "ALMO"
  )),
  Neurons = str_to_title(c(
    "REIN", "NEHR2", "SIC17A6", "TPR73", "UNO5", "LINX1AS", "DIXDOI1", "SÄI",
    "S30417022RIK", "MALZ211", "SENG11", "MAG2", "DB1", "TMEM0DA", "ISIL2",
    "IGFOP1", "GAL5", "SIMZ2", "ECE1", "ROBZ2", "DH1AS", "CEL4A", "CELL8",
    "NURH4", "GMZ", "NPY", "TBR1", "SIC32A1", "DB2", "NPAS4", "EB5", "BE11A",
    "CSENS82", "CISLN2", "DOYS5", "VSHM3", "TMEM130", "NSPC", "VGL", "BILHE2"
  )),
  OPCs = str_to_title(c(
    "POJARL", "LINT1", "DON", "FJM", "MIRIQ15", "CLOT1", "SAPOZ2", "KOKK1",
    "RASGF11", "POLTH15'", "CRIMSK1", "DB3", "COL1A2", "FAN70B", "SÄRT1",
    "PRBP", "CASP4", "LGR1", "PIPADEC1A", "NAPH1", "P61", "SILK1A1", "SHC4",
    "SMC1", "EMID1", "RIBP1", "LYND5", "MY1", "GLAS3", "CHGT1", "TMEM179",
    "MEG11", "NCSKL", "SDC3", "FPM", "CASR9D", "GIN1A", "FAM5C", "LNC42"
  )),
  NFO = str_to_title(c(
    "QARBO", "TMEM108", "FJM", "UST", "MOS3", "KIRT1SA", "18100411.15RIK",
    "96300134.20RIK", "SIV3", "PHJO3", "EMP6", "TRIX3", "EMP4", "MC1",
    "COLV3", "TMEM163", "RAG2A", "TMEM2", "OTKA3A", "CYHF2", "FND4A",
    "SIC12A2", "ITPR2", "RNTRZ2", "LIMS2", "SAM4B", "CINZ2", "PIP26A",
    "SIM", "GHP", "RNS2", "FRM2", "SEMLS3A", "FAM5C", "ODC271", "FAN72A",
    "EPO6G", "ALM", "LNC42", "CHU11"
  )),
  MO = str_to_title(c(
    "GOHL", "NIIGRI1", "PPORT15A8", "ADS41", "ASSA", "ACY3", "TRP53NP2",
    "PAL2G16", "ETHE1", "ILGSK1", "HSPH2", "MTP", "HCC2", "NMAT1",
    "CDCD2EP2", "MAI", "MO9", "SIC03A1", "APOD", "GAS", "PALM2", "PNT1A",
    "INR2", "TGGG3", "TBC1C9B", "NO3", "COMPB", "SIC45A3", "CAM1", "OPSLIN",
    "ARSG", "RTN1", "PILS80F", "TRI", "INSC", "CRYPTO", "KI5A8", "TRAK2",
    "CHU11", "BD2A1C"
  )),
  Microglia = str_to_title(c(
    "SIM2", "GAP64", "COR7", "BGZ81D", "TRIF", "CAT15", "OSM", "LINC25",
    "H2-CAA", "CAL63", "CA3", "SLANT8", "CS4", "GRA15", "IT10", "PLAU",
    "CO9", "TMEM119", "CR1A8", "1810011111RIK", "PILAG215", "CUET16",
    "CHZ3B", "COL12", "PLATT", "CO500A", "HI5", "STYL1", "SELE6G", "SASH3",
    "PLP", "TREM2", "TH2", "P2P6", "CO4", "BD2A1A", "FH1", "TIE1", "EMF1"
  )),
  Endothelial = str_to_title(c(
    "CELLS5", "TR", "LY6A", "MAKCAM1", "8430403022RIK", "AKT1614", "LY6C2",
    "MEXX1", "LY6C2", "CAR4", "REG", "AGTER", "SIGHR", "SIC01A4", "SIC18A1",
    "ICAM2", "KAIW3", "SIC19A3", "FAN101B", "SIC16A4", "NOSLIM", "SIGIR",
    "PIGS", "MYCT1", "VWAT", "ANKN37", "SOX18", "PIND", "SERPIN40B", "COS4",
    "PTYRYP1", "SIC35B2", "CABRDB", "FAN126A", "SYMS1", "FH1", "KOXE4",
    "CEB", "STPC3"
  )),
  Pericytes = str_to_title(c(
    "FMOD", "RES2", "IGF2", "OPC3", "OGN", "LINC32", "FINE", "GLB2", "INH2",
    "RBIN10", "BMP6", "ATBR1A2", "POSIM", "SATT1", "LAMC3", "SIC22A6",
    "OKC3B", "SIC43A13", "BBC1", "S100A10", "SENPING1", "COL1A1", "DON",
    "COL1A2", "PROBE", "OPP1B1", "CHEC1", "EMPL", "CEB", "AHNAK", "CO3A1",
    "FSL1", "CO4A5", "VIN", "LAMA2", "MITP4", "KOXE4", "STPC3", "CEB", "AHNAK"
  ))
)

#to continue with the analysis a generic mouse data set from ensembl was loaded
#create the function to convert mouse to human genes
homologs = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
find_orthologs <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (homologs %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (homologs %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}

#apply the function
human_markers <- find_orthologs(unlist(mouse_markers))

#calling all the cell types information from the sce object
cell_types <- unique(colData(sce)$cell_type)
cell_results <- list()
#loop  through each cell type
for (cell_type in cell_types){
  #creating group labels
  cell_groups <- factor(ifelse(colData(sce)$cell_type == cell_type, cell_type,'other'))
  names(cell_groups) <- colnames(cd)
  #build prior
  e_prior1 <- scde.expression.prior(models = error_model, counts = cd)
  #run DE
  g_diff1 <- scde.expression.difference(error_model,cd,e_prior1, groups = cell_groups, n.randomizations = 100, n.cores =1, verbose =1)
  #differental expression analysis
  p_values <- 2*pnorm(abs(g_diff1$Z),lower.tail=F) # 2-tailed p-value
  #adjusting p values to control for false discovery rate (FDR)
  p_values_adj <- 2*pnorm(abs(g_diff1$cZ),lower.tail=F)
  sign_genes <- which(p_values_adj<0.05)
  ord <- order(p_values_adj[sign_genes]) #order by p values
  c <- cbind(g_diff1[sign_genes, 1:3], p_values_adj[sign_genes])[ord,]
  colnames(c) <- c("Lower bound","log2 fold change","Upper bound","p-value")
  #store results
  cell_results[[cell_type]]<-c
}

#storing the DEGs in metadata of sce
metadata(sce)$DEG_results <- cell_results
#storing the log fold changes and p-value in the rowData
gene_names <- rownames(sce)
cells <- names(cell_results)
logFC  <- matrix(NA, nrow = length(gene_names), ncol = length(cells),
                 dimnames = list(gene_names,cells))
#fill the matrix per cell type
for (x in cell_types){
  fc <- cell_results[[x]]
  if (!is.null(fc)){
    #common_genes <- intersect(rownames(fc),gene_names)
    logFC[gene_names,x]<-fc[gene_names,'log2 fold change']
  }
}

#storing it into sce
rowData(sce)$logFC <- logFC

#filter human markers expression from sce object
#extract the average of cells per million
mean_cpm <- rowMeans(assay(sce,'cpm'))

#removing whitespace
rownames(sce) <- trimws(rownames(sce))
human_markers <- trimws(human_markers) %>% unlist()

#keep the markers with mean CPM > 0.1
overlap_genes <- intersect(human_markers, rownames(sce))
markers_filtered <- overlap_genes[mean_cpm[overlap_genes] > 0.1]
markers_filtered #as the output produces only 51 genes the analysis was stopped as clearly the output will not be enough to carry the rest of the analysis

#making a density of plot with the genes identified from the paper
#making a list of all the genes authors identified
genes_of_interest <- c('DAAM2', 'ASPA', 'MAL', 'SEC14L5', 'MAP6D1', 'DPYD', 'PPP1R14A', 'GJB1', 'FA2H', 
                       'MAG', 'CDK18', 'LGI3', 'SHC4', 'UGT8', 'KLK6', 'KCNH8', 'GPR37', 'MOBP', 'LPAR1', 
                       'ADAMTS4', 'ERMN', 'OPALIN', 'CLDN11', 'PLEKHB1', 'GSN', 'GRM3', 'CNP', 'MBP', 'PLP1', 
                       'SLC14A1', 'GLIS3', 'GLI3', 'PPP1R3C', 'CHRDL1', 'CYBRD1', 'CTH', 'SORCS2', 'ITGB4', 'RNF43',
                       'NWD1', 'PAQR6', 'C16orf89', 'ALDH1L1', 'TRIM66', 'HGF', 'CBS', 'ITGA7', 'SLC30A10', 'SLC4A4', 
                       'FGFR3', 'BMPR1B', 'ATP13A4', 'AQP4', 'GPR183', 'CCL4', 'CD83', 'LAPTM5', 'CSF1R', 'HLA-DRA', 
                       'BCL2A1','BGL2A1', 'CD14', 'CCL2', 'APOLD1', 'TM4SF1', 'FLT1', 'A2M', 'PDGFRA', 'LHFPL3', 'MEGF11', 'PCDH15', 
                       'KCNK1', 'KIAA1324', 'LNX1', 'NELL1', 'COBL', 'SLITRK1', 'DPYSL5', 'C14orf37', 'DLX1', 'DLX5', 
                       'GLRA2', 'DLX2', 'DLX5', 'SLC10A4', 'EGFR', 'SST', 'PNOC', 'NXPH1', 'BCL11A', 'DCN', 'TMEM130', 
                       'CNTN4', 'CDO1', 'NFASC','LRRTM3', 'GRIA3', 'RELN')

#etract logFC matrix
logFC_matrix <- rowData(sce)$logFC

#define cell types to exclude
exclude_types <- c("hybrid", "fetal_replicating", "fetal_quiescent")

#get all cell types in your data
all_cell_types <- unique(colData(sce)$cell_type)

#keep only the cell types we want to analyze
cell_types <- setdiff(all_cell_types, exclude_types)

#subset the matrix
logFC_matrix_filtered <- logFC_matrix[, !colnames(logFC_matrix) %in% exclude_types]

#calculate average logFC for each cell type
avg_logFC <- sapply(cell_types, function(ct) {
  rowMeans(logFC_matrix_filtered[, ct, drop = FALSE], na.rm = TRUE)
})

#for genes of interest calculate their average
#when run initially for genes of interest list there was an error message saying that it was out of subscript
#when checked there were three genes which were not found in the matrix, maybe due to them being droped in quality check?
valid_genes <- genes_of_interest[genes_of_interest %in% rownames(avg_logFC)]
goi_logFC <- avg_logFC[valid_genes, ]


#convert to fold change 
fold_changes <-log(2^goi_logFC)

#calculate enrichment relative to background (all other genes)
background_genes <- setdiff(rownames(sce), genes_of_interest)
background_avg <- colMeans(avg_logFC[background_genes, ], na.rm = TRUE)

n_fold_enrichment <- t(t(fold_changes) / log(2^background_avg))

#create the final dataframe
enrichment_df <- as.data.frame(n_fold_enrichment) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, 
               names_to = "cell_type", 
               values_to = "log10_fold_enrichment")

#plotting a density graph
ggplot(enrichment_df, aes(x = log10_fold_enrichment, fill = cell_type, color = cell_type)) +
  geom_density(alpha = 0.3, size = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Log10 fold enrichment",
    y = "density"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) 

#to perform supervised hierarchical clustering, requires to calculate distances between cells again
#authors did not include the fetal cells in this part of the analysis
excluded <- c('fetal_replicating', 'fetal_quiescent')

#create a new filtered sce object
sce_filtered <- sce[, !colData(sce)$cell_type %in% excluded]

#exctracting expression matrix
markers_expression <- assay(sce_filtered[valid_genes,],'logcounts')

#scaling the data
scaled_markers <- t(scale(t(markers_expression)))

#replace any remaining NaNs/Infs with 0
scaled_markers[is.na(scaled_markers)] <- 0
scaled_markers[is.infinite(scaled_markers)] <- 0

#compute distance matrix
markers_distance <- dist(t(scaled_markers), method = 'euclidean')

#perform hierarchical clustering
gene_hclust <- hclust(markers_distance, method = 'ward.D2')

#ploting a dendogram
plot(gene_hclust, labels = FALSE, main = "Hierarchical Clustering of Cells (Markers)")

#assigning a number of clusters to the dendrogram
new_clusters<- cutree(gene_hclust, k = 7)
colData(sce_filtered)$biased_cluster <- factor(new_clusters)

#add this new clusters to the sce object
biased_clusters <- colData(sce_filtered)$biased_cluster

#create annotation dataframe
anno_df <- data.frame(BiasedCluster = as.factor(biased_clusters),
                      UnbiasedCluster = as.factor(colData(sce_filtered)$cluster)
)
rownames(anno_df) <- colnames(scaled_markers)

#define colours for annotations
unbiased_levels <- levels(anno_df$UnbiasedCluster)
biased_levels <- levels(anno_df$BiasedCluster)

type_colors <- structure(hue_pal()(length(unbiased_levels)),
                         names = unbiased_levels)

cluster_colors <- structure(hue_pal()(length(biased_levels)),
                            names = biased_levels)

#build annotation object
obj_annot <- HeatmapAnnotation(
  df = anno_df,col = list(
    BiasedCluster = cluster_colors,
    UnbiasedCluster = type_colors
  ),show_annotation_name = TRUE)

#plot heatmap
Heatmap(scaled_markers,
        name = "Z-score",
        top_annotation = obj_annot,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        column_order = order.dendrogram(as.dendrogram(gene_hclust)),
        col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
        row_names_gp = gpar(fontsize = 8),
        column_title = "Cells",
        row_title = "Marker Genes")


##next steps involve performing PCA on adult and fetal cells
#extract the colData information as a data frame
meta_data <- as.data.frame(colData(sce))

#adding stages of development to the data
meta_data$stage <- dplyr::case_when(
  !is.na(meta_data$age) & grepl('prenatal', meta_data$age, ignore.case = TRUE) &
    meta_data$cell_type %in% 'fetal_replicating' ~ 'Fetal replicating',
  !is.na(meta_data$age) & grepl('prenatal', meta_data$age, ignore.case = TRUE) &
    meta_data$cell_type %in% 'fetal_quiescent' ~ 'Fetal quiescent',
  !is.na(meta_data$age) & grepl('postnatal', meta_data$age, ignore.case = TRUE)
  & grepl('neuron', meta_data$cell_type, ignore.case = TRUE) ~ 'Adult',
  TRUE ~ NA_character_)

#adding stage column in the metadata again
colData(sce)$stage <- meta_data$stage

#subseting the sce again and extract only neuronal cells for these
neurons <- sce[, !is.na(colData(sce)$stage)]


#extracting expression matrix
expr_pca<- as.matrix(logcounts(neurons))
expr_pca <- t(expr_pca)
expr_pca_df <- as.data.frame(expr_pca)

#add group info
expr_pca_df$group <- colData(neurons)$stage
expr_pca_df <- expr_pca_df[!is.na(expr_pca_df$group), ]

#active: Adult + Fetal_Replicating
#subset to only desired groups
expr_pca_df <- expr_pca_df[expr_pca_df$group %in% c("Adult", "Fetal replicating", "Fetal quiescent"), ]
group_labels <- expr_pca_df$group
active_data <- expr_pca_df[, -ncol(expr_pca_df)]

pca_res <- PCA(active_data, 
               scale.unit = TRUE, 
               ncp = 3,
               graph = FALSE)
#extract coordinates
active_coords <- as.data.frame(pca_res$ind$coord)
active_coords$group <- expr_pca_df$group

#plot
plot_ly(active_coords, 
        x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, 
        color = ~group, colors = c("blue", "steelblue", 'pink'),
        type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D PCA: Projecting of neuronal cells",
         scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

#extract the loadings (gene contributions to PCs)
gene_loadings <- as.data.frame(pca_res$var$coord[, 1:3])  # PC1 to PC3

#add gene names as a column
gene_loadings$gene <- rownames(gene_loadings)

#reorder columns if desired
gene_loadings <- gene_loadings[, c("gene", "Dim.1", "Dim.2", "Dim.3")]

#preview the table
head(gene_loadings)

#get top 10 genes contributing to PC1
top_pc1 <- gene_loadings[order(abs(gene_loadings$Dim.1), decreasing = TRUE), ][1:10, ]

#get top 10 genes contributing to PC2
top_pc2 <- gene_loadings[order(abs(gene_loadings$Dim.2), decreasing = TRUE), ][1:10, ]

#get top 10 genes contributing to PC3
top_pc3 <- gene_loadings[order(abs(gene_loadings$Dim.3), decreasing = TRUE), ][1:10, ]

write.csv(gene_loadings, file = "PCA_gene_loadings_original.csv", row.names = FALSE)

#making a minimum spanning tree of all cells#
#authors did a minimum spanning tree only on neuronal cells
#based on this only neuron cells were used and extracted 
neurons_common <- sce_filtered[, !is.na(sce_filtered$cluster) & 
                                 !is.na(sce_filtered$biased_cluster)]

#filter based on neuron-only cells
neurons_common <- neurons_common[, grepl("neuron", neurons_common$cell_type, ignore.case = TRUE)]

#extracting expression matrix
exp_matrix <- assay(neurons_common, 'logcounts')

#compute the distance matrix
neurons_dist <- dist(t(as.matrix(exp_matrix)), method='euclidean')

#converting the distance matrix into a graph object
g <- graph_from_adjacency_matrix(as.matrix(neurons_dist), mode = "undirected", weighted = TRUE, diag = FALSE)

#compute the Minimum Spanning Tree (MST)
mst <- mst(g)

# Detect communities using Walktrap algorithm
wc <- cluster_walktrap(mst)

#find community memberships
colData(neurons_common)$community <- factor(membership(wc))

#set colors for each community
community_colors <- rainbow(length(unique(membership(wc))))
V(mst)$color <- community_colors[membership(wc)]

#plot MST plot
plot(
  mst,
  vertex.size = 4,
  vertex.label = NA,
  edge.width = 1,
  layout = layout_with_kk(mst)
)

#saving the graph
png("neuronal_mst_plot_original.png", width = 800, height = 800, res = 150)

dev.off()
