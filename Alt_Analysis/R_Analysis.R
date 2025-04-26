#Alternative clustering pipeline for scRNA-seq analysis of the Darmanis et al. (2015) dataset

#loading necessary libraries
library(Seurat)
library(GEOquery)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shiny)
library(clustree)

#setting working directory
setwd('/home/hf180/Desktop/scRNAseq/scRNAseq-GroupB/finished_GeneMatrices/alternative_GeneCounts')

#loading data into a data frame
results <- read.csv("alt_results.csv", row.names = 1)
matrix <- as.matrix(results)

#creating a Seurat object
seurat_matrix <- CreateSeuratObject(counts = matrix, project = "scRNAseq")

#checking the Seurat object
str(seurat_matrix)
head(colnames(seurat_matrix))  # GSM and SRR IDs (sample names)
head(rownames(seurat_matrix))  # Gene names
dim(seurat_matrix)

#cleaning row and column names
rownames(seurat_matrix) <- trimws(rownames(seurat_matrix))
colnames(seurat_matrix) <- sub("^X\\.", "", colnames(seurat_matrix))

#taking first 466 cells
cells_to_keep <- colnames(seurat_matrix)[1:466]
seurat_matrix <- subset(seurat_matrix, cells = cells_to_keep)

#visualizing QC metrics
VlnPlot(seurat_matrix, features = c("nFeature_RNA", "nCount_RNA"), ncol=3)
FeatureScatter(seurat_matrix, feature = "nCount_RNA", feature2 = "nFeature_RNA")

#quality control - filtering low quality cells 
seurat_matrix <- subset(seurat_matrix, subset = nCount_RNA < 16000 & nFeature_RNA < 10000)
seurat_matrix <- NormalizeData(seurat_matrix) #default "LogNormalize"

#plotting top 10 highly variable genes
seurat_matrix <- FindVariableFeatures(seurat_matrix)
top10 <- head(VariableFeatures(seurat_matrix),10)
plot1 <- VariableFeaturePlot(seurat_matrix)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling
seurat_matrix <- ScaleData(seurat_matrix, features = rownames(seurat_matrix))

#PCA - reducing noise from less variable genes for clustering
seurat_matrix <- RunPCA(seurat_matrix, features = VariableFeatures(object = seurat_matrix))

#Plots for visualization after PCA
ElbowPlot(seurat_matrix) 
DimPlot(seurat_matrix, reduction = "pca")
DimHeatmap(seurat_matrix, dims = 1:10, balanced = TRUE)

#Pair-wise distancing
pca_embeddings <- Embeddings(seurat_matrix, "pca")[, 1:10] #extracts coordinates of cells in PCA
pairwise_distances <- dist(pca_embeddings) #pairwise distances calculated (Euclidean distance) between cells based on PCA embedding
pairwise_distances_matrix <- as.matrix(pairwise_distances)
heatmap(pairwise_distances_matrix, main = "Pairwise Distance Heatmap") #closer cells have similar gene expression profiles

#plotting pairwise distances with Plotly
heatmap_plot <- plot_ly(
  x = colnames(pairwise_distances_matrix),
  y = rownames(pairwise_distances_matrix),
  z = pairwise_distances_matrix,
  type = "heatmap",
  colors = colorRampPalette(c("blue", "red"))(100)
)
heatmap_plot

#reducing noise after PCA by identifying statistically significant PCs (permutation)
seurat_matrix <- JackStraw(seurat_matrix, num.replicate = 100)
seurat_matrix <- ScoreJackStraw(seurat_matrix, dims = 1:20)
JackStrawPlot(seurat_matrix, dims = 1:20)

#unbiased clustering for multiple resolutions
seurat_matrix <- FindNeighbors(seurat_matrix, dims = 1:10)
seurat_matrix <- FindClusters(seurat_matrix, resolution = c(0.4, 0.6, 0.8, 1.0))
clustree(seurat_matrix)

#UMAP visualisation of unbiased clusters
seurat_matrix <- RunUMAP(seurat_matrix, dims = 1:10)
umap_plot <- DimPlot(seurat_matrix, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("Seurat Clusters") + theme(plot.title = element_text(size = 20, face = "bold"))
umap_plot <- ggplotly(umap_plot)
umap_plot

#all genes in the clusters
unbiased_markers <- FindAllMarkers(seurat_matrix, only.pos = TRUE)
head(unbiased_markers, 50)

#group by the markers by 'cluster' column
grouped_markers <- group_by(unbiased_markers, cluster)

#select top 5 genes with highest avg_log2FC for each cluster
top5 <- top_n(grouped_markers, 5, avg_log2FC)

#feature plot for top genes across all clusters - too compact
#FeaturePlot(seurat_matrix, features = top_genes_per_cluster)

plots <- list() #storing feature plots for each cluster as a list

#loop over each unique cluster and generate a separate FeaturePlot for each cluster
for (cluster in unique(top5$cluster)) {
  #get top 5 genes for current cluster
  cluster_genes <- top5[top5$cluster == cluster, "gene", drop = TRUE]
  #feature plots
  plots[[as.character(cluster)]] <- FeaturePlot(seurat_matrix, features = cluster_genes) +
    plot_annotation(title = paste("Cluster", cluster),
                    theme = theme(plot.title = element_text(size = 20, face = "bold")))
}

#view feature plots of each cluster
print(plots)

#loading soft file downloaded from GEO
gse_soft <- getGEO(filename = "GSE67835_family.soft", GSEMatrix = TRUE)

#extract GSM metadata
gsm_list <- GSMList(gse_soft)

#initialising an empty data frame to store metadata
metadata_df <- data.frame(GSM = character(), title = character(), cell_type = character(), age = character(), stringsAsFactors = FALSE)

#loop through each GSM in the gsm_list
for (i in 1:length(gsm_list)) {
  #get metadata for current GSM
  md <- Meta(gsm_list[[i]])
  #extract characteristics
  characteristics <- md$characteristics_ch1
  #extract 'cell_type' and 'age' from characteristics
  cell_type <- sub("cell type:\\s*", "", grep("^cell type:", characteristics, value = TRUE)[1])
  age <- sub("age:\\s*", "", grep("^age:", characteristics, value = TRUE)[1])
  
  #create a data frame with current GSM metadata
  gsm_metadata <- data.frame(
    GSM = names(gsm_list)[i],
    title = md$title,
    cell_type = cell_type,
    age = age,
    stringsAsFactors = FALSE
  )
  
  #add the current GSM's metadata to the main data frame
  metadata_df <- rbind(metadata_df, gsm_metadata)
}

#view meta data
head(metadata_df)
table(metadata_df$cell_type)

#extracting GSM ID from cell name
gsm_ids <- sapply(strsplit(colnames(seurat_matrix), "_"), `[`, 1)

#creating a new cell_type _ector based on GSM match
cell_type_vector <- metadata_df$cell_type[match(gsm_ids, metadata_df$GSM)]

#adding cell_type to seurat object metadata
seurat_matrix$cell_type <- cell_type_vector

#UMAP by cell type
umap_by_cell_type <- DimPlot(seurat_matrix, group.by = "cell_type", label = TRUE, repel = TRUE)
umap_by_cell_type <- ggplotly(umap_by_cell_type)
umap_by_cell_type

#creating a new column to store numeric age
metadata_df$age_numeric <- as.numeric(gsub(".*?(\\d+).*", "\\1", metadata_df$age))

#checking if changes are reflected
head(metadata_df$age_numeric)

#assigning changes to seurat object
age_numeric_vector <- metadata_df$age_numeric[match(gsm_ids, metadata_df$GSM)]
seurat_matrix$age_numeric <- age_numeric_vector

head(seurat_matrix@meta.data)

#subsetting cell types into adult, fetal_quiescent, and fetal_replicating
#creating a new column for 'grouped_cell_type'
seurat_matrix$grouped_cell_type <- ifelse(
  seurat_matrix$cell_type %in% c("astrocytes", "endothelial", "hybrid", "microglia", "neurons", "oligodendrocytes", "OPC"),
  "Adult", 
  ifelse(seurat_matrix$cell_type == "fetal_quiescent", "Fetal_quiescent", "Fetal_replicating")
)

#normalisation, scaling, PCA
seurat_matrix <- NormalizeData(seurat_matrix)
seurat_matrix <- FindVariableFeatures(seurat_matrix)
seurat_matrix <- ScaleData(seurat_matrix)
seurat_matrix <- RunPCA(seurat_matrix, npcs = 10)

#UMAP plot coloured by grouped cell type
seurat_matrix <- RunUMAP(seurat_matrix, dims = 1:10)
umap_grouped_cell_type <- DimPlot(seurat_matrix, reduction = "umap", group.by = "grouped_cell_type")
umap_grouped_cell_type <- ggplotly(umap_grouped_cell_type)
umap_grouped_cell_type

#distribution of cell types in each cluster
confusion_matrix <- table(seurat_matrix$seurat_clusters, seurat_matrix$cell_type)
write.csv(as.data.frame(confusion_matrix), file = "cluster_vs_celltype.csv", row.names = FALSE)
#convert the confusion matrix to a data frame
confusion_matrix_df <- as.data.frame(confusion_matrix)

#bar plot
bar_plot_cell_type <- plot_ly(
  data = confusion_matrix_df,
  x = ~Cluster,
  y = ~Freq,
  color = ~CellType,
  type = "bar"
) %>%
  layout(
    title = "Distribution of Cell Types in Each Cluster",
    xaxis = list(title = "Cluster"),
    yaxis = list(title = "Cell Count")
  )
bar_plot_cell_type

#agreement between seurat clusters and darmanis clusters
agreement_plot <- plot_ly(
  data = confusion_matrix_df,
  x = ~CellType,
  y = ~Freq,
  color = ~as.factor(Cluster),
  type = "bar",
  colors = "Set2"
) %>%
  layout(
    title = "Agreement Between Cell-Type Assignment",
    xaxis = list(title = "Biased Cell Type"),
    yaxis = list(title = "Cell Count")
  )
agreement_plot
