if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
setwd("~/zebrafish-development/data/TDR70")

BiocManager::install(c('scuttle', 'scran', 'scater', 'uwot', 'rtracklayer', 'scRNAseq', 'monocle', 'Seurat'))

# Load packages
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scRNAseq)
library(monocle)
library(fs)

install.packages("hdf5r")

# Create a Seurat object

test <- Read10X_h5("~/zebrafish-development/data/TDR70/filtered_feature_bc_matrix.h5")
data <- CreateSeuratObject(counts = test)
data.updated = UpdateSeuratObject(object = data)
data_counts <- as.matrix(data.updated@assays$RNA@layers$counts)
data <- CreateSeuratObject(counts = data_counts)

# Quality control and filtering
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Normalize the data
# This step adjusts for differences in sequencing depth across cells.
data <- NormalizeData(data)
# Identify variable features
# Variable features are genes that exhibit high variability across cells.
# These features are used in downstream analyses to focus on informative genes.
data <- FindVariableFeatures(data)

# Scale the data
# Scaling centers and scales each gene, making them comparable for PCA and clustering.
# This is an important step for most dimensionality reduction techniques.
data <- ScaleData(data)

# Step 4: Run Principal Component Analysis (PCA)
# PCA reduces the dimensionality of the data, allowing us to identify major sources of variation.
# Here, we use the variable features identified earlier to perform PCA.
data <- RunPCA(data)

# Step 5: Find Neighbors
# This step identifies nearest neighbors for each cell based on their PCA scores.
# It is an essential step before clustering the cells.
data <- FindNeighbors(data, dims = 1:10)

# Step 6: Cluster the cells
# Clustering groups cells with similar expression profiles, aiding in cell-type identification.
data <- FindClusters(data, resolution = 0.5)

# Step 7: Run UMAP for visualization
# UMAP (Uniform Manifold Approximation and Projection) is a popular method for visualizing high-dimensional data.
# This step reduces the data to two dimensions, making it easier to visualize clusters.
data <- RunUMAP(data, dims = 1:10)

# Step 8: Plot the UMAP results
# This visualization shows the clusters identified by Seurat, each in a different color.
UMAP_70 <- DimPlot(data, reduction = "umap")

ggsave("UMAP_70.png", plot = UMAP_70)

#repeat for each umap file, we started with samples 67-70.

file_move("~/zebrafish-development/data/TDR70/UMAP_70.png", "~/zebrafish-development/UMAP")

#repeat for each umap file, we started with samples 67-70, just change the number in each step.


##Next step: characterize clusters by enriched gene markers
data.markers <- FindAllMarkers(data, only.pos = TRUE)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#compare feature# by gene# in gene_names list -> look for biomarkers in each tissue in literature and in zebrafish website
install.packages('devtools')
devtools::install_github('immunogenomics/presto')

gene_names <- data.updated@assays$RNA@features
gene_names <- as.data.frame(rownames(gene_names))
gene_names$feature <- paste0("Feature", seq(1, nrow(gene_names)))
gene_names <- merge(gene_names, data.markers, by.x = "feature", by.y = "gene", all.x = TRUE)
data.markers <- merge(data.markers, gene_names, by.x = "gene", by.y = "feature", all.x = TRUE)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

gene_names <- as_tibble(gene_names)