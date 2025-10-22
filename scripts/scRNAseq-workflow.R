if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
setwd("~/zebrafish-development/data/TDR70")

BiocManager::install(c('scuttle', 'scran', 'scater', 'uwot', 'rtracklayer', 'scRNAseq', 'monocle', 'Seurat'))

# Load packages
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

file_move("~/zebrafish-development/data/TDR70/UMAP_70.png", "~/zebrafish-development/UMAP")




