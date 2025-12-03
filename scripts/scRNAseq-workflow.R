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

gene_names <- as_tibble(gene_names)

data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(data, features = top10$gene) + NoLegend()
specific_gene <- gene_names %>% filter(grepl("Feature7096", feature))
print(specific_gene) # Feature3539
shha_feature <- gene_names %>% filter(`rownames(gene_names)` == "shha") %>% pull(feature)
print(shha_feature)
wnt1_feature <- gene_names %>% filter(`rownames(gene_names)` == "wnt1") %>% pull(feature)
print(wnt1_feature) # Feature25030
notch1a_feature <- gene_names %>% filter(`rownames(gene_names)` == "notch1a") %>% pull(feature)
print(notch1a_feature) # Feature25656
tgfb1a_feature <- gene_names %>% filter(`rownames(gene_names)` == "tgfb1a") %>% pull(feature)
print(tgfb1a_feature) # Feature23350
hoxb2a_feature <- gene_names %>% filter(`rownames(gene_names)` == "hoxb2a") %>% pull(feature)
print(hoxb2a_feature) # Feature6413
hoxb1b_feature <- gene_names %>% filter(`rownames(gene_names)` == "hoxb1b") %>% pull(feature)
print(hoxb1b_feature) # Feature21631
fgf1a_feature <- gene_names %>% filter(`rownames(gene_names)` == "fgf1a") %>% pull(feature)
print(fgf1a_feature) # Feature18373
fgfr2_feature <- gene_names %>% filter(`rownames(gene_names)` == "fgfr2") %>% pull(feature)
print(fgfr2_feature) # Feature19918

data_test.markers <- FindAllMarkers(data, only.pos = TRUE)
data_test.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

data_test.markers <- filter(data.markers, p_val_adj < 0.05)
data_test.markers <- filter(data_test.markers, avg_log2FC > 2)
cluster1.markers <- filter(data_test.markers, cluster == 1) # Feature7096 elavl3
cluster2.markers <- filter(data.markers, cluster == 2) # Feature25709 col5a1
cluster3.markers <- filter(data.markers, cluster == 3) # Feature29076 syt5b
cluster4.markers <- filter(data.markers, cluster == 4) # Feature7240 hbbe1.3
cluster5.markers <- filter(data.markers, cluster == 5) # Feature17099 fabp7a
cluster6.markers <- filter(data.markers, cluster == 6) # Feature22116 rplp1, this is a pretty weak marker though
cluster7.markers <- filter(data.markers, cluster == 7) # Feature8239 krt4
cluster8.markers <- filter(data.markers, cluster == 8) # Feature25460 pfn1
cluster9.markers <- filter(data.markers, cluster == 9) # Feature21113 gngt2b
cluster10.markers <- filter(data.markers, cluster == 10) # Feature7242 hbae3
cluster11.markers <- filter(data.markers, cluster == 11) # Feature16779 actc1b 
cluster12.markers <- filter(data.markers, cluster == 12) # Feature20205 CR318588.4 this one is not amazing
cluster13.markers <- filter(data.markers, cluster == 13) # Feature17099 fabp7a
cluster14.markers <- filter(data.markers, cluster == 14) # Feature23221 apoa1b 
cluster15.markers <- filter(data.markers, cluster == 15) # Feature7452 lgals2b 
cluster16.markers <- filter(data.markers, cluster == 16) # Feature13536 prss59.2
cluster17.markers <- filter(data.markers, cluster == 17) # Feature6386 col1a1a
cluster18.markers <- filter(data.markers, cluster == 18) # Feature12614 atp1a1b
cluster19.markers <- filter(data.markers, cluster == 19) # Feature3603 cdh5 
cluster20.markers <- filter(data.markers, cluster == 20) # Feature28325 and3 
cluster21.markers <- filter(data.markers, cluster == 21) # Feature13327 abcb5
cluster22.markers <- filter(data.markers, cluster == 22) # Feature25002 si:dkey-205h13.2

top10_test <- left_join(top10, gene_names, by = c("gene" = "feature")) 
top10_test <- top10_test %>% select(-gene) #use this

FeaturePlot(data, features = c("Feature25709"))

library(readxl)
reference_markers <- read_excel("~/Downloads/1-s2.0-S0012160619304919-mmc7.xlsx")
reference_markers2 <- read_excel("~/Downloads/mmc2.xlsx") %>% filter(timepoint == "3dpf")
marker_genes <- reference_markers2 %>% select(annotation_human_readable, `marker genes`)

#what i found so far? just based on manual comparisons between the marker sets and ours
#cluster 2 cranial neural crest
#cluster 3 CNS - neuron
#cluster 5 CNS - neuron // radial glia
#cluster 15 intestine
#cluster 10 blood
#cluster 20 fin basal cell // epidermis
#cluster 21 periderm // integument


