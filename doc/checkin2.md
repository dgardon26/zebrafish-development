# How you’ve addressed prior feedback 

We've looked to generate some genes to identify our clusters using AUC analysis. We will compare these
with reference datasets to try to validate our results.
Reference datasets so far: 
1. Sur, A., Wang, Y., Capar, P., Margolin, G., Prochaska, M. K., & Farrell, J. A. (2023). Single-cell analysis of shared signatures and transcriptional diversity during zebrafish development. Developmental Cell, 58(24), 3028-3047.e12. https://doi.org/10.1016/j.devcel.2023.11.001
2. Farnsworth, D. R., Saunders, L. M., & Miller, A. C. (2020). A single-cell transcriptome atlas for zebrafish development. Developmental Biology, 459(2), 100-108. https://doi.org/10.1016/j.ydbio.2019.11.008

# New progress since last submission

As mentioned above, we extracted individual genes within the dataset that could be used as markers to separate our clusters
and potentially categorize tissue, and we will compare these with reference datasets. Using these markers, we were able to generate a heatmap showing the differential expression of the top 10 markers for each cluster.
<img width="2602" height="1252" alt="image" src="https://github.com/user-attachments/assets/6857d734-be3c-48cc-ba23-6b4c20e3c4f0" />

Additionally, we were able to visualize top genes for each cluster and where they are located in our PCA plot. Once clusters have been labeled, we can use this same approach to look for which developmental genes we are interested in to see which clusters they are enriched in. 
Example: Feature25709 (col5a1) is the top marker for cluster 2, and we can see it is highly enriched there. 

<img width="1476" height="1050" alt="image" src="https://github.com/user-attachments/assets/f5143c26-eb87-4d41-b88f-f7d972dc4d1f" />
<img width="2100" height="2100" alt="image" src="https://github.com/user-attachments/assets/0cb39c19-1477-4805-a15c-bf8128f4d3d8" />


# Project Organization

We've kept our same system since the last check-in. I think our only caveat right now is how to organize any reference datasets
we use, but otherwise anything else we generate has been relatively easy to keep organized within our current directories.

# Struggles you are encountering and questions you would like advice on

This was an unsupervised clustering method based on the fact that we did not have any annotations in the original dataset,
but one issue we're running into now is that a few of the clusters do not "look good" even when using their strongest marker per
AUC calculations. For example, clusters 0, 1, and 6 have very weak markers, and you can see this reflected in the heatmap. This will make it difficult to characterize those clusters.  I think this will be dependent on pursuing further genes in our AUC/heatmap and assessing comparisons to reference datasets. Additionally, most papers we have seen that characterize scRNA-seq clusters have a different # of clusters than ours, and typically have some sort of expert in zebrafish development to confirm that those tissues are correct, which is something we lack. 

Should we try to recluster with the same number of clusters as those from our reference papers to try and make the process simpler? We could also try to map to the same # of clusters from the zebrahub maps.

