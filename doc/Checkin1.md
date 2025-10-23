# How you’ve addressed prior feedback 

Most of our initial meeting was regarding steps to take given our general idea of doing RNAseq to find tissue differences in developmental gene expressions. Specifically, we got the suggestion to pick a certain timepoint and start there, first by replicating the Zebrahub clustering and then finding gene differences. We also got a recommendation to use Bioconductor to set up our workflow. Because of this, a lot of our work this time thus went into just setting up a workflow using the recommended Bioconductor and testing our pipeline on an example of the Zebrahub dataset (our example is TDR67).

# New progress since last submission

Our big piece of progress was getting a pipeline that worked! We normalized and scaled the data, then ran a PCA, found neighbors, found clustered, and got a working UMAP. Our test dataset (TDR67) had 18 clusters using our clustering method. This pipeline should work for any dataset we download from Zebrahub now, so this is really great for starting to replicate their clustering method. We have run our pipeline on all four of the fish (TDR67-70) found in our chosen timepoint at 3 dpf and generated 4 UMAPs (which we have uploaded). 

# Project Organization

So far organization has not been a major problem. As we download more datasets and adjust our scripts keeping our folders organized and our GitHub up to date will be things we continue doing. We have a separate scripts and raw data folder to help keep things neat.

# Struggles you are encountering and questions you would like advice on

The big issue we're currently having is regarding annotations. Now that we have a working pipeline, replicating what the authors did as well as finding genes in different tissues would be made a lot easier by having their annotations to reference. However, we cannot find them on the Zebrahub website, and the data is only provided in a .h5 format. Advice on where to find their annotations or for other possible references we use would help us a lot with our next steps.

# Figures

Figure 1: UMAP of TDR67
![alt text](https://github.com/dgardon26/zebrafish-development/blob/main/UMAP/UMAP_67.png)

Figure 2: UMAP of TDR68
![alt text](https://github.com/dgardon26/zebrafish-development/blob/main/UMAP/UMAP_68.png)

Figure 3: UMAP of TDR69
![alt text](https://github.com/dgardon26/zebrafish-development/blob/main/UMAP/UMAP_69.png)

Figure 4: UMAP of TDR70
![alt text](https://github.com/dgardon26/zebrafish-development/blob/main/UMAP/UMAP_70.png)
