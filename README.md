## Title: Transcriptomic Cross-Tissue Analysis of Key Developmental Genes in Zebrafish

## Proposal:

We will use scRNA-seq data from Zebrahub to look at expression of groups of major developmental molecules across zebrafish embryo tissues, starting with replicating their tissue clustering to confirm our pipeline works and then examining our molecules of choice from there. Specifically, we were thinking of looking at Hedgehog, Wnt, Notch, TGF-B, Hox and/or FGF, all major developmental pathways. Additionally, we may do a timepoint comparison of one gene across one tissue, and expand a timepoint across tissues if time permits and one is highly expressed, or compare results between datasets.


## Example Published Figure: 
![Example of clustering](/gr1_lrg.jpg)
![Example of heat map](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41586-023-06720-2/MediaObjects/41586_2023_6720_Fig5_HTML.png?as=webp)

## Datasets: 

Primary Dataset https://zebrahub.sf.czbiohub.org/data ->
https://drive.google.com/drive/folders/1WQNr8RkiUqxm6m_SAIXKedPXBxhZ6FWu
Looking at h5 count matrices of 4 zebrafish samples at 3 days post fertilization (3dpf). Samples TDR67 - TDR70.


Secondary Dataset (for stretch goals) https://www.nature.com/articles/s41586-023-06720-2 

## Software/Tools:

R, version 4.5.1, https://posit.co/download/rstudio-desktop/

Bioconductor, version 1.18.0, https://bioconductor.org/books/release/OSCA/

Seurat, version 5, https://satijalab.org/seurat/

SingleCellExperiment, version 1.30.1, https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html

scater, version 1.36.0, https://www.bioconductor.org/packages/release/bioc/html/scater.html

scran, version 1.36.0, https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html

scuttle, version 1.18.0, https://bioconductor.org/packages/release/bioc/html/scuttle.html

uwot, version 0.2.3, https://cran.r-project.org/web/packages/uwot/index.html

monocle, version 2.36.0, https://www.bioconductor.org/packages/release/bioc/html/monocle.html

scRNAseq, version 2.22.0, https://www.bioconductor.org/packages/release/data/experiment/html/scRNAseq.html

Potential Additional Tools:
https://github.com/chanzuckerberg/cellxgene 
https://github.com/czbiohub-sf/zebrahub_analysis

## Goals: 

1) Pick one gene and do a cross-sectional tissue comparison within an embryonic zebrafish dataset. Represent this with clustering and a heat map. From there could compare multiple genes of interest.

2) (Stretch) Do a timepoint comparison of a gene within a tissue set and potentially expand at a certain timepoint across tissues if gene is of interest (highly expressed).

3) (Stretch-Stretch) Use another zebrafish data set and compare results between data sets.

## Target Genes:
### Notch Signalling
Notch1a - notch receptor 1a\
Notch1b - notch receptor 1b\
Notch2 - notch receptor 2\
Notch3 - notch receptor 3\
dla - deltaA\
dlb - deltaB\
jag1a - jagged canonical Notch ligand 1a\
jag1b - jagged canonical Notch ligand 1b\
jag2a - jagged canonical Notch ligand 2a\
jag2b - jagged canonical Notch ligand 2b\
jag3 - jagged canonical Notch ligand 3\

### BMP Signaling
bmpr1aa	- bone morphogenetic protein receptor, type IAa\
bmpr1ab	-	bone morphogenetic protein receptor, type IAb\
bmpr1ba	-	bone morphogenetic protein receptor, type IBa\
bmpr1bb	-	bone morphogenetic protein receptor, type IBb\
bmp1l	- tolloid	bone morphogenetic protein 1, like\
bmp10l	-	bone morphogenetic protein 10, like\
bmp10	-	bone morphogenetic protein 10\
bmp2a	-	bone morphogenetic protein 2a\
bmp2b	-	bone morphogenetic protein 2b\
bmp4	-	bone morphogenetic protein 4\
bmp5	-	bone morphogenetic protein 5\
bmp7a	-	bone morphogenetic protein 7a\
bmp8a	-	bone morphogenetic protein 8a\
brinp3a.1	-	bone morphogenetic protein/retinoic acid inducible neural-specific 3a, tandem duplicate 1\
brinp3a.2	-	bone morphogenetic protein/retinoic acid inducible neural-specific 3a, tandem duplicate 2\
brinp3b	-	bone morphogenetic protein/retinoic acid inducible neural-specific 3b\
brinp2	-	bone morphogenetic protein/retinoic acid inducible neural-specific 2\
bmp15	-	bone morphogenetic protein 15\
bmp1a	-	bone morphogenetic protein 1a\
bmp3	-	bone morphogenetic protein 3\
bmp6	-	bone morphogenetic protein 6\
bmp7b	-	bone morphogenetic protein 7b\
brinp1	-	bone morphogenetic protein/retinoic acid inducible neural-specific 1\

### Others: FGF, Wnt, Shh
