# zebrafish-development

## Title: Transcriptomic Cross Tissue Analysis of Key Developmental Genes in Zebrafish

## Proposal:

We want use the data from Zebrahug to look at expression of groups of major developmental molecules across zebrafish embryo tissues, starting with replicating their tissue clustering to confirm our pipeline works and then examining our molecules of choice from there. Specifically, we were thinking of looking at Hedgehog, Wnt, Notch, TGF-B, Hox and/or FGF, all major developmental pathways. Additionally, we may do a timepoint comparison of one gene across one tissue, and expand a timepoint across tissues if time permits and one is highly expressed, or compare results between datasets.


## Example Published Figure: 
![Example of clustering](/gr1_lrg.jpg)
![Example of heat map](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41586-023-06720-2/MediaObjects/41586_2023_6720_Fig5_HTML.png?as=webp)

## Datasets: 

Primary Dataset https://www.cell.com/cell/fulltext/S0092-8674(24)01147-4#fig1 

Secondary Dataset (for stretch goals) https://www.nature.com/articles/s41586-023-06720-2 

## Software/Tools:

Bioconductor, version 1.18.0, https://bioconductor.org/books/release/OSCA/

Potential Additional Tools:
https://github.com/chanzuckerberg/cellxgene 
https://github.com/czbiohub-sf/zebrahub_analysis

## Goals: 

1) Pick one gene and do a cross sectional tissue comparison within an embryonic zebrafish dataset. Represent this with clustering and a heat map. From there could compare multiple genes of interest.

2) (Stretch) Do a timepoint comparison of a gene within a tissue set and potentially expand at a certain timepoint across tissues if gene is of interest (highly expressed).

3) (Strech-Strech) Use another zebrafish data set and compare results between data sets.
