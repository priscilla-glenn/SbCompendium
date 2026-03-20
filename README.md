# SbCompendium

## Purpose
This package was created to maintain the Mullet Lab Sorghum RNA-seq compendium.
A manuscript with the meta data for each experiment is currently under review.
There are currently 40 experiments within this R package

This package allows users to analyze expression data from the RNA-seq compendium with the following functions.

- Importing transcript-level data from the Mullet Sorghum RNA-seq datasets and aggregating transcripts to genes
- UMAP quality control of sample relationships
- Differential expression analysis using DESeq2
- k-means clustering of expression patterns
- Exporting heatmap tables with optional annotations and clustering results

## Installation
# Install from GitHub
install.packages("remotes")
remotes::install_github("priscilla-glenn/SbCompendium")
library(SbCompendium)

# How to load the compendium and an example dataset
data("sorghum_compendium")
experiment <- importTable(sorghum_compendium[["cle_33"]])
head(experiment)

#Please check the vignette and manual to see full details and usage for each function
