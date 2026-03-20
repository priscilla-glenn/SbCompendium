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
### Git LFS Requirement for Large Files

This repository includes large files tracked with **Git LFS**.  
**Before cloning**, make sure Git LFS is installed:

#### Install Git LFS
- **Windows / macOS / Linux:**
  - https://git-lfs.com/  
- Or via command line:
  - git lfs install

## Clone the repo with LFS (cpmmand line)
git clone https://github.com/priscilla-glenn/SbCompendium.git

cd SbCompendium

git lfs pull

In Rstudio, once the repo is cloned: change the path to the folder where you cloned the repo
- devtools::install_local("path/to/SbCompendium", build_vignettes = FALSE)

library(SbCompendium)

# How to load the compendium and an example dataset
data("example_data")

experiment <- importTable(example_data[["nodal_buds_28"]])

head(experiment)


# To load the full compendium dataset
data("sorghum_compendium")
- it will load into your environment as sorghum_compendium and each experiment can be viewed. ex: sorghum_compendium$nodal_buds_28

#Please check the vignette and manual to see full details and usage for each function
