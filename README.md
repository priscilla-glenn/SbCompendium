# SbCompendium

If you use the SbCompendium package/data, please cite our paper, "Compilation and Utilization of a Sorghum Transcriptome Compendium for Gene Regulatory Network Analysis and Crop Trait Engineering". https://doi.org/10.1111/tpj.70922 

If you have any questions about the datasets themselves, please email priscilla.doupnik@ag.tamu.edu.

If you have any bugs/questions about this R package, please open a new issue and we will respond in a timely manner. 

## Purpose
This package was created to maintain the Mullet Lab Sorghum RNA-seq compendium.
A manuscript with the meta data for each experiment is currently under review.
There are currently 38 experiments within this R package

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

## Clone the repo with LFS (command line)
git clone https://github.com/priscilla-glenn/SbCompendium.git

cd SbCompendium

git lfs pull
```

In RStudio, once the repository is cloned, change the path below to the folder
where you cloned it:

```r
devtools::install_local("path/to/SbCompendium", build_vignettes = FALSE)
library(SbCompendium)
```

## How to load the compendium and an example dataset

```r
data("example_data")

experiment <- importTable(example_data[["nodal_buds_28"]])

head(experiment)
```


## To load the full compendium dataset

```r
data("sorghum_compendium")
```

This creates `sorghum_compendium` in your environment. Access an experiment,
for example, with `sorghum_compendium$nodal_buds_28`. Use
`list_compendium(sorghum_compendium)` to list all experiments or supply
`pattern` to search their names.

## Compendium metadata

Metadata for the compendium is provided in the Excel workbook
[Compendium_meta_data_08-31-2026.xlsx](inst/extdata/Compendium_meta_data_08-31-2026.xlsx).
It is distributed with the package as a separate file, rather than embedded in
`sorghum_compendium`. You can open the workbook in Excel or locate the installed
copy from R:

```r
metadata_file <- system.file(
  "extdata", "Compendium_meta_data_08-31-2026.xlsx",
  package = "SbCompendium", mustWork = TRUE
)
metadata_file
```

See the vignette's **Compendium metadata** section for instructions on reading
the workbook into R.

# Please check the vignette and manual to see full details and usage for each function
Experiments in the compendium include:

- TX08001_coreRind_35
- TX08001_highRes_34
- nitrogen_33
- roots_field_water_stress_32
- wray_LCM_31
- della_LCM_30
- roots_GH_water_stress_29
- nodal_buds_28
- LS_removal_27
- TX08001_density_26
- seedling_MLG_25
- x58M_BSI_23.1
- x100M_BSI_23.2
- tx_sas_22
- wray_stem_dev_21
- DellaDevV3_20
- DYM_19
- DDYM_18
- TX08001_radial_17
- TX08001_radial_17
- nroots_16
- SYM_15
- temp_12_and_13*
- T11_R07020_AER
- cle_10
- SASV3_TPM_9
- dleaves_8
- DellaDiurnalCycling_7.1
- KellerDiurnalCycling_7.2*
- Dw2dw2_6
- internode_growth_5
- ccBTx623_4
- cc100M_3.1
- ccSM100_3.2*
- buds_2
- BTx623_atlas_1

*These experiments do not include count data for DEseq analysis
