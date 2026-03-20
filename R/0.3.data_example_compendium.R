#' Example sorghum compendium dataset
#'
#' A small example RNA-seq expression dataset containing transcript-level
#' TPM, TPMmean, and count data for multiple samples.
#'
#' The dataset is designed to demonstrate the full SbCompendium workflow,
#' including:
#' - `importTable()` transcript-to-gene aggregation
#' - `build_deseq()` differential expression analysis
#' - k-means clustering
#' - heatmap export
#' - UMAP visualization
#'
#' Columns include:
#' - `TranscriptID`
#' - replicate-level `_TPM`
#' - replicate-level `_counts`
#' - condition-level `_TPMmean`
#'
#' @format A data.frame with 40 rows (transcripts) and expression columns.
#'
#' @source Simulated data for package examples.
"example_data"
