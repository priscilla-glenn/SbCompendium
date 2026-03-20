#' Export k-means clusters as DAVID-ready gene lists
#'
#' Exports one text file per cluster containing gene IDs converted from
#' `Sobic.*` format to DAVID-compatible `SORBI_3*` format (v3 naming),
#' suitable for GO/annotation tools like DAVID.
#'
#' Files are written as one gene ID per line.
#'
#' @param df A data.frame containing at least `GeneIDV3` and `cluster` columns
#'   (e.g., output from `build_heatmap_cluster_export()`).
#' @param clusters Integer vector of cluster numbers to export. If NULL (default),
#'   exports all unique clusters found in `df$cluster`.
#' @param out_dir Output directory for the text files.
#' @param gene_col Name of the Gene ID column. Default "GeneIDV3".
#' @param cluster_col Name of the cluster column. Default "cluster".
#' @param prefix Filename prefix. Default "KMeans_cluster".
#' @param quiet Logical. If FALSE (default), prints files written.
#'
#' @return Invisibly returns a character vector of file paths written.
#'
#' @export
#'
#' @examples
#' # Example cluster table (similar to build_heatmap() output)
#' df <- data.frame(
#'   GeneIDV3 = c(
#'     "Sobic.001G000100",
#'     "Sobic.001G000200",
#'     "Sobic.001G000300",
#'     "Sobic.001G000400"
#'   ),
#'   cluster = c(1, 1, 2, 3),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Write DAVID-ready gene lists to a temporary directory
#' out_dir <- tempdir()
#'
#' files_written <- export_clusters_david(
#'   df,
#'   out_dir = out_dir
#' )
#'
#' files_written
#'
#' # Export only selected clusters
#' files_written2 <- export_clusters_david(
#'   df,
#'   clusters = c(1, 2),
#'   out_dir = out_dir,
#'   quiet = TRUE
#' )
#'
#' files_written2
#'
export_clusters_david <- function(df,
                                  clusters = NULL,
                                  out_dir = ".",
                                  gene_col = "GeneIDV3",
                                  cluster_col = "cluster",
                                  prefix = "KMeans_cluster",
                                  quiet = FALSE) {

    stopifnot(is.data.frame(df))
    if (!gene_col %in% names(df)) stop("df is missing gene_col = '", gene_col, "'.", call. = FALSE)
    if (!cluster_col %in% names(df)) stop("df is missing cluster_col = '", cluster_col, "'.", call. = FALSE)

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Determine clusters to export
    if (is.null(clusters)) {
        clusters <- sort(unique(df[[cluster_col]]))
        clusters <- clusters[!is.na(clusters)]
    }

    # Convert gene IDs to DAVID format
    convert_to_david <- function(x) {
        x <- as.character(x)
        x <- sub("^Sobic\\.(\\d+)G(\\d+)$", "SORBI_3\\1G\\2", x)
        x <- sub("^Sobic\\.(K\\d+)$", "SORBI_3\\1", x)
        x
    }

    written <- character(0)

    for (k in clusters) {

        sub_df <- df[df[[cluster_col]] == k, , drop = FALSE]

        genes <- convert_to_david(sub_df[[gene_col]])
        genes <- genes[!duplicated(genes)]

        filename <- paste0(prefix, k, ".txt")
        output_file <- file.path(out_dir, filename)

        writeLines(genes, con = output_file)

        written <- c(written, output_file)

        if (!quiet) {
            message("Wrote: ", output_file, " (", length(genes), " genes)")
        }
    }

    invisible(written)
}
