#' Run a DESeq2 contrast and return results + filtered top DE table
#'
#' @param dds A DESeq2::DESeqDataSet with a `condition` column in colData.
#' @param groupA,groupB Condition levels to compare (groupA vs groupB).
#' @param alpha Adjusted p-value cutoff.
#' @param lfc_thresh Absolute log2FC cutoff.
#' column in `colData(dds)` used for contrasts.
#'
#' @return A list with:
#'   - res: DESeqResults
#'   - res_df: data.frame with GeneIDV3 column
#'   - top_de: filtered data.frame
#' @export
#'
#'
#' @examples
#' # Create a small count matrix
#' counts <- matrix(
#'   c(100,120,30,40,
#'     20,15,90,85,
#'     80,75,40,50),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' rownames(counts) <- c(
#'   "Sobic.001G000100",
#'   "Sobic.001G000200",
#'   "Sobic.001G000300"
#' )
#'
#' colnames(counts) <- c(
#'   "INT1_rep1",
#'   "INT1_rep2",
#'   "INT2_rep1",
#'   "INT2_rep2"
#' )
#'
#' # Sample metadata
#' colData <- S4Vectors::DataFrame(
#'   condition = factor(c("INT1","INT1","INT2","INT2"))
#' )
#'
#' # Build DESeqDataSet
#' dds <- DESeq2::DESeqDataSetFromMatrix(
#'   countData = counts,
#'   colData = colData,
#'   design = ~ condition
#' )
#'
#' # Run DESeq quickly for the example
#' dds <- DESeq2::DESeq(dds, quiet = TRUE)
#'
#' # Run contrast (INT2 vs INT1) with alpha = 0.05
#' # Only genes with adjusted p-value < 0.05 and |log2FC| > 2.32 are returned in top_de
#' res <- run_contrast(dds, "INT2", "INT1", alpha = 0.05, lfc_thresh = 2.32)
#'
#' # Inspect filtered DE genes
#' res$top_de
#'
run_contrast <- function(dds,
                         groupA,
                         groupB,
                         alpha = 0.05,
                         lfc_thresh = 2.32) {

    if (!inherits(dds, "DESeqDataSet")) {
        stop("dds must be a DESeqDataSet.", call. = FALSE)
    }

    # Optional friendly note if DESeq() hasn't been run yet
    if (is.null(SummarizedExperiment::mcols(dds)$dispersion)) {
        message("Note: if DESeq2::results() errors, run DESeq2::DESeq(dds) first (or set run_deseq=TRUE in build_deseq()).")
    }

    groups <- levels(dds$condition)
    if (!groupA %in% groups) {
        stop(sprintf("groupA '%s' not found. Available: %s",
                     groupA, paste(groups, collapse = ", ")), call. = FALSE)
    }
    if (!groupB %in% groups) {
        stop(sprintf("groupB '%s' not found. Available: %s",
                     groupB, paste(groups, collapse = ", ")), call. = FALSE)
    }

    res <- DESeq2::results(dds, contrast = c("condition", groupA, groupB))
    res_df <- as.data.frame(res)
    res_df <- tibble::rownames_to_column(res_df, var = "GeneIDV3")

    top_de <- res_df |>
        dplyr::filter(!is.na(.data$padj)) |>
        dplyr::filter(.data$padj < alpha, abs(.data$log2FoldChange) > lfc_thresh)

    list(
        res = res,
        res_df = res_df,
        top_de = top_de
    )
}


