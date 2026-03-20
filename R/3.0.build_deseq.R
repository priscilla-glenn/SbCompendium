#' Build a DESeq2 object from a compendium expression table
#'
#' Takes an expression table containing columns ending in `_TPM` and `_counts`,
#' filters genes by a TPM threshold (max TPM across samples), and constructs a
#' DESeq2::DESeqDataSet with a `~ condition` design where condition is
#' inferred from sample base names.
#'
#' The input table should contain genes as rows and expression columns ending
#' in `_TPM` and `_counts`. Gene identifiers should be stored as row names.
#'
#' @param df A data.frame/tibble containing expression columns.
#' @param TPM_thresh Numeric. Genes are kept if max TPM across `_TPM`
#'   columns is > TPM_thresh. Default = 5.
#' @param ref_level Optional character. If supplied, sets the reference
#'   level for the condition factor.
#' @param design_formula A formula for the DESeq2 design.
#'   Default `~ condition`.
#' @param run_deseq Logical. If TRUE (default), runs DESeq2::DESeq()
#'   before returning.
#'
#' @return A list containing:
#' \itemize{
#'   \item dds: DESeq2::DESeqDataSet
#'   \item subset_df: filtered expression table
#'   \item tpm_cols: TPM column names used
#'   \item count_cols: count column names used
#'   \item sample_names: inferred group names (aligned to count columns)
#'   \item colData: DataFrame used for DESeq2
#'   \item sample_table: data.frame mapping sample columns to groups
#'   \item groups: character vector of group names
#'   \item rep_counts: table of replicates per group
#'   \item TPM_thresh: threshold used
#' }
#' @export
#'
#' @importFrom matrixStats rowMaxs
#' @importFrom S4Vectors DataFrame
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#'
#'
#' @examples
#' \donttest{
#' example_df <- data.frame(
#'   INT1.1_TPM = c(10, 2, 8, 0),
#'   INT1.2_TPM = c(12, 1, 7, 0),
#'   INT2.1_TPM = c(3, 9, 4, 0),
#'   INT2.2_TPM = c(4, 8, 5, 0),
#'   INT1.1_counts = c(100, 20, 80, 0),
#'   INT1.2_counts = c(120, 15, 75, 0),
#'   INT2.1_counts = c(30, 90, 40, 0),
#'   INT2.2_counts = c(40, 85, 50, 0)
#' )
#' rownames(example_df) <- c(
#'   "Sobic.001G000100",
#'   "Sobic.001G000200",
#'   "Sobic.001G000300",
#'   "Sobic.001G000400"
#' )
#'
#' res <- build_deseq(example_df, TPM_thresh = 5, run_deseq = FALSE)
#' names(res)
#' res$sample_table
#' }
build_deseq <- function(df,
                        TPM_thresh = 5,
                        ref_level = NULL,
                        design_formula = ~ condition,
                        run_deseq = TRUE) {

    stopifnot(is.data.frame(df))

    if (!is.numeric(TPM_thresh) || length(TPM_thresh) != 1) {
        stop("TPM_thresh must be a single numeric value.", call. = FALSE)
    }

    # Keep only TPM + count columns
    expr_cols <- grep("(_TPM$|_counts$)", names(df), value = TRUE)
    if (length(expr_cols) == 0) {
        stop("No columns ending with '_TPM' or '_counts' found.", call. = FALSE)
    }

    df2 <- df[, expr_cols, drop = FALSE]

    # Identify TPM and count columns
    tpm_cols <- grep("_TPM$", names(df2), value = TRUE)
    count_cols <- grep("_counts$", names(df2), value = TRUE)

    if (length(tpm_cols) == 0) stop("No '_TPM' columns found.", call. = FALSE)
    if (length(count_cols) == 0) stop("No '_counts' columns found.", call. = FALSE)

    # Filter genes by max TPM
    tpm_mat <- as.matrix(df2[, tpm_cols, drop = FALSE])
    keep <- matrixStats::rowMaxs(tpm_mat, useNames = FALSE) > TPM_thresh
    subset_df <- df2[keep, , drop = FALSE]

    # Count matrix
    count_cols <- grep("_counts$", names(subset_df), value = TRUE)
    countData <- subset_df[, count_cols, drop = FALSE]

    # Infer group names from count columns (aligned to countData column order)
    sample_names <- sub("([._-][0-9]+|[._-]?R[0-9]+)_counts$", "", count_cols)
    condition <- factor(sample_names, levels = unique(sample_names))

    # Relevel if requested
    if (!is.null(ref_level)) {
        if (!ref_level %in% levels(condition)) {
            stop(paste0(
                "ref_level '", ref_level, "' not found. Available levels: ",
                paste(levels(condition), collapse = ", ")
            ), call. = FALSE)
        }
        condition <- stats::relevel(condition, ref = ref_level)
    }

    colData <- DataFrame(condition = condition)

    # Ensure integer count matrix
    countData_mat <- as.matrix(countData)
    mode(countData_mat) <- "numeric"
    storage.mode(countData_mat) <- "integer"

    # Basic sanity checks (optional, but helpful)
    if (any(is.na(countData_mat))) stop("countData contains NA values.", call. = FALSE)
    if (any(countData_mat < 0)) stop("countData contains negative values.", call. = FALSE)

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = countData_mat,
        colData = colData,
        design = design_formula
    )

    if (isTRUE(run_deseq)) {
        dds <- DESeq2::DESeq(dds)
    }

    sample_table <- data.frame(
        sample = count_cols,
        group  = as.character(condition),
        stringsAsFactors = FALSE
    )

    list(
        dds = dds,
        subset_df = subset_df,
        tpm_cols = tpm_cols,
        count_cols = count_cols,
        sample_names = sample_names,
        colData = colData,
        sample_table = sample_table,
        groups = levels(dds$condition),
        rep_counts = table(dds$condition),
        TPM_thresh = TPM_thresh
    )
}
