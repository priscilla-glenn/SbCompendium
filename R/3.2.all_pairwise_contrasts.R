#' Run all pairwise DESeq2 contrasts for a condition factor
#'
#' Given a DESeq2::DESeqDataSet, runs DESeq2::results() for every pairwise
#' comparison among the levels of a design factor (default: "condition").
#'
#' @param dds A DESeq2::DESeqDataSet. If not already run through DESeq(), set run_deseq = TRUE.
#' @param factor_name Character. The column in colData(dds) used for contrasts. Default "condition".
#' @param run_deseq Logical. If TRUE, runs DESeq2::DESeq(dds) before contrasts. Default FALSE.
#' @param alpha Numeric. Adjusted p-value cutoff used by DESeq2 results(). Default 0.05.
#' @param fc_thresh Numeric. Absolute fold-change cutoff on the **fold-change (FC) scale**, not log2.
#'   For example, `fc_thresh = 5` will filter for genes with at least 5-fold up- or down-regulation.
#'   Internally, this is applied as abs(log2FoldChange) > log2(fc_thresh), so both strong up- and down-regulated genes are retained.
#' @param filter_significant Logical. If TRUE (default), filters to padj < alpha and
#'   abs(log2FoldChange) > log2(fc_thresh) for the combined table and top-gene summary.
#' @param add_fc Logical. If TRUE (default), adds FC = 2^(log2FoldChange) to outputs.
#' @param include_all_rows Logical. If TRUE, returns full unfiltered per-contrast tables in `results_list`.
#'   Default TRUE.
#'
#' @return A list with:
#' \itemize{
#'   \item results_list: named list of data.frames, one per contrast (A_vs_B)
#'   \item all_results: combined long data.frame (optionally filtered)
#'   \item top_by_gene: per-gene summary selecting the contrast with max absolute change
#'   \item contrasts: data.frame of all tested pairs
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Small synthetic dataset: 3 conditions, 2 replicates each
#' set.seed(1)
#'
#' counts <- matrix(
#'   c(
#'     200, 210,  90,  95,  20,  22,
#'      10,  12, 180, 175,  70,  75,
#'      60,  65,  15,  14, 160, 155,
#'      50,  52,  49,  48,  51,  50
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#'
#' rownames(counts) <- c(
#'   "Sobic.001G000100",
#'   "Sobic.001G000200",
#'   "Sobic.001G000300",
#'   "Sobic.001G000400"
#' )
#'
#' colnames(counts) <- c("A_1","A_2","B_1","B_2","C_1","C_2")
#'
#' colData <- S4Vectors::DataFrame(
#'   condition = factor(c("A","A","B","B","C","C"))
#' )
#'
#' dds <- DESeq2::DESeqDataSetFromMatrix(
#'   countData = counts,
#'   colData = colData,
#'   design = ~ condition
#' )
#'
#' # Run DESeq once, then compute all pairwise contrasts
#' dds <- DESeq2::DESeq(dds, quiet = TRUE)
#'
#' out <- run_all_pairwise_contrasts(
#'   dds,
#'   factor_name = "condition",
#'   run_deseq = FALSE,
#'   alpha = 0.1,
#'   fc_thresh = 5,
#'   filter_significant = TRUE,
#'   add_fc = TRUE,
#'   include_all_rows = TRUE
#' )
#'
#' # What contrasts were tested?
#' out$contrasts
#'
#' # Combined filtered results across all contrasts
#' head(out$all_results)
#'
#' # Per-gene “best” contrast
#' out$top_by_gene
#' }
#'
run_all_pairwise_contrasts <- function(dds,
                                       factor_name = "condition",
                                       run_deseq = FALSE,
                                       alpha = 0.05,
                                       fc_thresh = 5,
                                       filter_significant = TRUE,
                                       add_fc = TRUE,
                                       include_all_rows = TRUE) {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.", call. = FALSE)
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Package 'dplyr' is required.", call. = FALSE)
    }

    if (!factor_name %in% names(SummarizedExperiment::colData(dds))) {
        stop("factor_name '", factor_name, "' not found in colData(dds).", call. = FALSE)
    }

    SummarizedExperiment::colData(dds)[[factor_name]] <-
        as.factor(SummarizedExperiment::colData(dds)[[factor_name]])

    lvls <- levels(SummarizedExperiment::colData(dds)[[factor_name]])

    if (length(lvls) < 2) {
        stop("Need at least 2 levels in ", factor_name, ".", call. = FALSE)
    }

    if (!is.numeric(fc_thresh) || length(fc_thresh) != 1 || fc_thresh <= 0) {
        stop("fc_thresh must be a single numeric value > 0.", call. = FALSE)
    }

    lfc_cutoff <- log2(fc_thresh)

    if (isTRUE(run_deseq)) {
        dds <- DESeq2::DESeq(dds)
    }

    pairs <- utils::combn(lvls, 2)

    contrasts_df <- data.frame(
        groupA = pairs[1, ],
        groupB = pairs[2, ],
        stringsAsFactors = FALSE
    )

    results_list <- vector("list", nrow(contrasts_df))
    names(results_list) <- paste0(contrasts_df$groupA, "_vs_", contrasts_df$groupB)

    for (i in seq_len(nrow(contrasts_df))) {

        A <- contrasts_df$groupA[i]
        B <- contrasts_df$groupB[i]
        contrast_name <- paste0(A, " vs ", B)

        res <- DESeq2::results(dds, contrast = c(factor_name, A, B), alpha = alpha)
        res_df <- as.data.frame(res)
        res_df$GeneIDV3 <- rownames(res_df)
        res_df$contrast <- contrast_name

        if (isTRUE(add_fc)) {
            res_df$FC <- 2^(res_df$log2FoldChange)
        }

        results_list[[i]] <- res_df
    }

    all_results <- dplyr::bind_rows(results_list)

    if (isTRUE(filter_significant)) {
        if (!"padj" %in% names(all_results)) {
            stop("Expected padj in results.", call. = FALSE)
        }
        if (!"log2FoldChange" %in% names(all_results)) {
            stop("Expected log2FoldChange in results.", call. = FALSE)
        }

        all_results <- all_results |>
            dplyr::filter(!is.na(.data$padj), .data$padj < alpha) |>
            dplyr::filter(!is.na(.data$log2FoldChange),
                          abs(.data$log2FoldChange) > lfc_cutoff)
    }

    if (isTRUE(add_fc)) {
        top_by_gene <- all_results |>
            dplyr::group_by(.data$GeneIDV3) |>
            dplyr::slice_max(order_by = abs(.data$log2FoldChange),
                             n = 1,
                             with_ties = FALSE) |>
            dplyr::ungroup() |>
            dplyr::select(
                .data$GeneIDV3,
                .data$FC,
                .data$padj,
                .data$contrast,
                .data$log2FoldChange,
                dplyr::everything()
            )
    } else {
        top_by_gene <- all_results |>
            dplyr::group_by(.data$GeneIDV3) |>
            dplyr::slice_max(order_by = abs(.data$log2FoldChange),
                             n = 1,
                             with_ties = FALSE) |>
            dplyr::ungroup() |>
            dplyr::select(
                .data$GeneIDV3,
                .data$log2FoldChange,
                .data$padj,
                .data$contrast,
                dplyr::everything()
            )
    }

    if (!isTRUE(include_all_rows) && isTRUE(filter_significant)) {
        results_list <- lapply(results_list, function(df) {
            df <- df[!is.na(df$padj) & df$padj < alpha, , drop = FALSE]
            df <- df[!is.na(df$log2FoldChange) & abs(df$log2FoldChange) > lfc_cutoff, , drop = FALSE]
            df
        })
    }

    list(
        results_list = results_list,
        all_results = all_results,
        top_by_gene = top_by_gene,
        contrasts = contrasts_df
    )
}
