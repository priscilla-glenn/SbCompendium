#' Plot PCA and heatmaps for a DESeq2 dataset
#'
#' Applies DESeq2's variance-stabilizing transformation (VST), creates a
#' labeled and points-only PCA plots, calculates a sample-to-sample correlation
#' matrix, and draws correlation and scaled-expression heatmaps. This function is intended
#' for quality control immediately after `build_deseq()`.
#'
#' @param x Output from `build_deseq()` or a `DESeqDataSet`.
#' @param intgroup Character. Single `colData` column used to color the PCA and
#'   annotate the expression heatmap. Default `"condition"`.
#' @param blind Logical. Passed to the DESeq2 variance-stabilizing
#'   transformation. Default `TRUE`.
#' @param cor_method Character. Correlation method passed to `stats::cor()`.
#'   One of `"pearson"`, `"kendall"`, or `"spearman"`.
#' @param heatmap_scale Character. Scaling applied to the expression heatmap.
#'   One of `"row"`, `"column"`, or `"none"`. Default `"row"`.
#' @param show_rownames Logical. Show gene names on the expression heatmap.
#'   Default `FALSE`.
#' @param cor_cellheight Numeric. Height of each row in the sample-correlation
#'   heatmap. Increase this value when sample labels overlap. Default 20.
#' @param pca_point_size Numeric. PCA point size. Default 4.
#' @param pca_text_vjust Numeric. Vertical adjustment for PCA sample labels.
#'   Default -0.8.
#' @param print_plots Logical. If `TRUE` (default), prints the PCA and both
#'   heatmaps. Set to `FALSE` to create and return them without printing.
#'
#' @return A list containing:
#' \itemize{
#'   \item `vsd`: the variance-stabilized `DESeqTransform` object
#'   \item `pca_plot`: labeled ggplot2 PCA plot
#'   \item `pca_plot_points`: ggplot2 PCA plot with points only
#'   \item `pca_data`: data used for the PCA plot
#'   \item `percent_var`: percent variance explained by PC1 and PC2
#'   \item `cor_matrix`: sample-to-sample correlation matrix
#'   \item `correlation_heatmap`: pheatmap object for `cor_matrix`
#'   \item `expression_heatmap`: pheatmap object for the scaled VST assay
#'   \item `annotation_col`: sample annotation used by the expression heatmap
#' }
#'
#' @export
#' @examples
#' \donttest{
#' data(example_data)
#' gene_expr <- importTable(example_data[["nodal_buds_28"]])
#' dds_out <- build_deseq(gene_expr, TPM_thresh = 5, run_deseq = TRUE)
#'
#' qc <- plot_deseq_qc(dds_out, print_plots = FALSE)
#' qc$pca_plot
#' qc$pca_plot_points
#' qc$cor_matrix
#' }
plot_deseq_qc <- function(x,
                          intgroup = "condition",
                          blind = TRUE,
                          cor_method = c("pearson", "kendall", "spearman"),
                          heatmap_scale = c("row", "column", "none"),
                          show_rownames = FALSE,
                          cor_cellheight = 20,
                          pca_point_size = 4,
                          pca_text_vjust = -0.8,
                          print_plots = TRUE) {

    cor_method <- match.arg(cor_method)
    heatmap_scale <- match.arg(heatmap_scale)

    dds <- if (inherits(x, "DESeqDataSet")) {
        x
    } else if (is.list(x) && inherits(x$dds, "DESeqDataSet")) {
        x$dds
    } else {
        stop("x must be output from build_deseq() or a DESeqDataSet.",
             call. = FALSE)
    }

    if (!is.character(intgroup) || length(intgroup) != 1L) {
        stop("intgroup must be a single colData column name.", call. = FALSE)
    }

    if (!intgroup %in% names(SummarizedExperiment::colData(dds))) {
        stop("intgroup '", intgroup, "' was not found in colData(dds).",
             call. = FALSE)
    }

    if (!is.numeric(cor_cellheight) || length(cor_cellheight) != 1L ||
        !is.finite(cor_cellheight) || cor_cellheight <= 0) {
        stop("cor_cellheight must be a single positive number.", call. = FALSE)
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.", call. = FALSE)
    }

    vsd <- tryCatch(
        DESeq2::vst(dds, blind = blind),
        error = function(e) {
            if (!grepl("less than.*nsub", conditionMessage(e))) {
                stop(e)
            }

            message(
                "Using DESeq2::varianceStabilizingTransformation() because ",
                "the dataset has too few eligible genes for the vst() shortcut."
            )
            DESeq2::varianceStabilizingTransformation(dds, blind = blind)
        }
    )
    pca_data <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))

    pca_base <- ggplot2::ggplot(
        pca_data,
        ggplot2::aes(
            x = .data$PC1,
            y = .data$PC2,
            color = .data[[intgroup]]
        )
    ) +
        ggplot2::xlab(paste0("PC1: ", percent_var[1], "% variance")) +
        ggplot2::ylab(paste0("PC2: ", percent_var[2], "% variance")) +
        ggplot2::theme_classic()

    pca_plot_points <- pca_base +
        ggplot2::geom_point(size = pca_point_size)

    pca_plot <- pca_plot_points +
        ggplot2::geom_text(
            ggplot2::aes(label = .data$name),
            vjust = pca_text_vjust
        )

    vsd_assay <- SummarizedExperiment::assay(vsd)
    cor_matrix <- stats::cor(vsd_assay, method = cor_method)

    annotation_col <- as.data.frame(
        SummarizedExperiment::colData(dds)[, intgroup, drop = FALSE]
    )

    if (isTRUE(print_plots)) {
        print(pca_plot)
        print(pca_plot_points)
    }

    correlation_heatmap <- pheatmap::pheatmap(
        cor_matrix,
        cellheight = cor_cellheight,
        silent = !isTRUE(print_plots)
    )

    expression_heatmap <- pheatmap::pheatmap(
        vsd_assay,
        scale = heatmap_scale,
        annotation_col = annotation_col,
        show_rownames = show_rownames,
        silent = !isTRUE(print_plots)
    )

    list(
        vsd = vsd,
        pca_plot = pca_plot,
        pca_plot_points = pca_plot_points,
        pca_data = pca_data,
        percent_var = percent_var,
        cor_matrix = cor_matrix,
        correlation_heatmap = correlation_heatmap,
        expression_heatmap = expression_heatmap,
        annotation_col = annotation_col
    )
}
