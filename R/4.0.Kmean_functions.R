#' Build a TPMmean matrix for k-means from top_DE gene list + compendium TPMmeans
#'
#' Subsets TPMmean expression (from `gene_expr`, with GeneIDV3 as rownames) to the
#' GeneIDV3s in `top_de`, then applies adaptive preprocessing:
#' - If >2 TPMmean columns: optional row scaling (typical for pattern clustering)
#' - If <=2 TPMmean columns: avoid row scaling (prevents ±0.707 collapse); optionally log2(x+1)
#'
#' @param top_de A data.frame containing a gene ID column (default "GeneIDV3").
#' @param gene_expr A data.frame/matrix with GeneIDV3s as rownames and TPMmean columns.
#' @param GeneIDV3_col Character. Name of the gene ID column in `top_de`. Default "GeneIDV3".
#' @param tpmmean_cols Optional character vector of TPMmean column names to use.
#'   If NULL, columns are detected via `.get_tpmmean_cols(gene_expr)`.
#' @param keep_order Logical. If TRUE (default), keep the gene order as it appears in top_de.
#' @param drop_missing Logical. If TRUE (default), drop genes found in top_de but not in gene_expr.
#'   If FALSE, error if any top_de genes are missing from gene_expr.
#' @param na_action What to do with NA in the TPMmean matrix: "omit" drops rows with any NA,
#'   "fail" errors. Default "omit".
#' @param scale_rows_if_gt2 Logical. If TRUE (default), scale each gene (row) across samples
#'   when there are >2 TPMmean columns.
#' @param log_transform_if_le2 Logical. If TRUE (default), apply log2(x+1) when there are <=2 columns.
#' @param scale_cols Logical. If TRUE, scale each sample (column) across genes (usually FALSE).
#' @param drop_zero_var Logical. If TRUE (default), drop genes with zero variance across TPMmean columns.
#' @param groups Description of groups
#'
#' @return A numeric matrix: rows = genes, cols = TPMmean samples/conditions.
#'
#' @export
#'
#' @examples
#' # top_de table with GeneIDV3 column
#' top_de <- data.frame(
#'   GeneIDV3 = c("Sobic.001G000100", "Sobic.001G000200", "Sobic.001G999999")
#' )
#'
#' # gene_expr must have GeneIDV3s as rownames + TPMmean columns
#' gene_expr <- data.frame(
#'   A_TPMmean = c(10, 0, 5),
#'   B_TPMmean = c(12, 1, 4),
#'   C_TPMmean = c(3,  8, 6)
#' )
#' rownames(gene_expr) <- c("Sobic.001G000100", "Sobic.001G000200", "Sobic.001G000300")
#'
#' # Build matrix (drops missing gene Sobic.001G999999 by default)
#' mat <- prepare_kmeans_matrix_topDE_TPMmean(
#'   top_de,
#'   gene_expr,
#'   tpmmean_cols = c("A_TPMmean", "B_TPMmean", "C_TPMmean"),
#'   scale_rows_if_gt2 = TRUE
#' )
#'
#' dim(mat)
#' mat[1:2, , drop = FALSE]
#'
#' @examples
#' # Two-column case: avoids row scaling, optionally log2(x+1)
#' gene_expr2 <- gene_expr[, c("A_TPMmean", "B_TPMmean"), drop = FALSE]
#' mat2 <- prepare_kmeans_matrix_topDE_TPMmean(
#'   top_de,
#'   gene_expr2,
#'   tpmmean_cols = c("A_TPMmean", "B_TPMmean"),
#'   log_transform_if_le2 = TRUE
#' )
#' mat2
prepare_kmeans_matrix_topDE_TPMmean <- function(top_de,
                                                gene_expr,
                                                GeneIDV3_col = "GeneIDV3",
                                                tpmmean_cols = NULL,
                                                groups = NULL,
                                                keep_order = TRUE,
                                                drop_missing = TRUE,
                                                na_action = c("omit", "fail"),
                                                scale_rows_if_gt2 = TRUE,
                                                log_transform_if_le2 = TRUE,
                                                scale_cols = FALSE,
                                                drop_zero_var = TRUE) {

    stopifnot(is.data.frame(top_de))

    if (!is.data.frame(gene_expr) && !is.matrix(gene_expr)) {
        stop("gene_expr must be a data.frame or matrix.", call. = FALSE)
    }

    na_action <- match.arg(na_action)

    if (!GeneIDV3_col %in% names(top_de)) {
        stop("top_de is missing GeneIDV3_col = '", GeneIDV3_col, "'.", call. = FALSE)
    }

    gene_expr_ids <- rownames(gene_expr)

    if (is.null(gene_expr_ids) || anyNA(gene_expr_ids) || any(gene_expr_ids == "")) {
        stop("gene_expr must have non-empty rownames containing GeneIDV3s.", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Determine TPMmean columns
    # ------------------------------------------------------------

    if (is.null(tpmmean_cols)) {
        tpmmean_cols <- .get_tpmmean_cols(gene_expr)
    } else {
        missing_cols <- setdiff(tpmmean_cols, colnames(gene_expr))
        if (length(missing_cols) > 0) {
            stop("These tpmmean_cols are not in gene_expr: ",
                 paste(missing_cols, collapse = ", "), call. = FALSE)
        }
    }

    # Restrict to selected groups if requested
    if (!is.null(groups)) {

        if (!is.character(groups) || length(groups) != 2) {
            stop(
                "groups must be a character vector of length 2 (e.g. c('INT1','INT2')).",
                call. = FALSE
            )
        }

        keep <- grepl(groups[1], tpmmean_cols, fixed = TRUE) |
            grepl(groups[2], tpmmean_cols, fixed = TRUE)

        tpmmean_cols <- tpmmean_cols[keep]

        if (length(tpmmean_cols) == 0) {
            stop(
                "No TPMmean columns matched the specified groups. ",
                "Looked for columns containing: ",
                paste(groups, collapse = ", "),
                call. = FALSE
            )
        }
    }

    if (length(tpmmean_cols) < 2) {
        stop("Need at least two TPMmean columns for clustering.", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Determine genes to include
    # ------------------------------------------------------------

    genes <- as.character(top_de[[GeneIDV3_col]])
    genes <- genes[!is.na(genes) & genes != ""]

    if (length(genes) == 0) {
        stop("No valid Gene IDs found in top_de.", call. = FALSE)
    }

    genes <- genes[!duplicated(genes)]

    # Match genes to gene_expr rownames
    idx <- match(genes, gene_expr_ids)

    missing_genes <- genes[is.na(idx)]

    if (length(missing_genes) > 0) {

        if (!drop_missing) {
            stop(
                "Some Gene IDs from top_de are missing in gene_expr rownames (showing up to 20): ",
                paste(head(missing_genes, 20), collapse = ", "),
                call. = FALSE
            )
        }

        keep <- !is.na(idx)
        genes <- genes[keep]
        idx <- idx[keep]
    }

    expr_df <- gene_expr[idx, tpmmean_cols, drop = FALSE]

    # ------------------------------------------------------------
    # Convert to numeric matrix
    # ------------------------------------------------------------

    expr_mat <- as.matrix(expr_df)

    suppressWarnings(storage.mode(expr_mat) <- "numeric")

    # ------------------------------------------------------------
    # Handle NA values
    # ------------------------------------------------------------

    if (anyNA(expr_mat)) {

        if (na_action == "fail") {
            stop(
                "TPMmean matrix contains NA values. ",
                "Set na_action='omit' to drop rows.",
                call. = FALSE
            )

        } else {

            keep <- stats::complete.cases(expr_mat)

            expr_mat <- expr_mat[keep, , drop = FALSE]
            genes <- genes[keep]
        }
    }

    rownames(expr_mat) <- make.unique(genes)

    # ------------------------------------------------------------
    # Drop zero-variance rows (helps kmeans)
    # ------------------------------------------------------------

    if (drop_zero_var && nrow(expr_mat) > 0) {

        vars <- apply(expr_mat, 1, stats::var, na.rm = TRUE)

        keep <- is.finite(vars) & vars > 0

        expr_mat <- expr_mat[keep, , drop = FALSE]
    }

    # ------------------------------------------------------------
    # Adaptive preprocessing
    # ------------------------------------------------------------

    p <- ncol(expr_mat)

    if (p <= 2) {

        # Avoid row scaling collapse; optional log transform

        if (log_transform_if_le2) {

            if (any(expr_mat < 0, na.rm = TRUE)) {

                warning(
                    "expr_mat contains negative values; skipping log2(x+1) transform.",
                    call. = FALSE
                )

            } else {

                expr_mat <- log2(expr_mat + 1)
            }
        }

    } else {

        # Typical gene-pattern clustering: scale each gene

        if (scale_rows_if_gt2) {

            expr_mat <- t(scale(t(expr_mat)))
        }
    }

    # Optional column scaling

    if (scale_cols) {

        expr_mat <- scale(expr_mat)
    }

    # ------------------------------------------------------------
    # Optional row ordering
    # ------------------------------------------------------------

    if (!keep_order && nrow(expr_mat) > 1) {

        ord <- order(rownames(expr_mat))

        expr_mat <- expr_mat[ord, , drop = FALSE]
    }

    return(expr_mat)
}



#' Evaluate optimal k for k-means clustering (prints plots)
#'
#' Generates and prints diagnostic plots to help determine the optimal number
#' of clusters (k) using the silhouette and WSS methods.
#'
#' @param mat Numeric matrix with genes as rows and samples as columns.
#' @param k_max Integer. Maximum number of clusters to evaluate. Default 20.
#' @param dedup Logical. If TRUE (default), remove duplicate rows for k-evaluation only.
#' @param nstart Integer. Passed to kmeans during evaluation. Default 50.
#' @param iter_max Integer. Passed to kmeans during evaluation. Default 1000.
#' @param seed Optional integer seed for reproducibility.
#' @param algorithm Character. kmeans algorithm ("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen").
#'   Default "Lloyd" for stability.
#' @param print_plots Logical. If TRUE (default), prints the plots immediately.
#'
#' @return Invisibly returns a list with silhouette_plot, wss_plot, and meta.
#' @export
#'
#' @examples
#' # Small numeric matrix (genes x samples)
#' mat <- matrix(
#'   c(1,2,3,
#'     2,1,3,
#'     9,8,7,
#'     8,9,7,
#'     4,4,4),
#'   nrow = 5,
#'   byrow = TRUE
#' )
#' rownames(mat) <- paste0("gene", 1:5)
#' colnames(mat) <- c("A_TPMmean", "B_TPMmean", "C_TPMmean")
#'
#' if (requireNamespace("factoextra", quietly = TRUE)) {
#'   out <- evaluate_kmeans_k(mat, k_max = 4, seed = 1, print_plots = FALSE)
#'   names(out)
#' } else {
#'   message("Skipping: 'factoextra' not installed.")
#' }
evaluate_kmeans_k <- function(mat,
                              k_max = 20,
                              dedup = TRUE,
                              nstart = 50,
                              iter_max = 1000,
                              seed = NULL,
                              algorithm = "Lloyd",
                              print_plots = TRUE) {

    if (!is.matrix(mat)) mat <- as.matrix(mat)
    if (!is.numeric(mat)) stop("mat must be a numeric matrix.", call. = FALSE)
    if (nrow(mat) < 3) stop("Need at least 3 rows for k evaluation.", call. = FALSE)

    if (!requireNamespace("factoextra", quietly = TRUE)) {
        stop("Package 'factoextra' is required for evaluate_kmeans_k().", call. = FALSE)
    }

    mat_eval <- mat

    # Optionally deduplicate rows (common after scaling / discretization)
    if (dedup) {
        keep <- !duplicated(mat_eval)
        mat_eval <- mat_eval[keep, , drop = FALSE]
    }

    distinct_n <- nrow(mat_eval)

    # kmeans requires centers < nrow
    k_cap <- min(as.integer(k_max), distinct_n - 1)
    if (k_cap < 2) {
        stop(
            "Not enough distinct rows to evaluate k. Distinct rows = ", distinct_n,
            ". Try turning off scaling, lowering filtering, or evaluate on more genes.",
            call. = FALSE
        )
    }

    # Custom kmeans wrapper with safer defaults for evaluation
    kmeans_fun <- function(x, centers) {
        if (!is.null(seed)) set.seed(seed)
        stats::kmeans(
            x,
            centers = centers,
            nstart = nstart,
            iter.max = iter_max,
            algorithm = algorithm
        )
    }

    silhouette_plot <- factoextra::fviz_nbclust(
        mat_eval,
        FUNcluster = kmeans_fun,
        method = "silhouette",
        k.max = k_cap
    )

    wss_plot <- factoextra::fviz_nbclust(
        mat_eval,
        FUNcluster = kmeans_fun,
        method = "wss",
        k.max = k_cap
    )

    if (isTRUE(print_plots)) {
        print(silhouette_plot)
        print(wss_plot)
    }

    out <- list(
        silhouette_plot = silhouette_plot,
        wss_plot = wss_plot,
        meta = list(
            input_rows = nrow(mat),
            eval_rows = nrow(mat_eval),
            dedup = dedup,
            k_max_requested = k_max,
            k_max_used = k_cap,
            nstart = nstart,
            iter_max = iter_max,
            algorithm = algorithm
        )
    )

    return(invisible(out))
}






#' Run k-means clustering on a prepared expression matrix
#'
#' Runs stats::kmeans on a numeric matrix (genes x samples), returning the kmeans
#' object and a tidy cluster assignment table.
#'
#' @param mat Numeric matrix with genes as rows and samples as columns.
#'   Typically the output of `prepare_kmeans_matrix_topDE_TPMmean()`.
#' @param centers Integer. Number of clusters (k).
#' @param nstart Integer. Number of random starts for kmeans. Default 25.
#' @param iter_max Integer. Maximum number of iterations. Default 100.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item kmeans: the stats::kmeans result
#'   \item clusters: data.frame with GeneIDV3 and cluster
#' }
#'
#' @export
#'
#'
#' @examples
#' # Example matrix (genes x samples)
#' mat <- matrix(
#'   c(0, 1, 2,
#'     0, 2, 3,
#'     8, 7, 6,
#'     9, 8, 7),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#' rownames(mat) <- c("Sobic.001G000100","Sobic.001G000200",
#'                    "Sobic.001G000300","Sobic.001G000400")
#' colnames(mat) <- c("A_TPMmean","B_TPMmean","C_TPMmean")
#'
#' km_out <- run_kmeans(mat, centers = 2, seed = 1)
#' km_out$clusters
#' km_out$kmeans$centers
#'
run_kmeans <- function(mat,
                       centers,
                       nstart = 25,
                       iter_max = 100,
                       seed = NULL) {

    if (!is.matrix(mat)) mat <- as.matrix(mat)
    if (!is.numeric(mat)) stop("mat must be a numeric matrix.", call. = FALSE)
    if (nrow(mat) < 2) stop("mat must have at least 2 rows (genes).", call. = FALSE)
    if (ncol(mat) < 2) stop("mat must have at least 2 columns (samples).", call. = FALSE)

    if (!is.numeric(centers) || length(centers) != 1 || centers < 1) {
        stop("centers must be a single integer >= 1.", call. = FALSE)
    }
    centers <- as.integer(centers)

    if (!is.null(seed)) set.seed(seed)

    km <- stats::kmeans(mat, centers = centers, nstart = nstart, iter.max = iter_max)

    clusters <- data.frame(
        GeneIDV3 = rownames(mat),
        cluster = km$cluster,
        stringsAsFactors = FALSE
    )

    return(list(kmeans = km, clusters = clusters))
}
