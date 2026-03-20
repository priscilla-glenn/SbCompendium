#' Run UMAP on compendium TPM expression (auto-infers groups + replicates)
#'
#' Builds a gene-by-sample matrix from replicate-level `_TPM` columns (ignores `_TPMmean`),
#' aggregates rows by GeneIDV3, optionally filters by TPM threshold, scales the data,
#' selects the top variable genes by MAD, then runs UMAP on samples.
#'
#' Gene identifiers may be stored either in a `GeneIDV3` column (default) OR as row names.
#'
#' @param df A data.frame containing gene IDs (column or rownames) and replicate-level `_TPM` columns.
#' @param GeneIDV3_col Character. Name of gene ID column if present. Default "GeneIDV3".
#' @param tpm_cols Optional character vector of `_TPM` column names. If NULL, auto-detects.
#' @param aggregate_fun Function used to aggregate to GeneIDV3. Default `sum`.
#' @param TPM_thresh Optional numeric. If not NULL, keeps genes with max TPM > TPM_thresh
#'   prior to scaling/MAD selection.
#' @param scale_data Logical. If TRUE (default), scales (z-scores) each gene across samples.
#' @param n_top_genes Integer. Number of top MAD genes to use. Default 1000.
#' @param seed Optional integer seed for reproducibility (user-picked).
#' @param n_neighbors Integer. UMAP n_neighbors. Default 15.
#' @param min_dist Numeric. UMAP min_dist. Default 0.1.
#' @param metric Character. UMAP distance metric. Default "euclidean".
#' @param make_labels Logical. If TRUE (default), label = \code{"{group}.r{rep}"}.
#' @param return_input Logical. If TRUE, includes the aggregated TPM matrix in output. Default FALSE.
#'
#' @return A list with:
#' \itemize{
#'   \item umap: the object returned by `umap::umap()`
#'   \item umap_df: data.frame with UMAP1, UMAP2, sample, group, label
#'   \item scaled_data_top: gene x sample matrix used for UMAP (after scaling + MAD filter)
#'   \item tpm_cols: TPM columns used
#'   \item sample_table: mapping of sample columns to group/rep
#'   \item groups: factor levels of inferred groups
#'   \item rep_counts: table of replicate counts per group
#'   \item meta: list of parameters used (including seed, n_neighbors, min_dist)
#'   \item input_matrix: (optional) aggregated TPM matrix (genes x samples)
#' }
#'
#' @export
#'
#' @importFrom stats mad
run_umap_tpm <- function(df,
                         GeneIDV3_col = "GeneIDV3",
                         tpm_cols = NULL,
                         aggregate_fun = sum,
                         TPM_thresh = NULL,
                         scale_data = TRUE,
                         n_top_genes = 1000,
                         seed = NULL,
                         n_neighbors = 15,
                         min_dist = 0.1,
                         metric = "euclidean",
                         make_labels = TRUE,
                         return_input = FALSE) {

    stopifnot(is.data.frame(df))

    # Accept GeneIDV3 either as a column or as rownames
    if (!GeneIDV3_col %in% names(df)) {
        rn <- rownames(df)
        if (!is.null(rn) && !anyNA(rn) && all(rn != "")) {
            df[[GeneIDV3_col]] <- rn
        } else {
            stop(
                "df must contain GeneIDV3 either as a column ('", GeneIDV3_col,
                "') or as non-empty rownames.", call. = FALSE
            )
        }
    }

    # Detect TPM columns (explicitly excludes *_TPMmean)
    if (is.null(tpm_cols)) {
        tpm_cols <- grep("_TPM$", names(df), value = TRUE)
        tpm_cols <- tpm_cols[!grepl("_TPMmean$", tpm_cols)]
    } else {
        missing <- setdiff(tpm_cols, names(df))
        if (length(missing) > 0) {
            stop("These tpm_cols are not in df: ", paste(missing, collapse = ", "), call. = FALSE)
        }
    }

    if (length(tpm_cols) < 2) {
        stop("Need at least 2 replicate-level '_TPM' columns for UMAP.", call. = FALSE)
    }

    # Validate UMAP params early
    if (!is.numeric(n_neighbors) || length(n_neighbors) != 1 || n_neighbors < 2) {
        stop("n_neighbors must be a single integer >= 2.", call. = FALSE)
    }
    n_neighbors <- as.integer(n_neighbors)

    if (!is.numeric(min_dist) || length(min_dist) != 1 || min_dist < 0) {
        stop("min_dist must be a single numeric value >= 0.", call. = FALSE)
    }

    # Keep only GeneIDV3 + TPM columns
    data2 <- df[, c(GeneIDV3_col, tpm_cols), drop = FALSE]
    names(data2)[names(data2) == GeneIDV3_col] <- "GeneIDV3"

    # Aggregate to GeneIDV3
    data_ag <- stats::aggregate(. ~ GeneIDV3, data = data2, FUN = aggregate_fun)

    row_names <- data_ag$GeneIDV3
    data_ag <- data_ag[, setdiff(names(data_ag), "GeneIDV3"), drop = FALSE]
    data_ag <- as.data.frame(data_ag)
    rownames(data_ag) <- row_names

    # Infer groups/replicates from TPM column names
    meta_samples <- .infer_groups_from_tpm_cols(colnames(data_ag))

    # Optional TPM threshold filter (max TPM across samples)
    if (!is.null(TPM_thresh)) {
        if (!is.numeric(TPM_thresh) || length(TPM_thresh) != 1) {
            stop("TPM_thresh must be a single numeric value.", call. = FALSE)
        }
        tpm_mat <- as.matrix(data_ag)
        suppressWarnings(storage.mode(tpm_mat) <- "numeric")
        keep <- apply(tpm_mat, 1, max, na.rm = TRUE) > TPM_thresh
        data_ag <- data_ag[keep, , drop = FALSE]
    }

    # Numeric matrix
    mat <- as.matrix(data_ag)
    suppressWarnings(storage.mode(mat) <- "numeric")
    if (anyNA(mat)) stop("TPM matrix contains NA values.", call. = FALSE)

    # Scale genes across samples (z-score per gene)
    if (isTRUE(scale_data)) {
        # scale() scales columns, so transpose to scale genes then transpose back
        mat_scaled <- t(scale(t(mat)))
    } else {
        mat_scaled <- mat
    }

    # MAD gene selection
    if (!is.numeric(n_top_genes) || length(n_top_genes) != 1 || n_top_genes < 1) {
        stop("n_top_genes must be a single integer >= 1.", call. = FALSE)
    }
    n_top <- min(as.integer(n_top_genes), nrow(mat_scaled))
    mad_scores <- apply(mat_scaled, 1, stats::mad)
    top_idx <- order(mad_scores, decreasing = TRUE)[seq_len(n_top)]
    mat_top <- mat_scaled[top_idx, , drop = FALSE]

    # UMAP runs on samples x features
    data_matrix <- t(mat_top)  # samples x genes

    if (!requireNamespace("umap", quietly = TRUE)) {
        stop("Package 'umap' is required for run_umap_tpm().", call. = FALSE)
    }

    # Use umap.defaults so users can reproduce runs
    config <- umap::umap.defaults
    config$n_neighbors <- n_neighbors
    config$min_dist <- min_dist
    config$metric <- metric

    if (!is.null(seed)) set.seed(seed)

    um <- umap::umap(data_matrix, config = config)

    # Build tidy output
    groups <- as.character(meta_samples$condition)
    labels <- if (isTRUE(make_labels)) {
        paste0(meta_samples$sample_names, ".r", meta_samples$rep_id)
    } else {
        colnames(mat_top)
    }

    umap_df <- data.frame(
        UMAP1  = um$layout[, 1],
        UMAP2  = um$layout[, 2],
        sample = colnames(mat_top),
        group  = groups,
        label  = labels,
        stringsAsFactors = FALSE
    )

    out <- list(
        umap = um,
        umap_df = umap_df,
        scaled_data_top = mat_top,
        tpm_cols = meta_samples$tpm_cols,
        sample_table = meta_samples$sample_table,
        groups = levels(meta_samples$condition),
        rep_counts = meta_samples$rep_counts,
        meta = list(
            GeneIDV3_col = GeneIDV3_col,
            TPM_thresh = TPM_thresh,
            scale_data = scale_data,
            n_top_genes_requested = n_top_genes,
            n_top_genes_used = n_top,
            seed = seed,
            n_neighbors = n_neighbors,
            min_dist = min_dist,
            metric = metric
        )
    )

    if (isTRUE(return_input)) {
        out$input_matrix <- mat
    }

    out
}

#' Print UMAP group names for easy copy/paste
#'
#' Prints the group names from a `run_umap_tpm()` result in a format that can
#' be directly copied into `set_umap_group_order()`.
#'
#' @param umap_obj Output from `run_umap_tpm()`
#'
#' @return Invisibly returns the character vector of group names.
#' @export
print_umap_groups <- function(umap_obj) {

    if (!is.list(umap_obj) || is.null(umap_obj$umap_df)) {
        stop("umap_obj must be output from run_umap_tpm()", call. = FALSE)
    }

    groups <- if (is.factor(umap_obj$umap_df$group)) levels(umap_obj$umap_df$group) else unique(umap_obj$umap_df$group)

    cat("c(\n")

    for (i in seq_along(groups)) {
        comma <- if (i < length(groups)) "," else ""
        cat('  "', groups[i], '"', comma, "\n", sep = "")
    }

    cat(")\n")

    invisible(groups)
}


#' Set the order of groups in a UMAP result
#'
#' Reorders the group factor in a `run_umap_tpm()` result so that plots
#' (e.g. `plot_umap()`) display groups in the specified order.
#'
#' @param umap_obj Output from `run_umap_tpm()`
#' @param group_order Character vector specifying desired group order
#'
#' @return Updated UMAP object with reordered groups
#' @export
set_umap_group_order <- function(umap_obj, group_order) {

    if (!is.list(umap_obj) || is.null(umap_obj$umap_df)) {
        stop("umap_obj must be output from run_umap_tpm()", call. = FALSE)
    }

    df <- umap_obj$umap_df

    if (!all(group_order %in% df$group)) {
        missing <- setdiff(group_order, df$group)
        stop("These groups are not present in the UMAP data: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    df$group <- factor(df$group, levels = group_order)

    umap_obj$umap_df <- df
    umap_obj$groups <- group_order

    return(umap_obj)
}


#' Plot UMAP results with numeric group labels and a clean legend
#'
#' Plots UMAP coordinates and labels each sample with a number representing its group.
#' The legend shows the mapping from number -> group name (no redundant "1=" text).
#' Group order follows factor levels if `group` is a factor (e.g., after
#' `set_umap_group_order()`); otherwise uses the order of first appearance.
#'
#' @param x Output from `run_umap_tpm()` (a list with `umap_df`) OR a data.frame
#'   containing columns `UMAP1`, `UMAP2`, and `group`.
#' @param text_size Numeric. Text size for point labels. Default 5.
#' @param theme_fn A ggplot2 theme function. Default `ggplot2::theme_minimal`.
#'
#' @return A ggplot object.
#' @export
plot_umap <- function(x,
                      text_size = 5,
                      theme_fn = ggplot2::theme_minimal) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    # Accept run_umap_tpm() output or a data.frame
    if (is.list(x) && !is.null(x$umap_df)) {
        df <- x$umap_df
    } else if (is.data.frame(x)) {
        df <- x
    } else {
        stop("x must be output from run_umap_tpm() or a data.frame.", call. = FALSE)
    }

    if (!all(c("UMAP1", "UMAP2", "group") %in% names(df))) {
        stop("UMAP data must contain columns 'UMAP1', 'UMAP2', and 'group'.", call. = FALSE)
    }

    # Determine group order:
    # - if factor: use levels (respects set_umap_group_order())
    # - else: use order of first appearance
    if (is.factor(df$group)) {
        group_levels <- levels(df$group)
    } else {
        group_levels <- unique(df$group)
    }

    group_map <- data.frame(
        group = group_levels,
        group_id = seq_along(group_levels),
        stringsAsFactors = FALSE
    )

    # Map group -> numeric id WITHOUT reordering rows (avoid merge())
    df$group_id <- group_map$group_id[match(df$group, group_map$group)]

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = UMAP1,
            y = UMAP2,
            label = group_id,
            color = group
        )
    ) +
        ggplot2::geom_text(size = text_size) +
        ggplot2::scale_color_manual(
            values = rep("black", length(group_levels)),  # keep labels black
            breaks = group_levels,                        # legend order
            labels = group_levels,                        # legend text (no "1=")
            name = "Group"
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(
                override.aes = list(label = group_map$group_id)  # legend keys show 1,2,3...
            )
        ) +
        theme_fn()

    return(p)
}
