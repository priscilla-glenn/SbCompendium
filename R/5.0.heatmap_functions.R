#' Build a heatmap export table with TPMmean expression, optional clusters,
#' and optional FC + annotations
#'
#' Merges TPMmean expression data (from `gene_expr`, GeneIDV3 in rownames),
#' optional k-means cluster assignments, and optionally FC/DE columns and
#' annotation columns into a single table suitable for heatmaps and exporting.
#'
#' You can export either:
#' - all TPMmean columns (default), OR
#' - only the TPMmean columns corresponding to two DESeq2 group names
#'   (pass `groups = c(groupA, groupB)`).
#'
#' @param gene_expr data.frame/matrix with GeneIDV3s as rownames and TPMmean columns.
#' @param clusters Optional. Output from run_kmeans(), a data.frame with GeneIDV3 + cluster,
#'   or a named vector of cluster assignments. If NULL, cluster column is omitted.
#' @param top_de Optional data.frame containing GeneIDV3 and FC/DE columns.
#' @param annot_df Optional data.frame containing GeneIDV3 and annotation columns.
#' @param GeneIDV3_col Column name of gene IDs in `top_de` / `annot_df`. Default "GeneIDV3".
#' @param fc_cols Optional vector of FC columns to include from `top_de`.
#'   If NULL, uses common DESeq2 columns that exist.
#' @param annot_cols Optional vector of annotation columns to include from `annot_df`.
#'   If NULL, includes all columns except GeneIDV3_col.
#' @param tpmmean_cols Optional TPMmean columns. If NULL uses `.get_tpmmean_cols(gene_expr)`.
#' @param groups Optional character vector length 2: c(groupA, groupB).
#'   If provided, restricts TPMmean columns to those whose column names contain
#'   either groupA or groupB (fixed substring match).
#' @param cluster_col_name Name of cluster column (if clusters supplied). Default "cluster".
#' @param sort_by_cluster Logical. If TRUE (default) and clusters provided, sorts by cluster.
#' @param secondary_sort Character. One of:
#'   - "GeneIDV3"
#'   - "abs_log2FoldChange" (only if log2FoldChange exists)
#'   - any existing column name in the output
#'
#' @return data.frame with GeneIDV3, optional cluster, optional FC, optional annotations,
#'   TPMmean columns (all or group-filtered).
#' @export
#'
#'
#' @examples
#' # gene_expr: GeneIDV3s in rownames + TPMmean columns
#' gene_expr <- data.frame(
#'   INT1_TPMmean = c(10, 0, 5, 2),
#'   INT2_TPMmean = c(12, 1, 4, 8),
#'   INT3_TPMmean = c(3,  8, 6, 2)
#' )
#' rownames(gene_expr) <- c(
#'   "Sobic.001G000100",
#'   "Sobic.001G000200",
#'   "Sobic.001G000300",
#'   "Sobic.001G000400"
#' )
#'
#' # 1) Expression-only heatmap table
#' hm1 <- build_heatmap(gene_expr)
#' hm1
#'
#' # Differential expression table (optional)
#' top_de <- data.frame(
#'   GeneIDV3 = c(
#'     "Sobic.001G000100",
#'     "Sobic.001G000200",
#'     "Sobic.001G000300"
#'   ),
#'   log2FoldChange = c(2.0, -1.2, 3.1),
#'   padj = c(0.01, 0.20, 0.03),
#'   stringsAsFactors = FALSE
#' )
#'
#' # 2) FC + expression
#' hm_fc <- build_heatmap(
#'   gene_expr,
#'   top_de = top_de
#' )
#' hm_fc
#'
#' # Cluster assignments
#' clusters_df <- data.frame(
#'   GeneIDV3 = c(
#'     "Sobic.001G000100",
#'     "Sobic.001G000200",
#'     "Sobic.001G000300"
#'   ),
#'   cluster = c(1, 1, 2),
#'   stringsAsFactors = FALSE
#' )
#'
#' # 3) Cluster + expression
#' hm_cluster <- build_heatmap(
#'   gene_expr,
#'   clusters = clusters_df
#' )
#' hm_cluster
#'
#' # Annotation table (optional)
#' annot_df <- data.frame(
#'   GeneIDV3 = c(
#'     "Sobic.001G000100",
#'     "Sobic.001G000200",
#'     "Sobic.001G000300",
#'     "Sobic.001G000400"
#'   ),
#'   NAME = c("GeneA","GeneB","GeneC","GeneD"),
#'   Functional_category = c("TF","Metabolism","Transport","Unknown"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # 4) FC + cluster + annotations + expression
#' hm_full <- build_heatmap(
#'   gene_expr,
#'   clusters = clusters_df,
#'   top_de = top_de,
#'   annot_df = annot_df,
#'   fc_cols = c("log2FoldChange", "padj"),
#'   annot_cols = c("NAME", "Functional_category"),
#'   sort_by_cluster = TRUE,
#'   secondary_sort = "abs_log2FoldChange"
#' )
#' hm_full
#'
#' # 5) Restrict TPMmean columns to two groups (substring match)
#' hm_groups <- build_heatmap(
#'   gene_expr,
#'   top_de = top_de,
#'   groups = c("INT1", "INT2")
#' )
#' hm_groups
#'
#' # 6) Clusters can also be provided as a named vector
#' clusters_vec <- c(1, 2, 2)
#' names(clusters_vec) <- c(
#'   "Sobic.001G000100",
#'   "Sobic.001G000200",
#'   "Sobic.001G000300"
#' )
#'
#' hm_vec <- build_heatmap(
#'   gene_expr,
#'   clusters = clusters_vec
#' )
#' hm_vec
#'
build_heatmap <- function(gene_expr,
                          clusters = NULL,
                          top_de = NULL,
                          annot_df = NULL,
                          GeneIDV3_col = "GeneIDV3",
                          fc_cols = NULL,
                          annot_cols = NULL,
                          tpmmean_cols = NULL,
                          groups = NULL,
                          cluster_col_name = "cluster",
                          sort_by_cluster = TRUE,
                          secondary_sort = c("GeneIDV3", "abs_log2FoldChange")) {

    secondary_sort <- match.arg(secondary_sort)

    if (!is.data.frame(gene_expr) && !is.matrix(gene_expr)) {
        stop("gene_expr must be a data.frame or matrix.", call. = FALSE)
    }

    GeneIDV3s <- rownames(gene_expr)
    if (is.null(GeneIDV3s) || anyNA(GeneIDV3s) || any(GeneIDV3s == "")) {
        stop("gene_expr must contain GeneIDV3s as non-empty rownames.", call. = FALSE)
    }

    # TPMmean columns (all by default)
    if (is.null(tpmmean_cols)) {
        tpmmean_cols <- .get_tpmmean_cols(gene_expr)
    } else {
        missing <- setdiff(tpmmean_cols, colnames(gene_expr))
        if (length(missing) > 0) {
            stop("These tpmmean_cols are not in gene_expr: ",
                 paste(missing, collapse = ", "), call. = FALSE)
        }
    }

    # Optional: restrict TPMmean columns to only those for the DESeq2 groups
    if (!is.null(groups)) {
        if (!is.character(groups) || length(groups) != 2) {
            stop("groups must be a character vector of length 2: c(groupA, groupB).", call. = FALSE)
        }
        groupA <- groups[1]
        groupB <- groups[2]

        keep <- grepl(groupA, tpmmean_cols, fixed = TRUE) | grepl(groupB, tpmmean_cols, fixed = TRUE)
        tpmmean_cols <- tpmmean_cols[keep]

        if (length(tpmmean_cols) == 0) {
            stop(
                "groups filtering removed all TPMmean columns. ",
                "I looked for TPMmean colnames containing: '", groupA, "' or '", groupB, "'. ",
                "Either pass tpmmean_cols manually, or adjust naming/matching.",
                call. = FALSE
            )
        }
    }

    # Expression df (base)
    out <- as.data.frame(gene_expr[, tpmmean_cols, drop = FALSE])
    out$GeneIDV3 <- rownames(out)

    # Optional clusters
    if (!is.null(clusters)) {

        # Normalize clusters input -> data.frame(GeneIDV3, cluster_col_name)
        if (is.list(clusters) && !is.null(clusters$clusters)) {
            cluster_df <- clusters$clusters
        } else if (is.data.frame(clusters)) {
            cluster_df <- clusters
        } else if (is.atomic(clusters) && !is.null(names(clusters))) {
            cluster_df <- data.frame(
                GeneIDV3 = names(clusters),
                cluster = as.integer(clusters),
                stringsAsFactors = FALSE
            )
        } else {
            stop("clusters must be run_kmeans() output, a data.frame, a named vector, or NULL.",
                 call. = FALSE)
        }

        if (!("GeneIDV3" %in% names(cluster_df))) {
            stop("clusters must contain a GeneIDV3 column.", call. = FALSE)
        }
        if (!("cluster" %in% names(cluster_df)) && !(cluster_col_name %in% names(cluster_df))) {
            stop("clusters must contain a cluster column named 'cluster' or cluster_col_name.",
                 call. = FALSE)
        }

        # Rename cluster column to cluster_col_name
        if ("cluster" %in% names(cluster_df) && cluster_col_name != "cluster") {
            names(cluster_df)[names(cluster_df) == "cluster"] <- cluster_col_name
        }

        cluster_df <- cluster_df[, c("GeneIDV3", cluster_col_name), drop = FALSE]

        # Merge clusters onto expression (inner join keeps only genes present in both)
        out <- merge(cluster_df, out, by = "GeneIDV3")
    }

    # Optional FC data
    if (!is.null(top_de)) {
        if (!is.data.frame(top_de)) stop("top_de must be a data.frame.", call. = FALSE)
        if (!GeneIDV3_col %in% names(top_de)) {
            stop("top_de missing GeneIDV3_col = '", GeneIDV3_col, "'.", call. = FALSE)
        }

        if (is.null(fc_cols)) {
            candidates <- c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
            fc_cols <- intersect(candidates, names(top_de))
        } else {
            missing <- setdiff(fc_cols, names(top_de))
            if (length(missing) > 0) {
                stop("fc_cols not in top_de: ", paste(missing, collapse = ", "), call. = FALSE)
            }
        }

        if (length(fc_cols) > 0) {
            de_sub <- top_de[, c(GeneIDV3_col, fc_cols), drop = FALSE]
            names(de_sub)[names(de_sub) == GeneIDV3_col] <- "GeneIDV3"
            out <- merge(out, de_sub, by = "GeneIDV3", all.x = TRUE)
        }
    }

    # Optional annotations
    if (!is.null(annot_df)) {
        if (!is.data.frame(annot_df)) stop("annot_df must be a data.frame.", call. = FALSE)
        if (!GeneIDV3_col %in% names(annot_df)) {
            stop("annot_df missing GeneIDV3_col = '", GeneIDV3_col, "'.", call. = FALSE)
        }

        if (is.null(annot_cols)) {
            annot_cols <- setdiff(names(annot_df), GeneIDV3_col)
        } else {
            missing <- setdiff(annot_cols, names(annot_df))
            if (length(missing) > 0) {
                stop("annot_cols not in annot_df: ", paste(missing, collapse = ", "), call. = FALSE)
            }
            annot_cols <- setdiff(annot_cols, GeneIDV3_col)
        }

        if (length(annot_cols) > 0) {
            ann_sub <- annot_df[, c(GeneIDV3_col, annot_cols), drop = FALSE]
            names(ann_sub)[names(ann_sub) == GeneIDV3_col] <- "GeneIDV3"
            out <- merge(out, ann_sub, by = "GeneIDV3", all.x = TRUE)
        }
    }

    # Column order
    has_cluster <- (!is.null(clusters) && cluster_col_name %in% names(out))
    core <- c("GeneIDV3", if (has_cluster) cluster_col_name else NULL)
    expr <- tpmmean_cols

    fc_present <- c()
    if (!is.null(top_de)) {
        fc_present <- if (is.null(fc_cols)) {
            intersect(c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), names(out))
        } else {
            intersect(fc_cols, names(out))
        }
    }

    annot_present <- c()
    if (!is.null(annot_df)) {
        annot_present <- if (is.null(annot_cols)) {
            intersect(setdiff(names(annot_df), GeneIDV3_col), names(out))
        } else {
            intersect(setdiff(annot_cols, GeneIDV3_col), names(out))
        }
    }

    keep_cols <- c(core, fc_present, annot_present, expr)
    keep_cols <- keep_cols[keep_cols %in% names(out)]
    out <- out[, keep_cols, drop = FALSE]

    # Sorting
    if (has_cluster && isTRUE(sort_by_cluster)) {
        if (secondary_sort == "abs_log2FoldChange" && "log2FoldChange" %in% names(out)) {
            ord <- order(out[[cluster_col_name]], -abs(out[["log2FoldChange"]]), out$GeneIDV3)
        } else if (secondary_sort %in% names(out)) {
            ord <- order(out[[cluster_col_name]], out[[secondary_sort]], out$GeneIDV3)
        } else {
            ord <- order(out[[cluster_col_name]], out$GeneIDV3)
        }
        out <- out[ord, , drop = FALSE]
    } else {
        # No clusters: sort by secondary (or GeneIDV3)
        if (secondary_sort == "abs_log2FoldChange" && "log2FoldChange" %in% names(out)) {
            ord <- order(-abs(out[["log2FoldChange"]]), out$GeneIDV3)
        } else if (secondary_sort %in% names(out) && secondary_sort != "GeneIDV3") {
            ord <- order(out[[secondary_sort]], out$GeneIDV3)
        } else {
            ord <- order(out$GeneIDV3)
        }
        out <- out[ord, , drop = FALSE]
    }

    rownames(out) <- out$GeneIDV3
    out
}

#' Write heatmap export table to Excel using Python (robust UTF-8).
#'
#' @param export_df data.frame to export.
#' @param file Output .xlsx file path.
#' @param sheet Sheet name.
#' @param expr_cols Expression columns to color. If NULL, auto-detect columns
#'   containing "TPMmean".
#' @param cluster_col_name Cluster column name, used only for reference (optional).
#' @param replace_zero_annotation_with_blank Replace 0 with blank in
#'   non-TPMmean columns.
#'
#' @return Invisibly returns `file`.
#' @export
write_heatmap_xlsx <- function(export_df,
                               file,
                               sheet = "Sheet1",
                               expr_cols = NULL,
                               cluster_col_name = "cluster",
                               replace_zero_annotation_with_blank = TRUE) {

    stopifnot(is.data.frame(export_df))

    # --- Safe UTF-8 conversion ---
    export_df <- as.data.frame(export_df, stringsAsFactors = FALSE)
    export_df[] <- lapply(export_df, function(x) {
        if (is.character(x)) iconv(x, from = "", to = "UTF-8", sub = "byte") else x
    })

    # Ensure GeneIDV3 exists
    if (!"GeneIDV3" %in% colnames(export_df)) {
        export_df$GeneIDV3 <- rownames(export_df)
    }

    # --- Check reticulate package ---
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required.", call. = FALSE)
    }

    # --- Check required Python modules ---
    required_modules <- c("pandas", "xlsxwriter", "openpyxl")
    missing_modules <- required_modules[!vapply(required_modules,
                                                reticulate::py_module_available,
                                                logical(1))]

    if (length(missing_modules) > 0) {
        stop(
            paste0(
                "Missing required Python modules: ",
                paste(missing_modules, collapse = ", "),
                ". Please install them in your Python environment."
            ),
            call. = FALSE
        )
    }

    # --- Locate Python script ---
    py_file <- system.file("python", "export_heatmap_xlsx.py",
                           package = "SbCompendium")
    if (py_file == "") {
        stop("Python script 'export_heatmap_xlsx.py' not found in package.",
             call. = FALSE)
    }

    # --- Load Python script ---
    reticulate::source_python(py_file)

    # --- Call Python function via reticulate ---
    reticulate::py$py_write_heatmap_xlsx(
        export_df = export_df,
        file = normalizePath(file, mustWork = FALSE),
        sheet = sheet,
        expr_cols = expr_cols,
        cluster_col_name = cluster_col_name,
        replace_zero_annotation_with_blank = replace_zero_annotation_with_blank
    )

    invisible(file)
}
