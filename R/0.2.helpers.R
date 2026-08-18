# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

# @keywords internal
.get_annotation <- function() {
    sorghum_compendium[["Sorghum_gene_annotation"]]
}

# @keywords internal
.get_tpmmean_cols <- function(df, tpm_pattern = "_TPMmean$") {
    cols <- grep(tpm_pattern, colnames(df), value = TRUE)
    if (length(cols) == 0) {
        stop("No TPM columns matched pattern: ", tpm_pattern, call. = FALSE)
    }
    cols
}

#' Calculate mean expression for sample groups
#'
#' Calculates row means from a named assay of a `SummarizedExperiment`,
#' grouping samples by a field in `colData`. This replaces stored `_TPMmean`
#' columns: condition means are derived reproducibly when they are needed.
#'
#' @param x A `SummarizedExperiment`.
#' @param group_by Name of a `colData` field defining the sample groups.
#' @param assay_name Name of the assay to summarize. Default is `"TPM"`.
#' @param na.rm Logical; remove missing values when calculating means.
#'
#' @return A numeric gene-by-group matrix.
#' @export
summarize_tpm <- function(x,
                          group_by = "condition",
                          assay_name = "TPM",
                          na.rm = FALSE) {
    if (!methods::is(x, "SummarizedExperiment")) {
        stop("x must be a SummarizedExperiment.", call. = FALSE)
    }
    if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
        stop("Assay not found: ", assay_name, call. = FALSE)
    }

    sample_data <- SummarizedExperiment::colData(x)
    if (!group_by %in% names(sample_data)) {
        stop("colData field not found: ", group_by, call. = FALSE)
    }

    groups <- as.character(sample_data[[group_by]])
    if (anyNA(groups) || any(!nzchar(groups))) {
        stop("Grouping field contains missing or empty values.",
             call. = FALSE)
    }
    groups <- factor(groups, levels = unique(groups))

    expression <- SummarizedExperiment::assay(x, assay_name)
    group_indices <- split(seq_len(ncol(expression)), groups, drop = TRUE)
    means <- vapply(
        group_indices,
        function(index) {
            rowMeans(expression[, index, drop = FALSE], na.rm = na.rm)
        },
        numeric(nrow(expression))
    )
    rownames(means) <- rownames(x)
    means
}

#' Infer groups and replicate IDs from *_TPM sample column names
#'
#' Internal helper used by UMAP and other functions. Removes replicate suffixes
#' like .1, _1, -1, .R1, _R2, -R3 that appear immediately before `_TPM`.
#'
#' @param tpm_cols Character vector of column names ending in `_TPM`.
#' @return A list with sample_names, condition (factor), rep_counts, rep_id, and sample_table.
#' @keywords internal
.infer_groups_from_tpm_cols <- function(tpm_cols) {

    # Group name = base sample name without replicate suffix and "_TPM"
    sample_names <- sub("([._-][0-9]+|[._-]?R[0-9]+)_TPM$", "", tpm_cols)

    # Replicate id (best-effort; useful for labels)
    rep_id <- sub("^.*([._-][0-9]+|[._-]?R[0-9]+)_TPM$", "\\1", tpm_cols)
    rep_id <- sub("^[._-]", "", rep_id)   # drop leading separator
    rep_id <- sub("^R", "", rep_id)       # drop leading R if present

    condition <- factor(sample_names, levels = unique(sample_names))

    sample_table <- data.frame(
        sample = tpm_cols,
        group  = as.character(condition),
        rep    = rep_id,
        stringsAsFactors = FALSE
    )

    list(
        tpm_cols = tpm_cols,
        sample_names = sample_names,
        condition = condition,
        rep_counts = table(sample_names),
        rep_id = rep_id,
        sample_table = sample_table
    )
}
