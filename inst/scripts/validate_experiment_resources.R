# Validation helpers for SbCompendium ExperimentHub resources.
#
# These functions are used while preparing data and are not part of the
# installed SbCompendium user API. Source this file from the package root:
#
# source("inst/scripts/validate_experiment_resources.R")

.nonempty_scalar <- function(x) {
    length(x) == 1L && !is.na(x) && nzchar(trimws(as.character(x)))
}

.assert_resource <- function(ok, message) {
    if (!isTRUE(ok)) {
        stop(message, call. = FALSE)
    }
}

#' Validate one prepared SbCompendium experiment
#'
#' @param se A SummarizedExperiment object.
#' @param expected_experiment_id Optional experiment identifier expected in
#'   `metadata(se)$experiment_id` and `colData(se)$experiment_id`.
#' @param publication_ready If `TRUE`, also reject placeholder or missing
#'   genome, source, and processing metadata.
#' @return `se`, invisibly, when all checks pass.
validate_experiment_resource <- function(
    se,
    expected_experiment_id = NULL,
    publication_ready = FALSE
) {
    .assert_resource(
        methods::is(se, "SummarizedExperiment"),
        "Resource must inherit from SummarizedExperiment."
    )

    assay_names <- SummarizedExperiment::assayNames(se)
    .assert_resource(
        "TPM" %in% assay_names,
        "Resource must contain an assay named 'TPM'."
    )

    .assert_resource(nrow(se) > 0L, "Resource contains no genes.")
    .assert_resource(ncol(se) > 0L, "Resource contains no samples.")
    .assert_resource(
        !is.null(rownames(se)) && !anyNA(rownames(se)) &&
            all(nzchar(rownames(se))),
        "Every gene must have a non-empty identifier."
    )
    .assert_resource(
        !anyDuplicated(rownames(se)),
        "Gene identifiers must be unique."
    )
    .assert_resource(
        !is.null(colnames(se)) && !anyNA(colnames(se)) &&
            all(nzchar(colnames(se))),
        "Every sample must have a non-empty identifier."
    )
    .assert_resource(
        !anyDuplicated(colnames(se)),
        "Sample identifiers must be unique."
    )

    sample_data <- SummarizedExperiment::colData(se)
    .assert_resource(
        identical(colnames(se), rownames(sample_data)),
        "Assay columns must match colData row names and order."
    )

    required_sample_fields <- c(
        "sample_id",
        "experiment_id",
        "condition",
        "replicate"
    )
    missing_sample_fields <- setdiff(
        required_sample_fields,
        names(sample_data)
    )
    .assert_resource(
        length(missing_sample_fields) == 0L,
        paste0(
            "Missing required colData fields: ",
            paste(missing_sample_fields, collapse = ", ")
        )
    )
    .assert_resource(
        identical(as.character(sample_data$sample_id), colnames(se)),
        "colData(se)$sample_id must match assay column names and order."
    )

    gene_data <- SummarizedExperiment::rowData(se)
    .assert_resource(
        "gene_id" %in% names(gene_data),
        "rowData(se) must contain a 'gene_id' field."
    )
    .assert_resource(
        identical(as.character(gene_data$gene_id), rownames(se)),
        "rowData(se)$gene_id must match assay row names and order."
    )

    tpm <- SummarizedExperiment::assay(se, "TPM")
    .assert_resource(is.numeric(tpm), "TPM assay must be numeric.")
    .assert_resource(
        all(is.finite(tpm)),
        "TPM assay contains missing or non-finite values."
    )
    .assert_resource(all(tpm >= 0), "TPM assay contains negative values.")

    if ("counts" %in% assay_names) {
        counts <- SummarizedExperiment::assay(se, "counts")
        .assert_resource(is.numeric(counts), "Counts assay must be numeric.")
        .assert_resource(
            identical(dim(counts), dim(tpm)),
            "Counts and TPM assays must have identical dimensions."
        )
        .assert_resource(
            identical(dimnames(counts), dimnames(tpm)),
            "Counts and TPM assays must have identical dimnames."
        )
        .assert_resource(
            all(is.finite(counts)),
            "Counts assay contains missing or non-finite values."
        )
        .assert_resource(
            all(counts >= 0),
            "Counts assay contains negative values."
        )
        .assert_resource(
            all(abs(counts - round(counts)) < sqrt(.Machine$double.eps)),
            "Counts assay must contain integer-like values."
        )
    }

    object_metadata <- S4Vectors::metadata(se)
    experiment_id <- object_metadata$experiment_id
    .assert_resource(
        .nonempty_scalar(experiment_id),
        "metadata(se)$experiment_id must be one non-empty value."
    )
    .assert_resource(
        all(as.character(sample_data$experiment_id) == experiment_id),
        "Every colData experiment_id must match metadata(se)$experiment_id."
    )

    if (!is.null(expected_experiment_id)) {
        .assert_resource(
            identical(as.character(experiment_id),
                      as.character(expected_experiment_id)),
            "The resource experiment ID does not match the expected ID."
        )
    }

    if (isTRUE(publication_ready)) {
        required_metadata <- c(
            "genome_build",
            "source_project_id",
            "processing_version"
        )
        missing_metadata <- required_metadata[
            !vapply(
                object_metadata[required_metadata],
                .nonempty_scalar,
                logical(1)
            )
        ]
        placeholder_metadata <- required_metadata[
            vapply(
                object_metadata[required_metadata],
                function(x) identical(toupper(trimws(as.character(x))),
                                      "TODO"),
                logical(1)
            )
        ]
        .assert_resource(
            length(c(missing_metadata, placeholder_metadata)) == 0L,
            paste0(
                "Publication metadata is incomplete: ",
                paste(
                    unique(c(missing_metadata, placeholder_metadata)),
                    collapse = ", "
                )
            )
        )
    }

    invisible(se)
}

#' Validate an RDS resource after serialization
#'
#' @param path Path to an RDS file.
#' @param expected_experiment_id Optional expected experiment identifier.
#' @param publication_ready Passed to `validate_experiment_resource()`.
#' @return The reloaded resource, invisibly.
validate_resource_file <- function(
    path,
    expected_experiment_id = NULL,
    publication_ready = FALSE
) {
    .assert_resource(file.exists(path), paste0("File does not exist: ", path))
    resource <- readRDS(path)
    validate_experiment_resource(
        resource,
        expected_experiment_id = expected_experiment_id,
        publication_ready = publication_ready
    )
    invisible(resource)
}
