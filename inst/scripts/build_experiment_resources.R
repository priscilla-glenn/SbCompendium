# Build SbCompendium ExperimentHub resource files.
#
# Run from the package root. The first argument is an experiment ID or "all";
# the second argument is the output directory. Examples:
#
# Rscript inst/scripts/build_experiment_resources.R nodal_buds_28 resources
# Rscript inst/scripts/build_experiment_resources.R all resources
#
# Large generated RDS files are written to resources/, which is excluded from
# both Git and R package builds.

arguments <- commandArgs(trailingOnly = TRUE)
requested_id <- if (length(arguments) >= 1L) arguments[[1L]] else
    "nodal_buds_28"
resource_dir <- if (length(arguments) >= 2L) arguments[[2L]] else
    "resources"

required_packages <- c(
    "readxl",
    "S4Vectors",
    "SummarizedExperiment"
)
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop(
        "Install required preparation packages before building resources: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

package_description <- file.path("DESCRIPTION")
package_name <- if (file.exists(package_description)) {
    unname(read.dcf(package_description, fields = "Package")[1L, 1L])
} else {
    NA_character_
}
if (!identical(package_name, "SbCompendium")) {
    stop("Run this script from the SbCompendium package root.", call. = FALSE)
}

source(file.path("inst", "scripts", "validate_experiment_resources.R"))

manifest_path <- file.path(
    "inst",
    "schema",
    "Draft_ExperimentHub_Manifest.xlsx"
)
compendium_path <- file.path("data", "sorghum_compendium.rda")

if (!file.exists(manifest_path)) {
    stop("Manifest not found: ", manifest_path, call. = FALSE)
}
if (!file.exists(compendium_path)) {
    stop("Compendium data not found: ", compendium_path, call. = FALSE)
}

hub_manifest <- as.data.frame(
    readxl::read_excel(manifest_path, sheet = "Hub Manifest"),
    stringsAsFactors = FALSE
)
sample_manifest <- as.data.frame(
    readxl::read_excel(manifest_path, sheet = "Sample Manifest"),
    stringsAsFactors = FALSE
)

load_environment <- new.env(parent = emptyenv())
loaded_names <- load(compendium_path, envir = load_environment)
if (!"sorghum_compendium" %in% loaded_names) {
    stop(
        "The compendium RDA does not contain 'sorghum_compendium'.",
        call. = FALSE
    )
}
compendium <- load_environment$sorghum_compendium
if (!is.list(compendium)) {
    stop("sorghum_compendium must be a named list.", call. = FALSE)
}

.clean_placeholder <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x)) | toupper(trimws(x)) == "TODO"] <- NA
    x
}

.first_existing_column <- function(data, candidates) {
    match <- intersect(candidates, names(data))
    if (length(match)) match[[1L]] else NULL
}

.gene_level_matrix <- function(data, value_columns) {
    transcript_column <- .first_existing_column(
        data,
        c("TranscriptID", "TranscriptIDV3")
    )
    if (is.null(transcript_column)) {
        stop(
            "Expression table lacks TranscriptID or TranscriptIDV3.",
            call. = FALSE
        )
    }
    if (!length(value_columns)) {
        stop("No expression columns were supplied.", call. = FALSE)
    }

    expression_data <- data[value_columns]
    numeric_columns <- vapply(expression_data, is.numeric, logical(1))
    if (!all(numeric_columns)) {
        stop(
            "Non-numeric expression columns: ",
            paste(value_columns[!numeric_columns], collapse = ", "),
            call. = FALSE
        )
    }

    transcript_ids <- as.character(data[[transcript_column]])
    if (anyNA(transcript_ids) || any(!nzchar(transcript_ids))) {
        stop("Transcript identifiers contain missing values.", call. = FALSE)
    }
    gene_ids <- sub("\\.[0-9]+$", "", transcript_ids)
    matrix_data <- as.matrix(expression_data)
    storage.mode(matrix_data) <- "double"
    rowsum(matrix_data, group = gene_ids, reorder = FALSE, na.rm = FALSE)
}

.sample_components <- function(tpm_columns) {
    sample_id <- sub("_TPM$", "", tpm_columns)
    replicate_pattern <- "([._-][0-9]+|[._-]?R[0-9]+)$"
    has_replicate <- grepl(replicate_pattern, sample_id)
    condition <- sub(replicate_pattern, "", sample_id)
    replicate <- rep(NA_character_, length(sample_id))
    replicate[has_replicate] <- sub(
        paste0("^.*", replicate_pattern),
        "\\1",
        sample_id[has_replicate]
    )
    replicate <- sub("^[._-]", "", replicate)
    replicate <- sub("^R", "", replicate)
    data.frame(
        sample_id = sample_id,
        condition = condition,
        replicate = replicate,
        stringsAsFactors = FALSE,
        row.names = sample_id
    )
}

.prepare_gene_metadata <- function(gene_ids, annotation) {
    base <- data.frame(
        gene_id = gene_ids,
        stringsAsFactors = FALSE,
        row.names = gene_ids
    )
    if (is.null(annotation)) {
        return(S4Vectors::DataFrame(base))
    }

    annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
    gene_column <- .first_existing_column(
        annotation,
        c("GeneIDV3", "GeneID", "gene_id")
    )
    annotation_ids <- if (is.null(gene_column))
        rownames(annotation) else as.character(annotation[[gene_column]])
    if (is.null(annotation_ids)) {
        warning("Gene annotation has no recognizable gene identifiers.")
        return(S4Vectors::DataFrame(base))
    }

    keep <- !is.na(annotation_ids) & nzchar(annotation_ids) &
        !duplicated(annotation_ids)
    annotation <- annotation[keep, , drop = FALSE]
    annotation_ids <- annotation_ids[keep]
    if (!is.null(gene_column)) {
        annotation[[gene_column]] <- NULL
    }
    aligned <- annotation[match(gene_ids, annotation_ids), , drop = FALSE]
    rownames(aligned) <- gene_ids
    S4Vectors::DataFrame(cbind(base, aligned, stringsAsFactors = FALSE))
}

.prepare_sample_metadata <- function(
    tpm_columns,
    experiment_id,
    sample_manifest
) {
    sample_data <- .sample_components(tpm_columns)
    sample_data$experiment_id <- experiment_id

    source_sample_column <- .first_existing_column(
        sample_manifest,
        c("SampleID", "Sample Name")
    )
    if (is.null(source_sample_column)) {
        warning("Sample manifest has no SampleID or Sample Name column.")
        return(S4Vectors::DataFrame(sample_data))
    }

    experiment_rows <- sample_manifest[
        as.character(sample_manifest$ExperimentID) == experiment_id,
        ,
        drop = FALSE
    ]
    source_ids <- as.character(experiment_rows[[source_sample_column]])
    matched <- match(sample_data$condition, source_ids)
    descriptive_columns <- setdiff(
        names(experiment_rows),
        c("SampleID", "ExperimentID", "MappingStatus", "Sample Name")
    )
    if (length(descriptive_columns)) {
        descriptive_data <- experiment_rows[
            matched,
            descriptive_columns,
            drop = FALSE
        ]
        names(descriptive_data) <- make.names(
            names(descriptive_data),
            unique = TRUE
        )
        rownames(descriptive_data) <- rownames(sample_data)
        sample_data <- cbind(sample_data, descriptive_data)
    }

    sample_data$manifest_match <- !is.na(matched)
    if (any(!sample_data$manifest_match)) {
        warning(
            sum(!sample_data$manifest_match),
            " assay samples did not match a condition in Sample Manifest."
        )
    }
    S4Vectors::DataFrame(sample_data)
}

.build_one_experiment <- function(
    experiment_id,
    expression_table,
    experiment_manifest_row,
    sample_manifest,
    annotation = NULL
) {
    tpm_columns <- grep("_TPM$", names(expression_table), value = TRUE)
    count_columns <- grep("_counts$", names(expression_table), value = TRUE)
    if (!length(tpm_columns)) {
        stop(experiment_id, " has no replicate-level _TPM columns.",
             call. = FALSE)
    }

    tpm <- .gene_level_matrix(expression_table, tpm_columns)
    sample_ids <- sub("_TPM$", "", tpm_columns)
    colnames(tpm) <- sample_ids
    assays <- list(TPM = tpm)

    if (length(count_columns)) {
        count_ids <- sub("_counts$", "", count_columns)
        missing_counts <- setdiff(sample_ids, count_ids)
        extra_counts <- setdiff(count_ids, sample_ids)
        if (length(missing_counts) || length(extra_counts)) {
            stop(
                experiment_id,
                " has unmatched TPM/count sample IDs. Missing counts: ",
                paste(missing_counts, collapse = ", "),
                "; extra counts: ",
                paste(extra_counts, collapse = ", "),
                call. = FALSE
            )
        }
        ordered_count_columns <- count_columns[match(sample_ids, count_ids)]
        counts <- .gene_level_matrix(
            expression_table,
            ordered_count_columns
        )
        colnames(counts) <- sample_ids
        if (!identical(rownames(counts), rownames(tpm))) {
            stop("TPM and count gene order differs for ", experiment_id,
                 call. = FALSE)
        }
        assays$counts <- counts
    }

    gene_data <- .prepare_gene_metadata(rownames(tpm), annotation)
    sample_data <- .prepare_sample_metadata(
        tpm_columns,
        experiment_id,
        sample_manifest
    )

    source_project_id <- .clean_placeholder(
        experiment_manifest_row$JGIProjectId[[1L]]
    )
    genome_build <- .clean_placeholder(
        experiment_manifest_row$Genome[[1L]]
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        rowData = gene_data,
        colData = sample_data,
        metadata = list(
            experiment_id = experiment_id,
            title = experiment_manifest_row$Title[[1L]],
            source_project_id = source_project_id,
            genome_build = genome_build,
            processing_version = experiment_manifest_row$SourceVersion[[1L]],
            source_manifest = experiment_manifest_row
        )
    )
}

annotation_name <- intersect(
    c("Sorghum_gene_annotation", "sorghum_gene_annotation"),
    names(compendium)
)
annotation <- if (length(annotation_name))
    compendium[[annotation_name[[1L]]]] else NULL

available_ids <- intersect(
    as.character(hub_manifest$ResourceId),
    names(compendium)
)
requested_ids <- if (identical(tolower(requested_id), "all")) {
    unavailable <- setdiff(as.character(hub_manifest$ResourceId),
                           names(compendium))
    if (length(unavailable)) {
        warning(
            "Skipping manifest resources absent from the current compendium: ",
            paste(unavailable, collapse = ", ")
        )
    }
    available_ids
} else {
    requested_id
}

if (!length(requested_ids)) {
    stop("No buildable experiment resources were selected.", call. = FALSE)
}
missing_from_compendium <- setdiff(requested_ids, names(compendium))
if (length(missing_from_compendium)) {
    stop(
        "Experiment is absent from the current compendium: ",
        paste(missing_from_compendium, collapse = ", "),
        call. = FALSE
    )
}

dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)

for (experiment_id in requested_ids) {
    message("Building ", experiment_id, "...")
    manifest_row <- hub_manifest[
        as.character(hub_manifest$ResourceId) == experiment_id,
        ,
        drop = FALSE
    ]
    if (nrow(manifest_row) != 1L) {
        stop(
            "Expected one Hub Manifest row for ",
            experiment_id,
            "; found ",
            nrow(manifest_row),
            ".",
            call. = FALSE
        )
    }

    resource <- .build_one_experiment(
        experiment_id = experiment_id,
        expression_table = compendium[[experiment_id]],
        experiment_manifest_row = manifest_row,
        sample_manifest = sample_manifest,
        annotation = annotation
    )
    validate_experiment_resource(
        resource,
        expected_experiment_id = experiment_id,
        publication_ready = FALSE
    )

    output_path <- file.path(resource_dir, paste0(experiment_id, ".rds"))
    saveRDS(resource, output_path, compress = "xz")
    validate_resource_file(
        output_path,
        expected_experiment_id = experiment_id,
        publication_ready = FALSE
    )
    message("Saved and validated ", normalizePath(output_path))
}

message("Finished building ", length(requested_ids), " resource(s).")
