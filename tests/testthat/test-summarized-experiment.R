test_that("summarize_tpm averages samples by condition", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            TPM = matrix(c(1, 3, 5, 7, 2, 4, 6, 8), nrow = 2),
            counts = matrix(1:8, nrow = 2)
        ),
        colData = S4Vectors::DataFrame(
            condition = c("control", "control", "treated", "treated")
        )
    )
    rownames(se) <- c("gene1", "gene2")
    colnames(se) <- paste0("sample", 1:4)

    result <- summarize_tpm(se)

    expect_equal(colnames(result), c("control", "treated"))
    expect_equal(unname(result[, "control"]), c(3, 5))
    expect_equal(unname(result[, "treated"]), c(4, 6))
})

test_that("condition means feed clustering and heatmap preparation", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            TPM = matrix(c(1, 3, 5, 7, 2, 4, 6, 8), nrow = 2),
            counts = matrix(1:8, nrow = 2)
        ),
        rowData = S4Vectors::DataFrame(symbol = c("A", "B")),
        colData = S4Vectors::DataFrame(
            condition = c("control", "control", "treated", "treated")
        )
    )
    rownames(se) <- c("gene1", "gene2")
    colnames(se) <- paste0("sample", 1:4)
    top_de <- data.frame(GeneIDV3 = rownames(se))

    prepared <- prepare_kmeans_matrix_topDE_TPMmean(
        top_de, se, scale_rows_if_gt2 = FALSE, log_transform_if_le2 = FALSE
    )
    expect_equal(colnames(prepared), c("control", "treated"))

    heatmap_data <- build_heatmap(se, top_de = top_de)
    expect_equal(attr(heatmap_data, "expression_columns"), c("control", "treated"))
    expect_true("symbol" %in% names(heatmap_data))
})
