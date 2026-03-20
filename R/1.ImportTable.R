#' Summarize a compendium expression table to GeneIDV3
#'
#' Takes an expression table with a TranscriptID or TranscriptIDV3 column,
#' creates GeneIDV3 by stripping the transcript suffix, and sums all other
#' numeric columns by GeneIDV3. Returns a GeneIDV3-rowname data.frame.
#'
#' @param df A data.frame or tibble expression table with TranscriptID or TranscriptIDV3.
#' @return A data.frame with GeneIDV3s as rownames and summed expression columns.
#' @export
#' @examples
#' example_df <- data.frame(
#'   TranscriptID = c("Sobic.001G000100.1", "Sobic.001G000100.2", "Sobic.001G000200.1"),
#'   sample1 = c(5, 3, 2),
#'   sample2 = c(7, 1, 4)
#' )
#' importTable(example_df)

importTable <- function(df) {

    stopifnot(is.data.frame(df))

    # Detect transcript column automatically
    transcript_col <- intersect(
        c("TranscriptID", "TranscriptIDV3"),
        names(df)
    )

    if (length(transcript_col) == 0) {
        stop("df must contain either 'TranscriptID' or 'TranscriptIDV3'.",
             call. = FALSE)
    }

    transcript_col <- transcript_col[1]  # use whichever is found

    # Create GeneIDV3 by stripping transcript suffix
    df$GeneIDV3 <- sub("\\.[0-9]+$", "", df[[transcript_col]])

    # Reorder columns: transcript column, GeneIDV3, then everything else
    df <- df[, c(transcript_col, "GeneIDV3",
                 setdiff(names(df), c(transcript_col, "GeneIDV3")))]

    # Drop transcript column (keep GeneIDV3 + expression columns)
    df <- df[, setdiff(names(df), transcript_col)]

    # Sum all columns by GeneIDV3
    out <- df |>
        dplyr::group_by(GeneIDV3) |>
        dplyr::summarise(
            dplyr::across(dplyr::everything(), sum),
            .groups = "drop"
        )

    # Move GeneIDV3 into rownames
    out <- tibble::column_to_rownames(out, var = "GeneIDV3")

    # Remove rows with any NA values
    out <- out[complete.cases(out), , drop = FALSE]

    out
}
