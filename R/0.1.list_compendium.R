#' List or search datasets in the sorghum compendium
#'
#' Prints available datasets in `sorghum_compendium` and optionally
#' filters them by a search pattern.
#'
#' @param compendium A compendium list (e.g., `sorghum_compendium`).
#' @param pattern Optional character string used to filter dataset names
#'   using `grep()`. Case-insensitive.
#' @param return Logical. If TRUE, returns the matching dataset names.
#'
#' @return Invisibly returns the matching dataset names.
#' @export
#' @examples
#' data(sorghum_compendium)
#'
#' # List all datasets
#' list_compendium(sorghum_compendium)
#'
#' # Search for diel datasets
#' list_compendium(sorghum_compendium, pattern = "diel")
#'
#' # Return the matching names
#' matches <- list_compendium(sorghum_compendium, "100M", return = TRUE)
#'
#' # Create a small example compendium
#' example_compendium <- list(
#'   diel_100M = data.frame(x = 1:3),
#'   diel_430  = data.frame(x = 4:6),
#'   drought_100M = data.frame(x = 7:9)
#' )
#'
#' # List all datasets
#' list_compendium(example_compendium)
#'
#' # Search for diel datasets
#' list_compendium(example_compendium, pattern = "diel")
#'
#' # Return the matching dataset names
#' matches <- list_compendium(example_compendium, "100M", return = TRUE)
#' matches
#'
list_compendium <- function(compendium, pattern = NULL, return = FALSE) {

    if (!is.list(compendium)) {
        stop("compendium must be a list.", call. = FALSE)
    }

    dataset_names <- names(compendium)

    if (is.null(dataset_names)) {
        stop("Compendium does not contain named datasets.", call. = FALSE)
    }

    if (!is.null(pattern)) {
        dataset_names <- grep(pattern, dataset_names,
                              value = TRUE, ignore.case = TRUE)
    }

    if (length(dataset_names) == 0) {
        message("No datasets matched.")
        return(invisible(NULL))
    }

    cat(paste(dataset_names, collapse = "\n"), "\n")

    if (return) {
        return(dataset_names)
    }

    invisible(dataset_names)
}
