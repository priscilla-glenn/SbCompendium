.onLoad <- function(libname, pkgname) {
    if (requireNamespace("reticulate", quietly = TRUE)) {
        reticulate::py_require(
            c("pandas", "xlsxwriter", "openpyxl")
        )
    }
}
