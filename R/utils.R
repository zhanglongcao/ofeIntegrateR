#' Check that asreml is available, with an informative error
#'
#' `asreml` is a commercial package from VSNi and is not available on CRAN,
#' so it is listed only in Suggests. This helper gives a clear error
#' pointing collaborators to a licence/installation source rather than the
#' generic "could not find function" error.
#'
#' @keywords internal
#' @noRd
check_asreml <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop(
      "This function requires the 'asreml' package, which is commercial ",
      "software from VSNi and is not installed.\n",
      "See https://vsni.co.uk/software/asreml-r for licensing and ",
      "installation, or use `engine = \"lm\"` for an asreml-free ",
      "(non-spatial-residual) fallback where available.",
      call. = FALSE
    )
  }
}
