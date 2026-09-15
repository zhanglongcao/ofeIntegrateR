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

#' Check that sommer is available, with an informative error
#'
#' `sommer` is only needed for the open-source joint-model engine, and it
#' compiles C++ at install time, so it is a Suggests rather than a hard
#' dependency: users who only want the kriged-covariate route should not pay
#' for a build they will never call.
#'
#' @keywords internal
#' @noRd
check_sommer <- function() {
  if (!requireNamespace("sommer", quietly = TRUE)) {
    stop(
      "`engine = \"sommer\"` requires the 'sommer' package, which is not ",
      "installed.\n",
      "Install it with install.packages(\"sommer\"), or use ",
      "`fit_integrated_kriged()` with engine = \"gls\" for an open-source ",
      "route that needs no extra dependency.",
      call. = FALSE
    )
  }
}
