#' Build a candidate grid for a planned trial
#'
#' Lays out a rectangular trial as a lattice of cells with treatment strips, so
#' that sampling locations can be chosen with [place_point_samples()] before any
#' yield data exist. The result has the same shape as the output of
#' [grid_dense_layer()] — `row`, `col`, `x`, `y` and a treatment column — so the
#' same code works before and after harvest.
#'
#' This lays out **strips on a rectangle**. It does not choose the number of
#' treatments, the replication, or the strip width; those come from the trial
#' design, which this package does not perform.
#'
#' @param width,height Trial extent in metres. `width` runs along `x`.
#' @param cell_size Cell size in metres. A natural choice is the harvester
#'   swath, since that is the resolution the dense layer will be gridded to.
#' @param treatments Character vector of treatment labels.
#' @param strip_width Width of each treatment strip in metres. Defaults to the
#'   extent across the strips divided by the number of treatments, i.e. one
#'   strip per treatment. Give a smaller value for replicated strips, which then
#'   cycle through the treatments.
#' @param along Axis the strips run along: `"x"` (strips stacked up the `y`
#'   axis, the usual layout for machinery passes) or `"y"`.
#'
#' @return A data frame with one row per cell: `row`, `col`, `x`, `y` (cell
#'   centres in metres) and `treat`. The layout is recorded in
#'   `attr(, "layout")`.
#'
#' @examples
#' # A 240 x 108 m trial, three treatments, replicated strips two swaths wide
#' g <- make_trial_grid(width = 240, height = 108, cell_size = 9,
#'                      treatments = c("A", "B", "C"), strip_width = 18)
#' dim(g)
#' table(g$treat)
#'
#' # Choose 30 coring locations spread across the treatments
#' set.seed(1)
#' cores <- place_point_samples(g, n = 30, design = "stratified",
#'                              coords = c("x", "y"), strata = "treat")
#' table(cores$treat)
#'
#' @export
make_trial_grid <- function(width,
                             height,
                             cell_size = 9,
                             treatments = c("A", "B", "C"),
                             strip_width = NULL,
                             along = c("x", "y")) {
  along <- match.arg(along)
  if (width <= 0 || height <= 0) {
    stop("`width` and `height` must be positive.", call. = FALSE)
  }
  if (cell_size <= 0 || cell_size > min(width, height)) {
    stop("`cell_size` must be positive and no larger than the trial.",
         call. = FALSE)
  }
  n_treat <- length(treatments)
  if (n_treat < 2L) {
    stop("`treatments` must name at least two treatments.", call. = FALSE)
  }

  n_col <- floor(width / cell_size)
  n_row <- floor(height / cell_size)
  if (n_col < 1L || n_row < 1L) {
    stop("`cell_size` leaves no whole cells in the trial.", call. = FALSE)
  }

  g <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  g$x <- (g$col - 0.5) * cell_size
  g$y <- (g$row - 0.5) * cell_size

  # Strips run along `along`, so they are indexed by the other axis.
  across <- if (along == "x") g$y else g$x
  across_extent <- if (along == "x") height else width
  if (is.null(strip_width)) strip_width <- across_extent / n_treat
  if (strip_width <= 0) stop("`strip_width` must be positive.", call. = FALSE)

  strip_idx <- floor(across / strip_width)
  g$treat <- factor(treatments[(strip_idx %% n_treat) + 1L],
                    levels = treatments)

  g <- g[order(g$row, g$col), c("row", "col", "x", "y", "treat")]
  rownames(g) <- NULL
  attr(g, "layout") <- list(width = width, height = height,
                            cell_size = cell_size, strip_width = strip_width,
                            along = along, treatments = treatments)
  g
}
