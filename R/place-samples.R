#' Choose where to put point samples
#'
#' Picks `n` locations from a trial grid under a named sampling design. The
#' number of samples fixes the cost, so this function answers the complementary
#' question: given that budget, where should the cores go?
#'
#' Simulation work for the AAGI-CU-RD-OFE project found that spreading samples
#' over the trial beats both simple random placement and clustered placement
#' when the samples are there to build a covariate surface, though the effect is
#' second-order next to the sample count. Clustered `"nested"` designs
#' concentrate pairs at short lags, which is what identifies a variogram range;
#' they are a reconnaissance tool and are a poor choice for interpolation.
#'
#' @param grid Data frame of candidate locations, typically the output of
#'   [grid_dense_layer()] or the `grid` element of [simulate_ofe_trial()].
#' @param n Integer; number of samples to place.
#' @param design One of:
#'   \describe{
#'     \item{`"random"`}{simple random sample of grid cells.}
#'     \item{`"grid"`}{systematic lattice spanning the trial, spaced to respect
#'       its aspect ratio.}
#'     \item{`"stratified"`}{one random cell per stratum, strata formed by
#'       crossing `strata` (if given) with blocks along the long axis. Spatially
#'       balanced and executable without a lattice survey.}
#'     \item{`"nested"`}{clustered stages at geometrically increasing spacings
#'       around a few centres, for variogram reconnaissance rather than
#'       interpolation.}
#'   }
#' @param coords Character vector of length 2 naming the coordinate columns of
#'   `grid` (defaults to the `col`/`row` lattice indices this package uses).
#' @param strata Optional character; a column of `grid` (e.g. the treatment)
#'   whose levels should each be represented under `design = "stratified"`.
#' @param spacings Numeric vector of stage spacings for `design = "nested"`, in
#'   the units of `coords`. Defaults to `NULL`, which derives three stages from
#'   the trial extent so the design works at any scale; supply your own when the
#'   lags of interest are known from a prior variogram.
#'
#' @return The selected rows of `grid`, with the design recorded in
#'   `attr(, "design")`.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 5, seed = 1)
#'
#' set.seed(1)
#' pts_random <- place_point_samples(sim$grid, n = 12, design = "random")
#' pts_grid   <- place_point_samples(sim$grid, n = 12, design = "grid")
#' pts_strat  <- place_point_samples(sim$grid, n = 12, design = "stratified",
#'                                   strata = "treat")
#' nrow(pts_grid)
#' table(pts_strat$treat)
#'
#' @export
place_point_samples <- function(grid,
                                 n,
                                 design = c("random", "grid", "stratified", "nested"),
                                 coords = c("col", "row"),
                                 strata = NULL,
                                 spacings = NULL) {
  design <- match.arg(design)
  if (!all(coords %in% names(grid))) {
    stop("`coords` columns must be present in `grid`: ",
         paste(setdiff(coords, names(grid)), collapse = ", "), call. = FALSE)
  }
  if (!is.null(strata) && !strata %in% names(grid)) {
    stop("`strata` column '", strata, "' not found in `grid`.", call. = FALSE)
  }
  n <- as.integer(n)
  if (is.na(n) || n < 1L) stop("`n` must be a positive integer.", call. = FALSE)
  if (n > nrow(grid)) {
    stop("`n` (", n, ") exceeds the number of candidate locations (",
         nrow(grid), ").", call. = FALSE)
  }

  cx <- grid[[coords[1]]]
  cy <- grid[[coords[2]]]

  # Top an under-filled selection up at random so every design returns n rows.
  top_up <- function(sel) {
    sel <- unique(sel)
    if (length(sel) > n) sel <- sample(sel, n)
    if (length(sel) < n) {
      pool <- setdiff(seq_len(nrow(grid)), sel)
      sel <- c(sel, sample(pool, n - length(sel)))
    }
    sel
  }

  idx <- switch(
    design,

    random = sample(nrow(grid), n),

    grid = {
      # Split n between the axes in proportion to the trial's extent.
      span_x <- diff(range(cx)); span_y <- diff(range(cy))
      if (span_x <= 0 || span_y <= 0) {
        top_up(sample(nrow(grid), n))
      } else {
        n_x <- max(2L, round(sqrt(n * span_x / span_y)))
        n_y <- max(2L, ceiling(n / n_x))
        tx <- seq(min(cx), max(cx), length.out = n_x)
        ty <- seq(min(cy), max(cy), length.out = n_y)
        want <- expand.grid(x = tx, y = ty)
        # Snap each lattice node to its nearest candidate location.
        sel <- vapply(seq_len(nrow(want)), function(i) {
          which.min((cx - want$x[i])^2 + (cy - want$y[i])^2)
        }, integer(1))
        top_up(sel)
      }
    },

    stratified = {
      n_blk <- max(1L, floor(n / max(1L, if (is.null(strata)) 1L
                                     else length(unique(grid[[strata]])))))
      blk <- cut(cy, breaks = max(1L, n_blk), labels = FALSE,
                 include.lowest = TRUE)
      key <- if (is.null(strata)) blk else interaction(grid[[strata]], blk,
                                                        drop = TRUE)
      parts <- split(seq_len(nrow(grid)), key)
      sel <- unname(vapply(parts, function(ix) ix[sample(length(ix), 1L)],
                           integer(1)))
      top_up(sel)
    },

    nested = {
      span_x <- diff(range(cx)); span_y <- diff(range(cy))
      extent <- max(span_x, span_y)
      if (is.null(spacings)) {
        # Three stages spanning short to medium lags, scaled to the trial so
        # the design behaves the same on a lattice or in metres.
        spacings <- extent * c(0.05, 0.12, 0.28)
      }
      n_stage <- length(spacings)
      n_cent <- max(2L, ceiling(n / (n_stage + 1L)))

      # Centres sit inside the trial, not on its boundary, so that satellites
      # at every stage have somewhere to go and genuinely cluster.
      margin <- max(spacings) * 0.6
      k_x <- max(1L, round(sqrt(n_cent * max(span_x, 1) / max(span_y, 1))))
      k_y <- max(1L, ceiling(n_cent / k_x))
      seq_in <- function(lo, hi, k, m) {
        lo2 <- min(lo + m, (lo + hi) / 2); hi2 <- max(hi - m, (lo + hi) / 2)
        if (k == 1L) (lo2 + hi2) / 2 else seq(lo2, hi2, length.out = k)
      }
      centres <- expand.grid(x = seq_in(min(cx), max(cx), k_x, margin),
                             y = seq_in(min(cy), max(cy), k_y, margin))
      centres <- centres[seq_len(min(n_cent, nrow(centres))), , drop = FALSE]

      nearest <- function(tx, ty) which.min((cx - tx)^2 + (cy - ty)^2)
      sel <- vapply(seq_len(nrow(centres)),
                    function(i) nearest(centres$x[i], centres$y[i]), integer(1))
      stage <- 1L
      while (length(unique(sel)) < n && stage <= n_stage * 4L) {
        d <- spacings[((stage - 1L) %% n_stage) + 1L]
        ang <- stats::runif(nrow(centres), 0, 2 * pi)
        sat <- vapply(seq_len(nrow(centres)), function(i) {
          nearest(centres$x[i] + d * cos(ang[i]),
                  centres$y[i] + d * sin(ang[i]))
        }, integer(1))
        sel <- c(sel, sat)
        stage <- stage + 1L
      }
      top_up(sel)
    }
  )

  out <- grid[idx, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "design") <- design
  out
}
