#' Aggregate an irregular dense layer onto a regular lattice
#'
#' Yield-monitor and proximal-sensor data arrive as an irregular cloud of
#' GPS-referenced observations along harvester or sprayer passes, not as the
#' rectangular lattice that a separable AR1\eqn{\times}AR1 residual structure
#' requires. This function is the first step of the integration pipeline: it
#' snaps the dense layer onto a regular grid of a chosen cell size, aggregates
#' the response within each cell, and returns a **complete** lattice with
#' integer `row` and `col` indices ready for [fit_integrated_kriged()] or
#' [fit_integrated_joint()].
#'
#' Cells containing no observations are retained with an `NA` response when
#' `keep_empty = TRUE` (the default). This matters: `asreml`'s `ar1():ar1()`
#' residual is defined over the full row-by-column lattice, so gaps must be
#' present as missing values rather than dropped rows.
#'
#' @param data Data frame of dense-layer observations.
#' @param x,y Character; names of the projected coordinate columns in `data`
#'   (metres, e.g. UTM/MGA easting and northing). Longitude/latitude in degrees
#'   will produce meaningless cell sizes -- project first.
#' @param response Character; name of the column to aggregate (e.g. yield).
#' @param treat Character or `NULL`; name of a treatment column. When supplied,
#'   each cell is assigned the treatment of the majority of its observations,
#'   and `treat_purity` reports the majority share so that cells straddling a
#'   strip boundary can be identified and removed.
#' @param cell_size Numeric; the length of a (square) grid cell in the units of
#'   `x` and `y`. Choose it from the trial geometry: small enough to preserve
#'   the spatial pattern, large enough that most cells hold several
#'   observations.
#' @param fun Function used to aggregate `response` within a cell
#'   (default [mean()]); it is called with `na.rm = TRUE`.
#' @param n_min Integer; cells with fewer than `n_min` observations have their
#'   response set to `NA`. Use this to discard cells supported by one or two
#'   noisy yield-monitor pings.
#' @param keep_empty Logical; return the complete lattice including cells with
#'   no observations (default `TRUE`).
#'
#' @return A data frame, one row per lattice cell, with columns `row`, `col`
#'   (integer indices starting at 1), `x_centre`, `y_centre`, the aggregated
#'   `response`, `n_obs`, and -- when `treat` is given -- the majority `treat`
#'   and its `treat_purity`. The grid origin and cell size are attached as
#'   `attr(, "grid")`.
#'
#' @examples
#' # An irregular cloud of passes across a 100 x 60 m trial
#' set.seed(1)
#' pass_y <- rep(seq(2, 58, by = 4), each = 120)
#' cloud <- data.frame(
#'   x = rep(seq(1, 99, length.out = 120), times = length(unique(pass_y))),
#'   y = pass_y + stats::rnorm(length(pass_y), 0, 0.3)
#' )
#' cloud$treat <- cut(cloud$x, breaks = c(0, 33, 66, 100), labels = c("A", "B", "C"))
#' cloud$yield <- as.numeric(cloud$treat) + stats::rnorm(nrow(cloud), 0, 0.5)
#'
#' g <- grid_dense_layer(cloud, response = "yield", treat = "treat", cell_size = 5)
#' head(g)
#' table(is.na(g$yield))
#'
#' @export
grid_dense_layer <- function(data,
                              x = "x",
                              y = "y",
                              response,
                              treat = NULL,
                              cell_size = 5,
                              fun = mean,
                              n_min = 1L,
                              keep_empty = TRUE) {
  need <- c(x, y, response, treat)
  missing_cols <- setdiff(need, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(cell_size) || length(cell_size) != 1L || cell_size <= 0) {
    stop("`cell_size` must be a single positive number.", call. = FALSE)
  }

  xs <- data[[x]]
  ys <- data[[y]]
  ok <- is.finite(xs) & is.finite(ys)
  if (!any(ok)) stop("No observations with finite coordinates.", call. = FALSE)
  data <- data[ok, , drop = FALSE]
  xs <- xs[ok]; ys <- ys[ok]

  # Lattice indices, counting from the south-west corner of the data extent.
  x0 <- min(xs); y0 <- min(ys)
  col_i <- as.integer(floor((xs - x0) / cell_size)) + 1L
  row_i <- as.integer(floor((ys - y0) / cell_size)) + 1L
  n_col <- max(col_i); n_row <- max(row_i)

  key <- (row_i - 1L) * n_col + col_i

  # Aggregate the response and count the support of each cell.
  agg_val <- tapply(data[[response]], key, function(v) fun(v, na.rm = TRUE))
  agg_n   <- tapply(data[[response]], key, function(v) sum(!is.na(v)))

  out <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  out_key <- (out$row - 1L) * n_col + out$col

  idx <- match(as.character(out_key), names(agg_val))
  out[[response]] <- unname(agg_val[idx])
  out$n_obs <- unname(agg_n[idx])
  out$n_obs[is.na(out$n_obs)] <- 0L

  # Cells with too little support are unreliable; blank the response but keep
  # the cell so the lattice stays complete.
  out[[response]][out$n_obs < n_min] <- NA_real_

  out$x_centre <- x0 + (out$col - 0.5) * cell_size
  out$y_centre <- y0 + (out$row - 0.5) * cell_size

  if (!is.null(treat)) {
    tv <- as.character(data[[treat]])
    modal <- tapply(tv, key, function(v) {
      v <- v[!is.na(v)]
      if (!length(v)) return(NA_character_)
      tb <- table(v)
      names(tb)[which.max(tb)]
    })
    purity <- tapply(tv, key, function(v) {
      v <- v[!is.na(v)]
      if (!length(v)) return(NA_real_)
      max(table(v)) / length(v)
    })
    out[[treat]] <- factor(unname(modal[idx]),
                            levels = sort(unique(stats::na.omit(tv))))
    out$treat_purity <- unname(purity[idx])
  }

  if (!keep_empty) out <- out[out$n_obs > 0L, , drop = FALSE]

  ord <- order(out$row, out$col)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  front <- c("row", "col", "x_centre", "y_centre")
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]

  attr(out, "grid") <- list(x0 = x0, y0 = y0, cell_size = cell_size,
                             n_row = n_row, n_col = n_col)
  out
}
