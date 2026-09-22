#' Choose a cell size for gridding
#'
#' Evaluates a range of cell sizes for [grid_dense_layer()] and reports the
#' trade-off each one makes, because there is no default that survives contact
#' with a real trial: the right size depends on how fast the monitor logged, how
#' wide the plots are, and how the passes line up with the treatment boundaries.
#'
#' Two things pull in opposite directions. Small cells keep spatial detail but
#' rest on few observations each, so cell means are noisy and more cells
#' straddle a treatment boundary and have to be thrown away. Large cells are
#' well supported and mostly pure, but average away the variation the spatial
#' model exists to capture, and eventually average across the treatments
#' themselves.
#'
#' The default recommendation is the **smallest** size meeting both floors --
#' enough observations per cell, and enough cells lying inside a single
#' treatment. Smallest, because detail lost to a coarse grid cannot be
#' recovered, while the floors are what keep the cells trustworthy. It is a
#' stated rule, not a fact: the table is returned so you can apply a different
#' one.
#'
#' @param data Data frame of dense-layer observations.
#' @param response Character; the column to be aggregated.
#' @param x,y Character; projected coordinate columns, in metres.
#' @param treat Character or `NULL`; a treatment column. Without it, purity
#'   cannot be assessed and only the support columns are reported.
#' @param sizes Numeric vector of cell sizes to try. Defaults to a spread from
#'   roughly the typical spacing between observations up to twenty times it.
#' @param min_obs Numeric; the median observations per occupied cell a size
#'   must reach.
#' @param min_purity Numeric; the share of occupied cells that must lie wholly
#'   within one treatment.
#'
#' @return A data frame with one row per size: `cell_size`, `cells` in the
#'   lattice, `occupied` cells, `coverage` (the share occupied),
#'   `median_obs` per occupied cell, `pure` (share of occupied cells with
#'   `treat_purity == 1`), and `usable` (occupied, supported and pure). The
#'   recommendation and the rule behind it are attached as
#'   `attr(, "recommended")` and `attr(, "rule")`.
#'
#' @seealso [grid_dense_layer()], and [clean_yield_monitor()] to run first --
#'   deciding a cell size from a cloud that still contains pass-end artefacts
#'   sizes the grid around noise.
#'
#' @examples
#' sim <- simulate_yield_monitor(seed = 11)
#' tab <- choose_cell_size(sim$cloud, response = "yield", treat = "treat")
#' tab
#' attr(tab, "recommended")
#'
#' @export
choose_cell_size <- function(data, response, x = "x", y = "y", treat = NULL,
                             sizes = NULL, min_obs = 3, min_purity = 0.7) {
  data <- as.data.frame(data)
  needed <- c(response, x, y, treat)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(sizes)) {
    # Start near the spacing between neighbouring observations: below that,
    # cells cannot hold more than one point however the grid is drawn.
    sub <- if (nrow(data) > 400L) sample.int(nrow(data), 400L) else
      seq_len(nrow(data))
    dm <- as.matrix(stats::dist(cbind(data[[x]][sub], data[[y]][sub])))
    diag(dm) <- Inf
    step <- stats::median(apply(dm, 1, min))
    sizes <- unique(round(step * c(1, 2, 3, 5, 8, 12, 20), 1))
  }
  sizes <- sort(unique(sizes[sizes > 0]))
  if (!length(sizes)) stop("`sizes` must contain a positive value.",
                           call. = FALSE)

  rows <- lapply(sizes, function(cs) {
    g <- grid_dense_layer(data, x = x, y = y, response = response,
                          treat = treat, cell_size = cs)
    occ <- g$n_obs > 0
    pure <- if (is.null(treat)) NA_real_ else
      mean(g$treat_purity[occ] == 1, na.rm = TRUE)
    usable <- if (is.null(treat)) sum(occ & g$n_obs >= min_obs) else
      sum(occ & g$n_obs >= min_obs & !is.na(g$treat_purity) &
            g$treat_purity == 1)
    data.frame(cell_size = cs, cells = nrow(g), occupied = sum(occ),
               coverage = round(sum(occ) / nrow(g), 3),
               median_obs = stats::median(g$n_obs[occ]),
               pure = round(pure, 3), usable = usable)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  ok <- out$median_obs >= min_obs &
    (is.na(out$pure) | out$pure >= min_purity)
  attr(out, "recommended") <- if (any(ok)) out$cell_size[which(ok)[1]] else
    NA_real_
  attr(out, "rule") <- sprintf(
    "smallest size with median_obs >= %g and pure >= %g", min_obs, min_purity)
  if (!any(ok)) {
    warning("No size tried met both floors. Widen `sizes`, or lower ",
            "`min_obs` / `min_purity`.", call. = FALSE)
  }
  out
}
