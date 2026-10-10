#' Plan where to count plants, weeds or disease in each strip
#'
#' Lays out sampling stations for hand counts (crop establishment, weed counts,
#' disease scores) in every plot of a trial design: stations evenly spaced
#' along each plot between two end buffers, at the same positions in every
#' plot, with a few quadrats around each station kept clear of the plot sides.
#'
#' Crop counts are responses measured after the treatment, so each plot (strip)
#' is the experimental unit and the counts estimate its mean. Simulations for
#' the AAGI-CU-RD-OFE project (Milestone 5; 250 m x 12 m strips) found that:
#'
#' * spreading the stations along the whole strip matters most: putting all of
#'   them within 40 m of one end doubled the error of the strip mean, in the
#'   simulations and on 7700 real 260 m transects of commercial fields
#'   (black-grass surveys, Goodsell et al. 2023), where weed patches extend
#'   over 100--150 m;
#' * evenly spaced stations ("systematic") did slightly better than random
#'   positions, and using the same positions in every strip costs nothing;
#' * about 20 m from the strip ends and 1--2 m from its sides avoids headland
#'   effects, the build-up of the treatment at the start of a run and spill-over
#'   from the neighbouring strip;
#' * for a given number of quadrats, more stations beat more quadrats per
#'   station.
#'
#' Use [count_precision()] to choose `n_stations` and `n_quadrats`, and
#' [analyse_strip_counts()] to analyse the counts.
#'
#' @param design A trial layout from [make_trial_design()], or any data frame
#'   of lattice cells with columns `x`, `y` (cell centres, metres) and `plot`,
#'   plus optionally `treat`, `rep` and `tier`. Plots are assumed to run along
#'   `y`, as in [make_trial_design()].
#' @param n_stations Number of stations per plot.
#' @param n_quadrats Number of quadrats per station.
#' @param end_buffer Distance in metres kept clear at each end of a plot.
#' @param side_buffer Distance in metres kept clear at each side of a plot.
#' @param quadrat_spread Quadrats at a station are placed up to this distance
#'   (metres) either side of the station along the plot, and anywhere across
#'   the usable width.
#' @param layout `"systematic"` (evenly spaced, the default) or `"random"`
#'   (independent random positions between the end buffers).
#' @param same_positions For `"systematic"`, use the same station positions in
#'   every plot (`TRUE`, recommended) or a separate random start per plot.
#' @param cell_size Size of the lattice cells in `design`, in metres. Taken
#'   from `attr(design, "design")$cell_size` when available.
#' @param seed Optional integer seed.
#'
#' @return A data frame of class `ofe_stations` with one row per quadrat:
#'   `plot`, `tier`, `treat`, `rep` (where available), `station`, `quadrat`,
#'   `x`, `y` (quadrat position, metres) and `along` (distance from the start
#'   of the plot, metres). The settings and the plot outlines are kept in
#'   `attr(, "sampling")` and `attr(, "plots")`. It has a [plot()] method that
#'   draws the sampling map.
#'
#' @seealso [count_precision()], [analyse_strip_counts()],
#'   [place_point_samples()] for soil cores.
#'
#' @examples
#' d <- make_trial_design(c("Control", "Treated"), n_rep = 4,
#'                        plot_width = 12, plot_length = 250, seed = 1)
#' st <- place_count_stations(d, n_stations = 8, n_quadrats = 2, seed = 1)
#' head(st)
#' table(st$plot)
#'
#' # Same positions along every strip
#' unique(round(st$along[st$quadrat == 1], 1))
#'
#' plot(st)
#'
#' @export
place_count_stations <- function(design,
                                 n_stations = 8,
                                 n_quadrats = 2,
                                 end_buffer = 20,
                                 side_buffer = 2,
                                 quadrat_spread = 2.5,
                                 layout = c("systematic", "random"),
                                 same_positions = TRUE,
                                 cell_size = NULL,
                                 seed = NULL) {
  layout <- match.arg(layout)
  if (!is.data.frame(design) || !all(c("x", "y", "plot") %in% names(design))) {
    stop("`design` must be a data frame with columns `x`, `y` and `plot`, ",
         "such as the output of make_trial_design().", call. = FALSE)
  }
  n_stations <- as.integer(n_stations)
  n_quadrats <- as.integer(n_quadrats)
  if (is.na(n_stations) || n_stations < 1L) {
    stop("`n_stations` must be a positive integer.", call. = FALSE)
  }
  if (is.na(n_quadrats) || n_quadrats < 1L) {
    stop("`n_quadrats` must be a positive integer.", call. = FALSE)
  }
  if (end_buffer < 0 || side_buffer < 0 || quadrat_spread < 0) {
    stop("Buffers and `quadrat_spread` must be non-negative.", call. = FALSE)
  }
  if (is.null(cell_size)) cell_size <- attr(design, "design")$cell_size
  if (is.null(cell_size)) {
    ux <- sort(unique(design$x))
    cell_size <- if (length(ux) > 1L) min(diff(ux)) else 1
  }
  if (!is.null(seed)) set.seed(seed)

  cells <- design[!is.na(design$plot), , drop = FALSE]
  if ("treat" %in% names(cells)) cells <- cells[!is.na(cells$treat), , drop = FALSE]
  if (!"tier" %in% names(cells)) cells$tier <- 1L
  key <- interaction(cells$plot, cells$tier, drop = TRUE, lex.order = TRUE)
  parts <- split(seq_len(nrow(cells)), key)

  # Outline of each plot (one per plot and tier). For a make_trial_design()
  # layout use its nominal geometry: the lattice cells are snapped to
  # `cell_size` and need not match the plot edges. Otherwise use the extent
  # of the plot's cells.
  dsg <- attr(design, "design")
  ox <- min(design$x) - cell_size / 2
  oy <- min(design$y) - cell_size / 2
  plots <- do.call(rbind, lapply(parts, function(ix) {
    p <- cells$plot[ix[1]]; tr <- cells$tier[ix[1]]
    if (!is.null(dsg$plot_width) && !is.null(dsg$plot_length)) {
      per_tier <- if (identical(dsg$layout, "stack")) dsg$n_plots / 2 else dsg$n_plots
      ci <- if (tr == 1L) p else p - per_tier
      y0 <- oy + (tr - 1L) * (dsg$plot_length + dsg$gap)
      bb <- c(ox + (ci - 1) * dsg$plot_width, ox + ci * dsg$plot_width,
              y0, y0 + dsg$plot_length)
    } else {
      bb <- c(min(cells$x[ix]) - cell_size / 2, max(cells$x[ix]) + cell_size / 2,
              min(cells$y[ix]) - cell_size / 2, max(cells$y[ix]) + cell_size / 2)
    }
    data.frame(plot = p, tier = tr,
               treat = if ("treat" %in% names(cells)) as.character(cells$treat[ix[1]]) else NA,
               rep = if ("rep" %in% names(cells)) cells$rep[ix[1]] else NA,
               xmin = bb[1], xmax = bb[2], ymin = bb[3], ymax = bb[4])
  }))
  rownames(plots) <- NULL

  len <- plots$ymax - plots$ymin
  wid <- plots$xmax - plots$xmin
  if (any(len - 2 * end_buffer <= 0)) {
    stop("`end_buffer` leaves no room along the shortest plot (",
         min(len), " m long).", call. = FALSE)
  }
  if (any(wid - 2 * side_buffer <= 0)) {
    stop("`side_buffer` leaves no room across the narrowest plot (",
         min(wid), " m wide).", call. = FALSE)
  }

  u_common <- stats::runif(1)
  out <- vector("list", nrow(plots))
  for (k in seq_len(nrow(plots))) {
    usable <- len[k] - 2 * end_buffer
    along <- if (layout == "systematic") {
      u <- if (same_positions) u_common else stats::runif(1)
      end_buffer + (seq_len(n_stations) - 1 + u) * usable / n_stations
    } else {
      sort(end_buffer + stats::runif(n_stations) * usable)
    }
    a_q <- rep(along, each = n_quadrats)
    if (n_quadrats > 1L) {
      a_q <- a_q + stats::runif(length(a_q), -quadrat_spread, quadrat_spread)
    }
    a_q <- pmin(pmax(a_q, end_buffer), len[k] - end_buffer)
    x_q <- stats::runif(length(a_q), plots$xmin[k] + side_buffer,
                        plots$xmax[k] - side_buffer)
    out[[k]] <- data.frame(plot = plots$plot[k], tier = plots$tier[k],
                           treat = plots$treat[k], rep = plots$rep[k],
                           station = rep(seq_len(n_stations), each = n_quadrats),
                           quadrat = rep(seq_len(n_quadrats), n_stations),
                           x = x_q, y = plots$ymin[k] + a_q, along = a_q)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  if (all(is.na(res$treat))) res$treat <- NULL
  if (all(is.na(res$rep))) res$rep <- NULL
  attr(res, "sampling") <- list(n_stations = n_stations, n_quadrats = n_quadrats,
                                end_buffer = end_buffer, side_buffer = side_buffer,
                                quadrat_spread = quadrat_spread, layout = layout,
                                same_positions = same_positions)
  attr(res, "plots") <- plots
  attr(res, "treatments") <- if (!is.null(dsg$treatments)) dsg$treatments else
    unique(plots$treat[!is.na(plots$treat)])
  class(res) <- c("ofe_stations", "data.frame")
  res
}

#' Draw a count-sampling plan
#'
#' Draws the plots of a trial, shades the end and side buffers, and marks each
#' quadrat, so the plan can be checked by eye and taken to the field.
#'
#' @param x An `ofe_stations` object from [place_count_stations()].
#' @param col_treat Optional named vector of fill colours, one per treatment.
#' @param main Plot title.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#' @export
plot.ofe_stations <- function(x, col_treat = NULL, main = "Count-sampling plan", ...) {
  plots <- attr(x, "plots")
  smp <- attr(x, "sampling")
  trts <- attr(x, "treatments")
  if (is.null(trts)) trts <- unique(plots$treat[!is.na(plots$treat)])
  if (is.null(col_treat)) {
    pal <- c("#E8F0DC", "#DCE7F5", "#FFF3D9", "#F7DDD9", "#E6E0F0", "#DDEFEF")
    col_treat <- stats::setNames(rep(pal, length.out = max(1L, length(trts))),
                                 if (length(trts)) trts else "plot")
  }
  graphics::plot(range(c(plots$xmin, plots$xmax)), range(c(plots$ymin, plots$ymax)),
                 type = "n", asp = 1, xlab = "x (m)", ylab = "y (m)", main = main, ...)
  for (k in seq_len(nrow(plots))) {
    fill <- if (!is.na(plots$treat[k])) col_treat[[plots$treat[k]]] else col_treat[[1]]
    graphics::rect(plots$xmin[k], plots$ymin[k], plots$xmax[k], plots$ymax[k],
                   col = fill, border = "grey40")
    # end and side buffers
    graphics::rect(plots$xmin[k], plots$ymin[k], plots$xmax[k],
                   plots$ymin[k] + smp$end_buffer, col = "#00000022", border = NA)
    graphics::rect(plots$xmin[k], plots$ymax[k] - smp$end_buffer, plots$xmax[k],
                   plots$ymax[k], col = "#00000022", border = NA)
    graphics::rect(plots$xmin[k], plots$ymin[k], plots$xmin[k] + smp$side_buffer,
                   plots$ymax[k], col = "#00000014", border = NA)
    graphics::rect(plots$xmax[k] - smp$side_buffer, plots$ymin[k], plots$xmax[k],
                   plots$ymax[k], col = "#00000014", border = NA)
  }
  graphics::points(x$x, x$y, pch = 22, bg = "#414042", col = "white", cex = 0.9)
  if (length(trts) > 1L) {
    graphics::legend("topright", fill = col_treat[trts], legend = trts, bty = "n",
                     inset = c(0, -0.02), xpd = TRUE, cex = 0.8)
  }
  invisible(x)
}
