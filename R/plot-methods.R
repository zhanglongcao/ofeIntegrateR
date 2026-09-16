#' Draw the trial layout
#'
#' A plan of the plots, coloured and labelled by treatment. Look at it before
#' the trial goes in: a randomisation that happens to put the same rate at both
#' ends of a slope, or a buffer in the wrong place, is obvious on the map and
#' invisible in the design table.
#'
#' Every plot carries its treatment label as well as its colour, so the plan
#' stays readable in greyscale, to a colour-blind reader, and past the six
#' treatments the categorical palette is validated for.
#'
#' @param x An `ofe_design` from [make_trial_design()].
#' @param label Logical; write the treatment on each plot.
#' @param rep_borders Logical; outline the replicate blocks.
#' @param main Plot title.
#' @param ... Passed to [ofe_map()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' d <- make_trial_design(c("N0", "N60", "N120"), n_rep = 4, seed = 1)
#' plot(d)
#'
#' @export
plot.ofe_design <- function(x, label = TRUE, rep_borders = TRUE, main = NULL,
                            ...) {
  info <- attr(x, "design")
  if (is.null(main)) {
    main <- if (is.null(info)) "Trial design" else
      sprintf("%s layout, %d treatments x %d reps  (%.1f ha)", info$layout,
              length(info$treatments), info$n_rep, info$area_ha)
  }
  d <- as.data.frame(x)
  cols <- ofe_palette("treatment", nlevels(factor(d$treat)))
  cell <- .ofe_spacing(d$x, d$y)

  ofe_map(d, "treat", type = "categorical", main = main, overlay = function() {
    if (isTRUE(rep_borders) && "rep" %in% names(d)) {
      for (r in unique(stats::na.omit(d$rep))) {
        i <- which(d$rep == r)
        graphics::rect(min(d$x[i]) - cell / 2, min(d$y[i]) - cell / 2,
                       max(d$x[i]) + cell / 2, max(d$y[i]) + cell / 2,
                       border = .ofe_chrome[["ink"]], lwd = 1.2)
      }
    }
    if (isTRUE(label) && "plot" %in% names(d)) {
      for (p in unique(stats::na.omit(d$plot))) {
        i <- which(d$plot == p)
        tr <- factor(d$treat)[i][1]
        graphics::text(mean(range(d$x[i])), mean(range(d$y[i])),
                       labels = as.character(tr), cex = 0.7, font = 2,
                       srt = if (diff(range(d$y[i])) > diff(range(d$x[i])))
                         90 else 0,
                       col = .ofe_ink_on(cols[as.integer(tr)]))
      }
    }
  }, ...)
  invisible(x)
}

#' Draw a zone map
#'
#' The zones from [partition_paddock()] or [partition_pseudo_env()], as a map.
#' This is the check that matters on a partition: whether the zones are shapes
#' you could actually work, and whether the boundaries fall where the paddock
#' changes or somewhere arbitrary.
#'
#' Zones are ordered, so they are drawn on one hue running light to dark rather
#' than in unrelated colours, and each carries its number, so identity never
#' rests on the fill -- which matters past seven zones, where the steps stop
#' being separable.
#'
#' @param x An `ofe_zones` object from [partition_paddock()] or
#'   [partition_pseudo_env()].
#' @param zone Character; the zone column, if it was not named `"zone"`.
#' @param label Logical; write the zone number on each zone.
#' @param main Plot title.
#' @param ... Passed to [ofe_map()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' g <- expand.grid(row = 1:24, col = 1:18)
#' g$x <- g$col * 10
#' g$y <- g$row * 10
#' set.seed(1)
#' g$elevation <- 100 + 6 * exp(-((g$y - 120)^2) / 2000) + rnorm(nrow(g), 0, .2)
#' plot(partition_paddock(g, covariates = "elevation"))
#'
#' @export
plot.ofe_zones <- function(x, zone = "zone", label = TRUE, main = NULL, ...) {
  info <- attr(x, "partition")
  d <- as.data.frame(x)
  if (!zone %in% names(d)) {
    stop("Zone column `", zone, "` not found. Pass `zone = `.", call. = FALSE)
  }
  k <- nlevels(factor(d[[zone]]))
  if (is.null(main)) {
    main <- if (is.null(info)) paste(k, "zones") else
      sprintf("%s: %d zone%s from %s", info$method, info$n_zones %||% info$k,
              if ((info$n_zones %||% info$k) == 1L) "" else "s",
              paste(info$covariates %||% info$response, collapse = " + "))
  }
  cols <- ofe_palette("zone", k)
  xy <- .ofe_coords(d)
  ofe_map(d, zone, type = "categorical", main = main, palette = cols,
          overlay = function() {
            if (!isTRUE(label)) return(invisible(NULL))
            lvs <- levels(factor(d[[zone]]))
            for (j in seq_along(lvs)) {
              i <- which(as.character(d[[zone]]) == lvs[j])
              if (length(i) == 0L) next
              # Median rather than mean: on an L-shaped zone the mean can land
              # outside the zone entirely, labelling a neighbour.
              graphics::text(stats::median(d[[xy[1]]][i]),
                             stats::median(d[[xy[2]]][i]),
                             labels = lvs[j], cex = 1.1, font = 2,
                             col = .ofe_ink_on(cols[j]))
            }
          }, ...)
  invisible(x)
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Diagnostic plots for a fitted spatial model
#'
#' Four panels, and the first is the one that earns its place. A spatial model
#' is fitted precisely to absorb the field's pattern, so the question it has to
#' answer is whether any pattern is left: if the residual map still shows
#' patches, the residual structure has not done its job and the treatment
#' standard errors are optimistic. The residual variogram is the same question
#' as a curve -- a flat one is what you want, and one still climbing at short
#' lags says the correlation was not absorbed.
#'
#' The residual map uses a diverging scale forced symmetric about zero, so that
#' over- and under-prediction of the same size get the same weight of ink.
#'
#' @param x An `ofe_fit` from [fit_ofe()].
#' @param which Integer vector selecting panels: 1 residual map, 2 residuals
#'   against fitted values, 3 normal quantile plot, 4 residual variogram.
#' @param ... Ignored.
#'
#' @return A list with the residual variogram (when drawn), invisibly.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 1)
#' fit <- fit_ofe(dense_response ~ treat, data = sim$grid)
#' plot(fit)
#'
#' @export
plot.ofe_fit <- function(x, which = 1:4, ...) {
  d <- x$data
  r <- x$residuals
  if (length(r) != nrow(d)) {
    stop("The fit and its data are out of step; refit before plotting.",
         call. = FALSE)
  }
  pos <- intersect(all.vars(x$residual), names(d))
  pos <- pos[vapply(d[pos], function(z)
    is.numeric(z) || is.factor(z), logical(1))]
  coords <- if (length(pos) >= 2L) pos[1:2] else NULL
  # A lattice `row` runs up the paddock and `col` across it, so drawn as-given
  # the map comes out transposed. Put col on the horizontal axis whichever
  # order the residual formula happened to name them in.
  if (!is.null(coords) && all(c("row", "col") %in% coords)) {
    coords <- c("col", "row")
  }

  which <- intersect(which, 1:4)
  n_panel <- length(which)
  if (n_panel == 0L) return(invisible(NULL))
  op_layout <- graphics::par(no.readonly = TRUE)
  on.exit({graphics::layout(1); graphics::par(op_layout)}, add = TRUE)
  if (n_panel > 1L) {
    graphics::par(mfrow = c(ceiling(n_panel / 2), min(2L, n_panel)))
  }
  op <- .ofe_par(mar = c(4.3, 4.3, 3.0, 1.2))
  on.exit(graphics::par(op), add = TRUE)

  vg <- NULL
  dd <- d
  dd$.resid <- r

  if (1L %in% which) {
    if (is.null(coords)) {
      graphics::plot.new()
      graphics::title(main = "Residual map unavailable")
      graphics::mtext("no position columns in the residual formula", side = 1,
                      cex = 0.8, col = .ofe_chrome[["ink2"]])
    } else {
      # legend = FALSE: a colour key inside a 2x2 diagnostic panel steals more
      # room than it explains, and the sign is what is being read here.
      ofe_map(dd, ".resid", x = coords[1], y = coords[2], type = "diverging",
              main = "Residuals in space", legend = FALSE,
              xlab = coords[1], ylab = coords[2])
    }
  }

  if (2L %in% which) {
    graphics::plot(x$fitted, r, pch = 16, cex = 0.5,
                   col = grDevices::adjustcolor(.ofe_chrome[["ink2"]], 0.45),
                   main = "Residuals vs fitted", xlab = "fitted",
                   ylab = "residual")
    .ofe_grid(h = TRUE)
    graphics::abline(h = 0, col = .ofe_chrome[["axis"]], lty = 2)
    ok <- stats::complete.cases(x$fitted, r)
    if (sum(ok) > 10L) {
      lo <- stats::lowess(x$fitted[ok], r[ok])
      graphics::lines(lo, col = .ofe_cat[1], lwd = 2)
    }
  }

  if (3L %in% which) {
    stats::qqnorm(r, pch = 16, cex = 0.5,
                  col = grDevices::adjustcolor(.ofe_chrome[["ink2"]], 0.45),
                  main = "Normal Q-Q", xlab = "theoretical", ylab = "residual")
    .ofe_grid(h = TRUE)
    stats::qqline(r, col = .ofe_cat[1], lwd = 2)
  }

  if (4L %in% which) {
    vg <- tryCatch({
      tmp <- data.frame(.x = as.numeric(as.character(dd[[coords[1]]])),
                        .y = as.numeric(as.character(dd[[coords[2]]])),
                        .r = r)
      tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE]
      if (nrow(tmp) > 1500L) tmp <- tmp[sample.int(nrow(tmp), 1500L), ]
      ofe_variogram(tmp, value = ".r", x = ".x", y = ".y",
                    model = "exponential")
    }, error = function(e) NULL)
    if (is.null(vg)) {
      graphics::plot.new()
      graphics::title(main = "Residual variogram unavailable")
    } else {
      # The verdict goes in the axis label rather than a separate mtext: an
      # mtext added after plot() has returned is placed against restored
      # margins and lands on top of the axis title.
      plot(vg, main = "Residual variogram",
           xlab = if (vg$nugget_ratio > 0.75)
             "distance  (flat: the correlation was absorbed)" else
               "distance  (still rising: correlation left in the residuals)")
    }
  }
  invisible(list(variogram = vg))
}
