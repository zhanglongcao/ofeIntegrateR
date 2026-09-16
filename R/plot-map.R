#' Map a variable over the trial
#'
#' Draws one cell per row of `data`, coloured by `value`. It is the workhorse
#' behind the package's other plots and is exported because most of what a
#' trial team wants to look at is a map of something: the yield layer, a kriged
#' soil surface, an elevation covariate, the zones, the residuals of a fit.
#'
#' The colour scale follows the job the variable does, which is why `type`
#' exists and why `"auto"` guesses rather than always using one ramp:
#'
#' \describe{
#'   \item{`"categorical"`}{A factor or character column: identity, so distinct
#'     hues, in a fixed order.}
#'   \item{`"sequential"`}{A numeric column that is all one sign: magnitude, so
#'     one hue running light to dark.}
#'   \item{`"diverging"`}{Polarity, so two hues meeting at a neutral grey. Only
#'     for variables where zero is the reference -- residuals, differences,
#'     departures from a target. The scale is forced symmetric about zero so
#'     that equal departures either side get equal ink; an asymmetric diverging
#'     scale makes one direction look larger than it is.}
#' }
#'
#' @param data Data frame with one row per cell.
#' @param value Character; the column to colour by.
#' @param x,y Character; coordinate columns. Defaults to `"x"`/`"y"`, then
#'   `"x_centre"`/`"y_centre"`, then `"col"`/`"row"`.
#' @param type Colour scale to use. `"auto"` picks categorical for a factor and
#'   sequential otherwise; it never picks diverging, because a variable holding
#'   negative values is not thereby a signed one, and a two-hue scale on a
#'   variable with no meaningful zero invents a division that is not in the
#'   data. Ask for `"diverging"` when zero really is the reference.
#' @param main,xlab,ylab Labels. `main` defaults to the column name.
#' @param legend Logical; draw the colour key.
#' @param palette Optional character vector of colours overriding the default
#'   scale.
#' @param cell Numeric; cell size in coordinate units. Defaults to the spacing
#'   inferred from the coordinates.
#' @param na_col Colour for cells whose value is missing.
#' @param border Colour for cell borders, or `NA` for none.
#' @param overlay Optional function of no arguments, called once the cells are
#'   drawn and while the map panel is still the active coordinate system, to
#'   add labels, outlines or points. Drawing after `ofe_map()` has returned
#'   will not work: the legend key is a second panel, so the device's
#'   coordinates no longer belong to the map.
#' @param ... Passed to [graphics::plot()].
#'
#' @return The breaks and colours used, invisibly.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 24, n_col = 12, seed = 1)
#' ofe_map(sim$grid, "dense_response")
#' ofe_map(sim$grid, "treat")
#'
#' @export
ofe_map <- function(data, value, x = NULL, y = NULL,
                    type = c("auto", "sequential", "diverging", "categorical"),
                    main = NULL, xlab = NULL, ylab = NULL, legend = TRUE,
                    palette = NULL, cell = NULL, na_col = "#eeeeec",
                    border = NA, overlay = NULL, ...) {
  type <- match.arg(type)
  data <- as.data.frame(data)
  if (!value %in% names(data)) {
    stop("`value` column `", value, "` not found in `data`.", call. = FALSE)
  }
  xy <- .ofe_coords(data, x, y)
  x <- xy[1]; y <- xy[2]
  xv <- as.numeric(data[[x]])
  yv <- as.numeric(data[[y]])
  v <- data[[value]]

  # `auto` never picks diverging. A variable that happens to contain negative
  # values is not thereby a signed quantity: yield centred near zero, a
  # standardised covariate, a temperature -- none of them have a meaningful
  # midpoint, and a two-hue scale would invent one and split the paddock into
  # "red half / blue half" that means nothing. Diverging has to be asked for,
  # and the plot methods that draw residuals ask for it.
  if (type == "auto") {
    type <- if (is.factor(v) || is.character(v) || is.logical(v))
      "categorical" else "sequential"
  }

  if (is.null(cell)) cell <- .ofe_spacing(xv, yv)
  if (is.null(main)) main <- value
  if (is.null(xlab)) xlab <- x
  if (is.null(ylab)) ylab <- y

  if (type == "categorical") {
    f <- factor(v)
    cols <- if (is.null(palette)) ofe_palette("treatment", nlevels(f)) else
      rep_len(palette, nlevels(f))
    fill <- cols[as.integer(f)]
    key <- list(levels = levels(f), colours = cols)
  } else {
    vn <- as.numeric(v)
    n_col <- 64L
    if (type == "diverging") {
      lim <- max(abs(vn), na.rm = TRUE)
      zlim <- c(-lim, lim)                # symmetric: equal ink either side
      cols <- if (is.null(palette)) ofe_palette("residual", n_col) else
        grDevices::colorRampPalette(palette, space = "Lab")(n_col)
    } else {
      zlim <- range(vn, na.rm = TRUE)
      if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
      cols <- if (is.null(palette)) ofe_palette("surface", n_col) else
        grDevices::colorRampPalette(palette, space = "Lab")(n_col)
    }
    idx <- as.integer(cut(vn, seq(zlim[1], zlim[2], length.out = n_col + 1L),
                          include.lowest = TRUE))
    fill <- cols[idx]
    key <- list(zlim = zlim, colours = cols)
  }
  fill[is.na(fill)] <- na_col

  if (legend) {
    # no.readonly = TRUE: the full par() list includes read-only entries that
    # cannot be set back, and restoring it wholesale warns once per entry.
    op_layout <- graphics::par(no.readonly = TRUE)
    graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(5, 1))
    on.exit({
      graphics::layout(1)
      graphics::par(op_layout)
    }, add = TRUE)
  }
  op <- .ofe_par(mar = c(4.4, 4.4, 3.0, 0.8))
  on.exit(graphics::par(op), add = TRUE)

  # Limits at the cell edges, not the cell centres, so the outermost cells are
  # not half outside the panel.
  # asp = 1 keeps the paddock in proportion, which means the narrow axis gets
  # padded. Draw the axes ourselves, clamped to the data, so the padding does
  # not come with tick marks over ground the trial does not cover.
  xr <- range(xv) + c(-1, 1) * cell / 2
  yr <- range(yv) + c(-1, 1) * cell / 2
  graphics::plot(xv, yv, type = "n", asp = 1, main = main, xlab = xlab,
                 ylab = ylab, xlim = xr, ylim = yr, xaxs = "i", yaxs = "i",
                 axes = FALSE, ...)
  .ofe_axis(1, xr)
  .ofe_axis(2, yr)
  graphics::rect(xv - cell / 2, yv - cell / 2, xv + cell / 2, yv + cell / 2,
                 col = fill, border = border)
  if (is.function(overlay)) overlay()

  if (legend) {
    if (type == "categorical") {
      op2 <- graphics::par(mar = c(4.4, 0.4, 3.0, 0.4))
      on.exit(graphics::par(op2), add = TRUE)
      graphics::plot.new()
      graphics::legend("topleft", legend = key$levels, fill = key$colours,
                       border = NA, bty = "n", cex = 0.85,
                       text.col = .ofe_chrome[["ink2"]])
    } else {
      .ofe_colourbar(key$zlim, key$colours, label = value)
    }
  }
  key$fill <- fill
  invisible(key)
}

#' Work out which columns hold the coordinates
#' @keywords internal
#' @noRd
.ofe_coords <- function(data, x = NULL, y = NULL) {
  if (is.null(x)) {
    x <- if ("x" %in% names(data)) "x" else
      if ("x_centre" %in% names(data)) "x_centre" else
        if ("col" %in% names(data)) "col" else NULL
  }
  if (is.null(y)) {
    y <- if ("y" %in% names(data)) "y" else
      if ("y_centre" %in% names(data)) "y_centre" else
        if ("row" %in% names(data)) "row" else NULL
  }
  if (is.null(x) || is.null(y)) {
    stop("Coordinates not found: pass `x` and `y`.", call. = FALSE)
  }
  missing_cols <- setdiff(c(x, y), names(data))
  if (length(missing_cols) > 0L) {
    stop("Coordinate column(s) not found: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  c(x, y)
}

#' Cell size from the commonest gap between distinct coordinates
#' @keywords internal
#' @noRd
.ofe_spacing <- function(xv, yv) {
  gap <- function(z) {
    u <- sort(unique(z))
    if (length(u) < 2L) return(NA_real_)
    d <- diff(u)
    as.numeric(names(which.max(table(signif(d, 8)))))
  }
  s <- c(gap(xv), gap(yv))
  s <- s[is.finite(s) & s > 0]
  if (length(s) == 0L) return(1)
  min(s)
}

#' An axis drawn only over the range the data occupies
#' @keywords internal
#' @noRd
.ofe_axis <- function(side, lim) {
  at <- pretty(lim)
  at <- at[at >= lim[1] & at <= lim[2]]
  if (length(at) < 2L) at <- lim
  graphics::axis(side, at = at, col = .ofe_chrome[["axis"]],
                 col.ticks = .ofe_chrome[["axis"]],
                 col.axis = .ofe_chrome[["ink2"]], cex.axis = 0.85, las = 1)
}

#' Ink that stays legible on a given fill
#'
#' Label colour has to follow the fill it sits on: white on a pale yellow zone
#' is unreadable, and dark ink on a deep blue one equally so. Chooses by
#' relative luminance rather than by guessing per palette.
#' @keywords internal
#' @noRd
.ofe_ink_on <- function(bg) {
  lum <- function(col) {
    v <- grDevices::col2rgb(col) / 255
    lin <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
    0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
  }
  contrast <- function(a, b) (pmax(a, b) + 0.05) / (pmin(a, b) + 0.05)
  l_bg <- lum(bg)
  l_dark <- lum(.ofe_chrome[["ink"]])
  l_light <- lum(.ofe_chrome[["surface"]])
  # Pick whichever ink actually has more contrast, rather than switching at a
  # guessed luminance threshold: a mid-toned yellow is brighter than it looks
  # and takes dark ink, which a threshold near the middle gets wrong.
  out <- ifelse(contrast(l_bg, l_dark) >= contrast(l_bg, l_light),
                .ofe_chrome[["ink"]], .ofe_chrome[["surface"]])
  unname(out)
}
