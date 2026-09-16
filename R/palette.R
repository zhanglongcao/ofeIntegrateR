# Colours for the package's plots.
#
# These are not a matter of taste and were not chosen by eye. Each set below was
# run through a colour-vision-deficiency validator against the light surface an
# R graphics device actually draws on, and only sets that cleared every gate are
# used. The gates that bind, and what they rule out, are recorded with each set.

#' Colours used by the package's plots
#'
#' Exposed so that a report can match its figures to the rest of its styling,
#' and so the choices can be inspected rather than taken on trust.
#'
#' Three jobs, three kinds of scale, because the job decides the scale:
#'
#' \describe{
#'   \item{`"treatment"`}{Categorical -- treatments have identity, not
#'     magnitude. The order is fixed and nested, so adding a treatment never
#'     repaints the others. Validated pairwise for up to six levels: every pair
#'     that can appear on screen together is separable both to normal colour
#'     vision and under simulated deuteranopia, protanopia and tritanopia. Past
#'     six the guarantee lapses, and the plots fall back to labelling.}
#'   \item{`"zone"`}{Ordinal -- zones from [partition_paddock()] and
#'     [partition_pseudo_env()] are ordered (by covariate mean, or along the
#'     trial), so a single hue running light to dark says so, where a
#'     categorical set would imply they are unrelated. Steps are spaced to stay
#'     distinguishable up to seven zones -- past that the fill can no longer
#'     separate them and the zone label has to carry identity.}
#'   \item{`"surface"`}{Sequential, for a continuous map -- a kriged layer, a
#'     yield surface, elevation. The same hue as `"zone"` but running from a
#'     lighter start, because a continuous scale is read as a gradient rather
#'     than as a set of levels to tell apart, and the extra range buys
#'     resolution.}
#'   \item{`"residual"`}{Diverging -- a residual has a sign and a natural zero,
#'     so two hues meet at a neutral grey. Never a rainbow: a hue at the
#'     midpoint invents a category where the data has nothing.}
#' }
#'
#' @param what Which scale: `"treatment"`, `"zone"`, `"surface"`, `"residual"`,
#'   or `"chrome"` for the ink, grid and surface colours.
#' @param n Number of colours wanted. Ignored for `"chrome"`.
#'
#' @return A character vector of hex colours; for `"chrome"`, a named vector.
#'
#' @examples
#' ofe_palette("treatment", 3)
#' ofe_palette("zone", 4)
#' ofe_palette("surface", 6)
#' ofe_palette("chrome")
#'
#' @export
ofe_palette <- function(what = c("treatment", "zone", "surface", "residual",
                                 "chrome"),
                        n = 5) {
  what <- match.arg(what)
  if (what == "chrome") return(.ofe_chrome)
  n <- max(1L, as.integer(n))
  switch(
    what,
    treatment = {
      if (n <= length(.ofe_cat)) .ofe_cat[seq_len(n)] else
        grDevices::colorRampPalette(.ofe_cat, space = "Lab")(n)
    },
    zone = if (n == 1L) .ofe_zone_ends[2] else
      grDevices::colorRampPalette(.ofe_zone_ends, space = "Lab")(n),
    surface = grDevices::colorRampPalette(.ofe_surface_ends, space = "Lab")(n),
    residual = grDevices::colorRampPalette(
      c(.ofe_div[1], .ofe_chrome[["midpoint"]], .ofe_div[2]),
      space = "Lab")(n)
  )
}

# Categorical, in fixed nested order. This is not the source palette's own
# order: that one places orange second, and orange beside yellow fails the
# normal-vision separation floor when both are on screen at once -- which is the
# usual case on a map, where every pair is adjacent. Re-ordering to put aqua
# second and drop orange entirely clears all pairs up to six.
#' @keywords internal
#' @noRd
.ofe_cat <- c("#2a78d6", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7")

# Ordinal ramp ends. The light end is held above the surface (2.06:1) so the
# first zone is a visible fill rather than a hole in the page.
#' @keywords internal
#' @noRd
.ofe_zone_ends <- c("#86b6ef", "#0d366b")

# Continuous surfaces may start lighter than the discrete ordinal ramp: the
# lightest step means "near the bottom of the scale" and is allowed to recede
# toward the page, which a discrete level is not.
#' @keywords internal
#' @noRd
.ofe_surface_ends <- c("#cde2fb", "#0d366b")

# Diverging poles: cool and warm, so the sign reads without a legend.
#' @keywords internal
#' @noRd
.ofe_div <- c("#2a78d6", "#e34948")

#' @keywords internal
#' @noRd
.ofe_chrome <- c(surface = "#fcfcfb", ink = "#0b0b0b", ink2 = "#52514e",
                 muted = "#898781", grid = "#e1e0d9", axis = "#c3c2b7",
                 midpoint = "#f0efec")

# Recessive axes and a hairline grid: the data is the figure, the frame is not.
#' @keywords internal
#' @noRd
.ofe_par <- function(...) {
  graphics::par(col.axis = .ofe_chrome[["ink2"]], col.lab = .ofe_chrome[["ink"]],
                col.main = .ofe_chrome[["ink"]], fg = .ofe_chrome[["axis"]],
                bty = "n", las = 1, cex.axis = 0.85, ...)
}

#' @keywords internal
#' @noRd
.ofe_grid <- function(h = TRUE, v = FALSE) {
  usr <- graphics::par("usr")
  if (h) graphics::abline(h = graphics::axTicks(2), col = .ofe_chrome[["grid"]],
                          lwd = 1)
  if (v) graphics::abline(v = graphics::axTicks(1), col = .ofe_chrome[["grid"]],
                          lwd = 1)
  invisible(usr)
}

# A legend strip for a continuous scale, drawn in the margin.
#' @keywords internal
#' @noRd
.ofe_colourbar <- function(zlim, cols, label = NULL, digits = 3) {
  op <- graphics::par(mar = c(5.1, 0.6, 4.1, 3.4))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i")
  y <- seq(zlim[1], zlim[2], length.out = length(cols) + 1L)
  graphics::rect(0, y[-length(y)], 1, y[-1L], col = cols, border = NA)
  graphics::axis(4, col = NA, col.ticks = .ofe_chrome[["axis"]],
                 col.axis = .ofe_chrome[["ink2"]], cex.axis = 0.8,
                 at = pretty(zlim, 4),
                 labels = format(pretty(zlim, 4), digits = digits))
  if (!is.null(label)) {
    graphics::mtext(label, side = 3, line = 0.4, cex = 0.8,
                    col = .ofe_chrome[["ink2"]], adj = 0)
  }
  invisible(NULL)
}
