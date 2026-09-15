#' How far apart should the samples be, and how many are needed?
#'
#' Works out the grid interval at which ordinary kriging attains a target
#' precision, given a variogram, and converts that interval into a sample count
#' for a trial of known area. This is the question that has to be answered
#' before the sampling budget is set, and it cannot be answered by a rule of
#' thumb: the same sample count is generous on a paddock that varies over
#' hundreds of metres and useless on one that turns over every ten.
#'
#' Precision is expressed as the **kriging standard error relative to the field
#' standard deviation**, \eqn{\sigma_K / \sqrt{c_0 + c_1}}. A relative error of
#' 0.5 means the interpolated surface is twice as precise as guessing the field
#' mean everywhere.
#'
#' There is a floor. At an unsampled location the nugget component cannot be
#' filtered out however densely you sample, so no interval attains a relative
#' error below \eqn{\sqrt{c_0 / (c_0 + c_1)}}. A property with a nugget-to-sill
#' ratio of 0.3 cannot be mapped to better than 0.55 of the field standard
#' deviation at any cost. The returned `kse_floor` reports this, and it is
#' often the most useful number here: it says whether the target is worth
#' budgeting for at all.
#'
#' The calculation places four samples at the corners of a
#' \eqn{\Delta \times \Delta} cell and predicts its centre — the worst-supported
#' point on a regular grid — so the answer is conservative by construction.
#'
#' Variogram parameters normally come from a reconnaissance survey or from
#' ancillary data such as EM38 or several seasons of yield maps. Where no prior
#' exists, run the calculation across a plausible range and plan for the
#' shortest range you might encounter.
#'
#' @param nugget,psill,range Variogram parameters: nugget \eqn{c_0}, partial
#'   sill \eqn{c_1}, and practical range, in the units of the trial's
#'   coordinates (normally metres).
#' @param target_kse Target kriging standard error as a fraction of the field
#'   standard deviation. Must exceed `kse_floor` to be attainable.
#' @param area_ha Optional trial area in hectares. When given, the sample count
#'   implied by the interval is returned as well.
#' @param model Variogram model passed to [gstat::vgm()].
#' @param intervals Optional numeric vector of intervals to evaluate. Defaults
#'   to a sequence spanning a twentieth of the range to twice the range.
#'
#' @return A list with the largest `interval` meeting the target, the `rel_kse`
#'   achieved there, `n_samples` for `area_ha` (or `NA`), the `target_kse`, the
#'   `kse_floor` set by the nugget, and a `curve` data frame of interval against
#'   relative kriging standard error, for plotting the trade-off. `interval` is
#'   `NA` when the target lies below the floor.
#'
#' @examples
#' # A medium-range soil property over a 20 ha trial
#' s <- kriging_sample_interval(nugget = 0.2, psill = 0.8, range = 60,
#'                              target_kse = 0.6, area_ha = 20)
#' s$interval
#' s$n_samples
#'
#' # Short-range variation needs far denser sampling for the same precision
#' kriging_sample_interval(nugget = 0.2, psill = 0.8, range = 15,
#'                         target_kse = 0.6, area_ha = 20)$n_samples
#'
#' # The nugget sets a floor no sampling effort can beat
#' kriging_sample_interval(nugget = 0.5, psill = 0.5, range = 60)$kse_floor
#'
#' @export
kriging_sample_interval <- function(nugget,
                                     psill,
                                     range,
                                     target_kse = 0.5,
                                     area_ha = NULL,
                                     model = "Exp",
                                     intervals = NULL) {
  if (!is.numeric(range) || length(range) != 1L || range <= 0) {
    stop("`range` must be a single positive number.", call. = FALSE)
  }
  if (nugget < 0 || psill <= 0) {
    stop("`nugget` must be non-negative and `psill` positive.", call. = FALSE)
  }
  if (target_kse <= 0 || target_kse >= 1) {
    stop("`target_kse` must be between 0 and 1: it is the kriging standard ",
         "error as a fraction of the field standard deviation.", call. = FALSE)
  }

  if (is.null(intervals)) {
    intervals <- seq(range / 20, range * 2, length.out = 60)
  }
  intervals <- sort(unique(intervals[intervals > 0]))

  field_sd <- sqrt(nugget + psill)
  # The nugget cannot be filtered at an unsampled point, so this is the best
  # any interval can do.
  kse_floor <- sqrt(nugget) / field_sd
  vmod <- gstat::vgm(psill = psill, model = model, range = range,
                     nugget = nugget)

  # Four samples at the corners of a delta x delta cell, predicting the centre:
  # the least well supported point on a regular grid, so this is conservative.
  ok_var <- function(delta) {
    samp <- data.frame(x = c(0, delta, 0, delta),
                       y = c(0, 0, delta, delta), z = 0)
    pred <- data.frame(x = delta / 2, y = delta / 2)
    sp::coordinates(samp) <- ~ x + y
    sp::coordinates(pred) <- ~ x + y
    suppressMessages(
      gstat::krige(z ~ 1, samp, pred, model = vmod, nmax = 4,
                   debug.level = 0)$var1.var
    )
  }

  rel_kse <- sqrt(vapply(intervals, ok_var, numeric(1))) / field_sd
  curve <- data.frame(interval = intervals, rel_kse = rel_kse)

  ok <- which(rel_kse <= target_kse)
  if (!length(ok)) {
    if (target_kse <= kse_floor) {
      warning("A relative kriging standard error of ", target_kse,
              " is unattainable at any sampling density: the nugget sets a ",
              "floor of ", round(kse_floor, 3), ". Either accept a larger ",
              "target or reduce the nugget, e.g. by compositing cores.",
              call. = FALSE)
    } else {
      warning("No interval in the range tested reaches ", target_kse,
              " (floor is ", round(kse_floor, 3),
              "); widen the search with `intervals`.", call. = FALSE)
    }
    interval <- NA_real_
    achieved <- NA_real_
  } else {
    # The largest interval meeting the target is the cheapest one that does.
    interval <- intervals[max(ok)]
    achieved <- rel_kse[max(ok)]
  }

  n_samples <- if (!is.null(area_ha) && is.finite(interval)) {
    ceiling(area_ha * 10000 / interval^2)
  } else {
    NA_real_
  }

  list(interval = interval,
       rel_kse = achieved,
       n_samples = n_samples,
       target_kse = target_kse,
       kse_floor = kse_floor,
       curve = curve)
}
