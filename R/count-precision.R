#' How precisely do hand counts estimate a strip mean?
#'
#' Simulates many strips with a known mean density, counts them with a given
#' design (stations evenly spaced along the strip, a few quadrats per station),
#' and reports how far the strip mean from the counts falls from the truth. The
#' result is the **95% margin**: in 95% of simulated strips the counted mean
#' was within this fraction of the true strip mean.
#'
#' The model follows the AAGI-CU-RD-OFE Milestone 5 study. On a 0.5 m grid the
#' log density is the sum of an independent cell-to-cell term (`var_cell`, per
#' 0.25 m^2), patches (`var_patch`, practical range `patch_range`) and a
#' field-scale trend (`var_trend`, practical range `trend_range`), centred so
#' that `mean_density` is the mean. Counts in a quadrat are Poisson given the
#' density. Sown crops are usually more evenly spaced than random, so the
#' Poisson assumption errs on the safe side.
#'
#' The presets were set from public establishment trials (University of
#' Adelaide, GRDC UOA1803-009RTX, CC BY 4.0): repeat counts within plots for
#' `var_cell`, and for `var_patch` the range from plot-to-plot variation in
#' research trials (`"low"`) to half the run-to-run variation of commercial
#' seeders (`"high"`). For weeds, the correlation between 20 m x 20 m
#' quadrats was calibrated to black-grass surveys of 200 commercial fields
#' (Goodsell et al. 2023, CC0): patches extend over 100--150 m. Supply your own
#' variances where a paddock has been measured.
#'
#' The margin describes one strip. Whether a *treatment difference* can be
#' detected depends mostly on the number of replicate strips and on how much
#' strips differ from each other; see [analyse_strip_counts()].
#'
#' @param n_stations Integer vector; stations per strip. Several values are
#'   evaluated on the same simulated strips.
#' @param n_quadrats Quadrats per station.
#' @param quadrat_area Quadrat area in m^2. Rounded to whole 0.25 m^2 cells
#'   (0.25, 0.5, 1 are exact).
#' @param crop Preset: `"wheat"`, `"barley"` (both the cereal preset),
#'   `"canola"`, `"weeds"`, or `"custom"` to use only the values supplied.
#' @param patchiness For the crop presets, `"low"`, `"central"` or `"high"`
#'   patch variance (see Details). Ignored for `"weeds"`.
#' @param mean_density,var_cell,var_patch,patch_range,var_trend,trend_range
#'   Model settings; any value supplied overrides the preset. Densities per
#'   m^2, variances on the log scale, ranges in metres.
#' @param strip_length,strip_width Strip size in metres.
#' @param end_buffer,side_buffer Distances kept clear at the strip ends and
#'   sides, in metres.
#' @param n_sim Number of simulated strips.
#' @param seed Optional integer seed.
#'
#' @return A data frame with one row per value of `n_stations`: the design,
#'   `counts` (stations x quadrats), `margin95` (95% margin, as a fraction),
#'   `rmse` and `bias` of the relative error. The settings used are in
#'   `attr(, "settings")`; the relative error of every simulated strip
#'   (`n_sim` rows, one column per value of `n_stations`) is in
#'   `attr(, "errors")`.
#'
#' @seealso [count_sample_size()], [place_count_stations()]
#'
#' @examples
#' # Wheat establishment, 0.25 m2 quadrats, 2 per station
#' count_precision(n_stations = c(5, 8, 12), n_quadrats = 2,
#'                 quadrat_area = 0.25, crop = "wheat", n_sim = 100, seed = 1)
#'
#' # A patchier paddock
#' count_precision(8, 2, 0.25, crop = "wheat", patchiness = "high",
#'                 n_sim = 100, seed = 1)
#'
#' @export
count_precision <- function(n_stations = c(5, 8, 12, 20),
                            n_quadrats = 2,
                            quadrat_area = 0.25,
                            crop = c("wheat", "barley", "canola", "weeds", "custom"),
                            patchiness = c("central", "low", "high"),
                            mean_density = NULL, var_cell = NULL, var_patch = NULL,
                            patch_range = NULL, var_trend = NULL, trend_range = NULL,
                            strip_length = 250, strip_width = 12,
                            end_buffer = 20, side_buffer = 2,
                            n_sim = 400, seed = NULL) {
  crop <- match.arg(crop)
  patchiness <- match.arg(patchiness)
  set_ <- count_preset(crop, patchiness)
  given <- list(mean_density = mean_density, var_cell = var_cell,
                var_patch = var_patch, patch_range = patch_range,
                var_trend = var_trend, trend_range = trend_range)
  for (nm in names(given)) if (!is.null(given[[nm]])) set_[[nm]] <- given[[nm]]
  if (crop == "custom" && any(vapply(set_, is.null, logical(1)))) {
    stop("For crop = \"custom\" supply mean_density, var_cell, var_patch, ",
         "patch_range, var_trend and trend_range.", call. = FALSE)
  }
  n_stations <- sort(unique(as.integer(n_stations)))
  if (any(is.na(n_stations)) || any(n_stations < 1L)) {
    stop("`n_stations` must be positive integers.", call. = FALSE)
  }
  n_quadrats <- as.integer(n_quadrats)
  if (is.na(n_quadrats) || n_quadrats < 1L) {
    stop("`n_quadrats` must be a positive integer.", call. = FALSE)
  }
  if (quadrat_area <= 0) stop("`quadrat_area` must be positive.", call. = FALSE)
  if (strip_length - 2 * end_buffer <= 0 || strip_width - 2 * side_buffer <= 0) {
    stop("The buffers leave no room in the strip.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  cell <- 0.5
  ny <- round(strip_length / cell); nx <- round(strip_width / cell)
  # quadrat footprint in 0.5 m cells
  n_cells <- max(1L, round(quadrat_area / cell^2))
  fc <- if (n_cells <= 2L) 1L else 2L
  fr <- as.integer(ceiling(n_cells / fc))
  area_used <- fr * fc * cell^2

  k_patch <- count_kernel(ny, nx, set_$patch_range, cell)
  # the field-scale trend is smooth: simulate it on a 2 m grid and expand
  agg <- 4L
  nyc <- ceiling(ny / agg); nxc <- ceiling(nx / agg)
  k_trend <- count_kernel(nyc, nxc, set_$trend_range, cell * agg)
  ry <- rep(seq_len(nyc), each = agg)[seq_len(ny)]
  rx <- rep(seq_len(nxc), each = agg)[seq_len(nx)]
  mu0 <- log(set_$mean_density) - (set_$var_cell + set_$var_patch + set_$var_trend) / 2

  err <- matrix(NA_real_, n_sim, length(n_stations))
  for (s in seq_len(n_sim)) {
    eta <- mu0 + sqrt(set_$var_cell) * matrix(stats::rnorm(ny * nx), ny, nx) +
      sqrt(set_$var_patch) * count_field(k_patch, ny, nx) +
      sqrt(set_$var_trend) * count_field(k_trend, nyc, nxc)[ry, rx]
    lam <- exp(eta)
    truth <- mean(lam)
    for (j in seq_along(n_stations)) {
      n <- n_stations[j]
      usable <- strip_length - 2 * end_buffer
      along <- end_buffer + (seq_len(n) - 1 + stats::runif(1)) * usable / n
      a_q <- rep(along, each = n_quadrats)
      if (n_quadrats > 1L) a_q <- a_q + stats::runif(length(a_q), -2.5, 2.5)
      a_q <- pmin(pmax(a_q, 0), strip_length - fr * cell)
      x_q <- stats::runif(length(a_q), side_buffer, strip_width - side_buffer - fc * cell)
      r0 <- floor(a_q / cell) + 1L
      c0 <- floor(x_q / cell) + 1L
      expected <- numeric(length(r0))
      for (q in seq_along(r0)) {
        expected[q] <- sum(lam[r0[q] + 0:(fr - 1L), c0[q] + 0:(fc - 1L)]) * cell^2
      }
      y <- stats::rpois(length(expected), expected)
      est <- sum(y) / (length(y) * area_used)
      err[s, j] <- est / truth - 1
    }
  }
  out <- data.frame(n_stations = n_stations, n_quadrats = n_quadrats,
                    quadrat_area = area_used,
                    counts = n_stations * n_quadrats,
                    margin95 = apply(abs(err), 2, stats::quantile, probs = 0.95, names = FALSE),
                    rmse = sqrt(colMeans(err^2)),
                    bias = colMeans(err))
  colnames(err) <- paste0("n", n_stations)
  attr(out, "errors") <- err
  attr(out, "settings") <- c(list(crop = crop, patchiness = patchiness), set_,
                             list(strip_length = strip_length, strip_width = strip_width,
                                  end_buffer = end_buffer, side_buffer = side_buffer,
                                  n_sim = n_sim))
  out
}

#' Smallest number of stations for a target precision
#'
#' Finds the fewest stations per strip whose 95% margin (from
#' [count_precision()]) is at or below `target`, for a given number of
#' quadrats per station and quadrat size.
#'
#' @param target Target 95% margin as a fraction, e.g. `0.10` for +/-10%.
#' @param max_stations Largest number of stations to consider.
#' @param ... Passed to [count_precision()] (e.g. `crop`, `quadrat_area`,
#'   `n_quadrats`, `patchiness`, `n_sim`, `seed`).
#'
#' @return A one-row data frame as from [count_precision()] for the smallest
#'   adequate design, or `NA` stations (with a warning) when even
#'   `max_stations` does not reach the target.
#'
#' @examples
#' count_sample_size(0.15, crop = "wheat", quadrat_area = 0.25, n_quadrats = 2,
#'                   n_sim = 100, seed = 1)
#'
#' @export
count_sample_size <- function(target, max_stations = 40, ...) {
  if (!is.numeric(target) || length(target) != 1L || target <= 0) {
    stop("`target` must be a single positive number, e.g. 0.10.", call. = FALSE)
  }
  tab <- count_precision(n_stations = seq_len(max_stations)[-1], ...)
  ok <- which(tab$margin95 <= target)
  if (!length(ok)) {
    warning("No design with up to ", max_stations, " stations reached a margin of ",
            target, "; consider larger or more quadrats.", call. = FALSE)
    res <- tab[nrow(tab), , drop = FALSE]
    res$n_stations <- NA_integer_
    return(res)
  }
  res <- tab[min(ok), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "settings") <- attr(tab, "settings")
  res
}

# Presets from the public establishment trials (see count_precision()).
count_preset <- function(crop, patchiness) {
  patch_cereal <- c(low = 0.003, central = 0.03, high = 0.087)
  patch_canola <- c(low = 0, central = 0.04, high = 0.087)
  switch(crop,
         wheat = , barley = list(mean_density = 120, var_cell = 0.015,
                                 var_patch = patch_cereal[[patchiness]], patch_range = 20,
                                 var_trend = 0.01, trend_range = 200),
         canola = list(mean_density = 40, var_cell = 0.02,
                       var_patch = patch_canola[[patchiness]], patch_range = 20,
                       var_trend = 0.02, trend_range = 200),
         weeds = list(mean_density = 10, var_cell = 0.6, var_patch = 0.6,
                      patch_range = 20, var_trend = 0.6, trend_range = 400),
         custom = list(mean_density = NULL, var_cell = NULL, var_patch = NULL,
                       patch_range = NULL, var_trend = NULL, trend_range = NULL))
}

# Exponential kernel for an FFT-convolution Gaussian field with unit variance.
# The practical range (correlation 0.05) of the resulting field is about 5.4
# times the kernel scale theta.
count_kernel <- function(ny, nx, range_m, cell) {
  th <- max(range_m / 5.4 / cell, 1e-6)
  pad <- ceiling(4 * th)
  NR <- stats::nextn(ny + 2 * pad); NC <- stats::nextn(nx + 2 * pad)
  dr <- pmin(0:(NR - 1), NR - 0:(NR - 1)); dc <- pmin(0:(NC - 1), NC - 0:(NC - 1))
  k <- exp(-sqrt(outer(dr^2, dc^2, "+")) / th)
  k <- k / sqrt(sum(k^2))
  list(K = stats::fft(k), NR = NR, NC = NC)
}

count_field <- function(kf, ny, nx) {
  z <- matrix(stats::rnorm(kf$NR * kf$NC), kf$NR, kf$NC)
  f <- Re(stats::fft(stats::fft(z) * kf$K, inverse = TRUE)) / (kf$NR * kf$NC)
  f[seq_len(ny), seq_len(nx)]
}
