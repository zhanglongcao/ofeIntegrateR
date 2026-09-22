#' Simulate a trial as a yield monitor actually records it
#'
#' [simulate_ofe_trial()] produces a tidy lattice: one observation per cell,
#' treatments aligned to whole columns. Real dense layers look nothing like
#' that, and code that only ever sees the tidy version tends to break on first
#' contact with a harvester file. This function generates the awkward version
#' --- GPS-referenced points along machinery passes, position error, a paddock
#' that is not a rectangle, occasional missing passes, and point samples taken
#' wherever the sampler could reach --- so a workflow can be tested end to end
#' through [grid_dense_layer()].
#'
#' Treatment strips run along the direction of travel, as they do in practice,
#' each covering `band_swaths` adjacent passes.
#'
#' @param field_x,field_y Paddock extent, in metres.
#' @param swath Harvester swath width (m); passes run along the `x` axis.
#' @param along_spacing Distance between recorded points within a pass (m).
#' @param band_swaths Number of adjacent swaths per treatment strip.
#' @param treat_effects Named or unnamed numeric vector of true treatment
#'   effects; treatments are labelled `A`, `B`, ... in order.
#' @param point_range,point_psill,point_nugget Exponential variogram parameters
#'   for the point-source surface, in metres.
#' @param dense_var_weight Coefficient relating the point-source surface to the
#'   dense response.
#' @param noise_sd Standard deviation of measurement noise on the dense response.
#' @param gps_jitter Standard deviation of across-track position error (m).
#' @param skip_pass_prob Probability that a pass is missing entirely.
#' @param cut_corner Clip a triangle off one corner, so the paddock is not a
#'   rectangle.
#' @param n_point_samples Number of point-source samples to draw.
#' @param point_obs_noise_sd Measurement noise on the point samples.
#' @param point_design `"random"` spreads samples over the trial;
#'   `"strip_ends"` places them at both ends of every strip, the pattern growers
#'   commonly use when each sample is expensive to collect.
#' @param defects Logical or a named list; add the faults a real monitor file
#'   carries, and label each damaged point in a `defect` column so a cleaning
#'   rule can be scored rather than eyeballed. `FALSE` (default) leaves the
#'   cloud sound. `TRUE` uses settings in the range reported for commercial
#'   monitors: the opening and closing 12 m of each pass read low and high as
#'   grain reaches the sensor (`ramp_m`, `ramp_low`, `ramp_high`), the machine
#'   slows into each turn (`slow_m`, `slow_factor`), and 2% of readings are
#'   dropouts or spikes (`outlier_prob`, `spike_factor`). A named list
#'   overrides any of those.
#' @param truth_res Resolution (m) of the grid on which the point-source surface
#'   is simulated before being read off at the recorded locations.
#' @param seed Optional integer seed.
#'
#' @return A list with `cloud` (the irregular dense observations, with `x`, `y`,
#'   `treat`, `point_true` and `yield`), `point_samples` (sparse observations at
#'   arbitrary coordinates) and `true_effects`.
#'
#' @examples
#' sim <- simulate_yield_monitor(n_point_samples = 20, seed = 2)
#' nrow(sim$cloud)
#'
#' # The whole point: this needs gridding before it can be modelled
#' g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
#'                       cell_size = 9)
#' table(empty = g$n_obs == 0)
#'
#' @export
simulate_yield_monitor <- function(field_x = 240,
                                    field_y = 108,
                                    swath = 9,
                                    along_spacing = 1.5,
                                    band_swaths = 2L,
                                    treat_effects = c(0, 0.8, 1.6),
                                    point_range = 25,
                                    point_psill = 1.2,
                                    point_nugget = 0.1,
                                    dense_var_weight = 0.6,
                                    noise_sd = 0.4,
                                    gps_jitter = 0.5,
                                    skip_pass_prob = 0.08,
                                    cut_corner = TRUE,
                                    n_point_samples = 30L,
                                    point_obs_noise_sd = 0.3,
                                    point_design = c("random", "strip_ends"),
                                    defects = FALSE,
                                    truth_res = 3,
                                    seed = NULL) {
  point_design <- match.arg(point_design)
  if (!is.null(seed)) set.seed(seed)

  n_treat <- length(treat_effects)
  names(treat_effects) <- LETTERS[seq_len(n_treat)]

  # Point-source surface on a coarse truth grid, read off at any location.
  tg <- expand.grid(x = seq(truth_res / 2, field_x, by = truth_res),
                    y = seq(truth_res / 2, field_y, by = truth_res))
  vg <- gstat::vgm(psill = point_psill, model = "Exp", range = point_range,
                   nugget = point_nugget)
  gsim <- gstat::gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE,
                       beta = 0, model = vg, nmax = 30)
  tg$point_true <- stats::predict(gsim, newdata = tg, nsim = 1,
                                  debug.level = 0)$sim1
  nx <- length(unique(tg$x))
  read_truth <- function(x, y) {
    ix <- pmin(pmax(round((x - truth_res / 2) / truth_res) + 1, 1), nx)
    iy <- pmin(pmax(round((y - truth_res / 2) / truth_res) + 1, 1),
               length(unique(tg$y)))
    tg$point_true[(iy - 1) * nx + ix]
  }

  # Passes along x, some missing.
  pass_centres <- seq(swath / 2, field_y, by = swath)
  pass_centres <- pass_centres[stats::runif(length(pass_centres)) > skip_pass_prob]
  if (!length(pass_centres)) {
    stop("All passes were skipped; lower `skip_pass_prob`.", call. = FALSE)
  }
  along <- seq(0, field_x, by = along_spacing)
  # A harvester drives alternate passes in opposite directions, and the file
  # arrives in recording order. Both matter to anything that has to work out
  # where a pass starts and ends, so they are simulated rather than assumed.
  cloud <- do.call(rbind, lapply(seq_along(pass_centres), function(i) {
    yc <- pass_centres[i]
    xs <- if (i %% 2L == 0L) rev(along) else along
    data.frame(pass = i, x = xs,
               y = yc + stats::rnorm(length(xs), 0, gps_jitter))
  }))
  cloud$speed <- stats::rnorm(nrow(cloud), 1.6, 0.08)

  band_width <- swath * band_swaths
  band_idx <- floor(cloud$y / band_width)
  cloud$treat <- factor(LETTERS[(band_idx %% n_treat) + 1L],
                        levels = LETTERS[seq_len(n_treat)])

  if (cut_corner) {
    x_c <- 0.8 * field_x; y_c <- 0.8 * field_y
    y_line <- field_y - (field_y - y_c) * (cloud$x - x_c) / (field_x - x_c)
    cloud <- cloud[!(cloud$x > x_c & cloud$y > y_line), , drop = FALSE]
  }
  cloud <- cloud[cloud$y >= 0 & cloud$y <= field_y, , drop = FALSE]

  cloud$point_true <- read_truth(cloud$x, cloud$y)
  cloud$yield <- treat_effects[as.character(cloud$treat)] +
    dense_var_weight * cloud$point_true +
    stats::rnorm(nrow(cloud), 0, noise_sd)

  # Every point is sound until a defect is applied, and each defect records
  # itself. Without that label a cleaning routine can only be eyeballed; with
  # it, how much of the damage a rule actually finds can be measured.
  cloud$defect <- NA_character_
  if (!isFALSE(defects)) {
    spec <- .ym_defect_spec(defects)
    cloud <- cloud[order(cloud$pass, seq_len(nrow(cloud))), , drop = FALSE]

    for (p in unique(cloud$pass)) {
      i <- which(cloud$pass == p)
      if (length(i) < 6L) next
      d_along <- abs(cloud$x[i] - cloud$x[i][1])
      # Grain takes a few seconds to reach the sensor, so the opening metres of
      # a pass read low and the closing metres read high.
      lead <- i[d_along <= spec$ramp_m]
      tail_i <- i[d_along >= max(d_along) - spec$ramp_m]
      if (length(lead)) {
        cloud$yield[lead] <- cloud$yield[lead] * spec$ramp_low
        cloud$defect[lead] <- "pass_start"
      }
      if (length(tail_i)) {
        cloud$yield[tail_i] <- cloud$yield[tail_i] * spec$ramp_high
        cloud$defect[tail_i] <- "pass_end"
      }
      # The machine slows into the turn at each end of the pass.
      slow <- i[d_along <= spec$slow_m | d_along >= max(d_along) - spec$slow_m]
      cloud$speed[slow] <- cloud$speed[slow] * spec$slow_factor
    }

    n <- nrow(cloud)
    hit <- sample.int(n, max(0L, round(spec$outlier_prob * n)))
    if (length(hit)) {
      zero <- hit[seq_len(round(length(hit) / 2))]
      spike <- setdiff(hit, zero)
      cloud$yield[zero] <- 0
      cloud$yield[spike] <- cloud$yield[spike] * spec$spike_factor
      cloud$defect[hit] <- "outlier"
    }
    rownames(cloud) <- NULL
  }
  cloud$seq <- seq_len(nrow(cloud))
  rownames(cloud) <- NULL

  if (point_design == "strip_ends") {
    bands <- sort(unique(floor(cloud$y / band_width)))
    per_end <- max(1L, round(n_point_samples / (2 * length(bands))))
    si <- unlist(lapply(bands, function(b) {
      inb <- which(floor(cloud$y / band_width) == b)
      if (!length(inb)) return(integer(0))
      xs <- cloud$x[inb]
      lo <- inb[xs <= stats::quantile(xs, 0.10)]
      hi <- inb[xs >= stats::quantile(xs, 0.90)]
      c(if (length(lo)) sample(lo, min(per_end, length(lo))) else integer(0),
        if (length(hi)) sample(hi, min(per_end, length(hi))) else integer(0))
    }))
    si <- si[seq_len(min(length(si), n_point_samples))]
  } else {
    si <- sample(nrow(cloud), min(n_point_samples, nrow(cloud)))
  }

  pts <- data.frame(x = cloud$x[si], y = cloud$y[si],
                    point_true = cloud$point_true[si])
  pts$point_obs <- pts$point_true +
    stats::rnorm(nrow(pts), 0, point_obs_noise_sd)
  rownames(pts) <- NULL

  list(cloud = cloud, point_samples = pts, true_effects = treat_effects)
}

#' Defect settings for simulate_yield_monitor()
#'
#' `defects = TRUE` gives values in the range reported for commercial monitors;
#' a named list overrides any of them.
#' @keywords internal
#' @noRd
.ym_defect_spec <- function(defects) {
  d <- list(ramp_m = 12, ramp_low = 0.55, ramp_high = 1.35,
            slow_m = 8, slow_factor = 0.45,
            outlier_prob = 0.02, spike_factor = 4)
  if (isTRUE(defects)) return(d)
  if (!is.list(defects)) {
    stop("`defects` must be TRUE, FALSE, or a named list of settings.",
         call. = FALSE)
  }
  unknown <- setdiff(names(defects), names(d))
  if (length(unknown) > 0L) {
    stop("Unknown `defects` setting(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(d), collapse = ", "), ".", call. = FALSE)
  }
  utils::modifyList(d, defects)
}
