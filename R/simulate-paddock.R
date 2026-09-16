# A paddock with a known truth: a correlated yield-potential surface, covariate
# layers that are related to it, pseudo-environments with their own treatment
# response, and a trial laid into it.
#
# The package already simulated trials, but never a paddock: no covariate
# layers, so partition_paddock() had nothing to zone; no zones, so
# partition_pseudo_env() had no truth to be checked against; and no separation
# between yield potential and the treatment, so nothing to attribute.

#' Draw one spatially correlated surface
#' @keywords internal
#' @noRd
.sim_surface <- function(grid, range, psill, nugget = 0, beta = 0) {
  vg <- gstat::vgm(psill = psill, model = "Exp", range = range,
                   nugget = nugget)
  g <- gstat::gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE,
                    beta = beta, model = vg, nmax = 20)
  out <- stats::predict(g, newdata = grid[, c("x", "y")], nsim = 1,
                        debug.level = 0)
  as.numeric(out$sim1)
}

#' Simulate a paddock, its covariate layers, and a trial in it
#'
#' Generates a paddock whose truth is known, so that the rest of the package can
#' be checked rather than merely run: the zones it should find, the treatment
#' effect it should recover, and how much of the yield was never attributable to
#' treatment in the first place.
#'
#' # What is simulated
#'
#' A spatially correlated **yield potential** surface is drawn first; it is the
#' paddock the grower already has. Each named **covariate** is then drawn as a
#' mixture of that surface and its own independent structure, in the proportion
#' set by `covariate_cor` -- which is the point of them. A covariate unrelated
#' to yield is not worth zoning on, and one identical to it is an unrealistically
#' easy test; the default puts elevation and soil part-way, as they are.
#'
#' Optional **pseudo-environments** cut the paddock into bands across the
#' strips. Each band gets its own treatment response through `zone_effect`, so
#' a zone-by-treatment interaction is genuinely there to be found -- or, with
#' the default of one zone, genuinely is not, which is the more important case
#' to be able to test.
#'
#' Given `treatments`, a trial is laid out with [make_trial_design()] and
#' harvested: yield is potential, plus the treatment response for that cell's
#' zone, plus noise.
#'
#' # Using it
#'
#' The truth is attached as `attr(, "truth")` -- the zone breaks, the per-zone
#' treatment effects, the variance of each component -- so a check can be made
#' against what was simulated rather than against what looks plausible. That is
#' what separates this from a demonstration dataset.
#'
#' @param n_row,n_col Lattice size of the paddock.
#' @param cell_size Cell size in metres.
#' @param yield_mean Mean yield potential, in whatever unit you are working in
#'   (t/ha, say).
#' @param yield_range Practical range of the yield-potential surface, in metres.
#'   Short relative to the paddock gives a patchy field; long gives a gradient.
#' @param yield_sd Standard deviation of the yield-potential surface.
#' @param covariates Named list of covariate layers to draw. Each element is a
#'   list which may name `mean`, `sd`, `range` and `cor` (its correlation with
#'   yield potential). Defaults supply `elevation` and `clay`. Use `NULL` for
#'   none.
#' @param covariate_cor Default correlation between a covariate and yield
#'   potential, for layers that do not set their own `cor`.
#' @param n_zones Integer; number of pseudo-environments across the paddock.
#'   `1` (the default) means no zone structure.
#' @param zone_props Optional proportions of the paddock length given to each
#'   zone; defaults to equal bands. Real boundaries are not tidy, so unequal
#'   proportions are the more realistic test.
#' @param treatments Optional treatment labels. When given, a trial is laid out
#'   and harvested; when `NULL`, only the paddock and its layers are returned.
#' @param n_rep Replicates, passed to [make_trial_design()].
#' @param treat_effect Numeric vector of treatment effects, one per treatment,
#'   in yield units. Defaults to an evenly spaced response.
#' @param zone_effect Numeric vector of multipliers, one per zone, scaling the
#'   treatment effect in that zone. `rep(1, n_zones)` means the treatment works
#'   equally everywhere.
#' @param noise_sd Standard deviation of the independent harvest noise.
#' @param randomise Passed to [make_trial_design()].
#' @param seed Optional random seed.
#'
#' @return A data frame with one row per cell: `row`, `col`, `x`, `y`,
#'   `yield_potential`, one column per covariate, `zone`, and -- when
#'   `treatments` was given -- `treat`, `rep`, `plot` and `yield`. Classed as
#'   `ofe_zones` so `plot()` maps the zones. `attr(, "truth")` holds the
#'   settings and the realised component variances.
#'
#' @seealso [simulate_ofe_trial()] for a bare lattice and
#'   [simulate_yield_monitor()] for the irregular shape a harvester records.
#'
#' @examples
#' p <- simulate_paddock(n_row = 30, n_col = 20, n_zones = 2,
#'                       treatments = c("N0", "N60", "N120"), seed = 1)
#' head(p)
#' attr(p, "truth")$zone_effect
#'
#' # The zoning functions can now be checked against a known answer
#' z <- partition_paddock(p, covariates = c("elevation", "clay"))
#' table(z$zone, p$zone)
#'
#' @export
simulate_paddock <- function(n_row = 40L, n_col = 24L, cell_size = 10,
                             yield_mean = 3, yield_range = 80, yield_sd = 0.6,
                             covariates = NULL, covariate_cor = 0.6,
                             n_zones = 1L, zone_props = NULL,
                             treatments = NULL, n_rep = 4L,
                             treat_effect = NULL, zone_effect = NULL,
                             noise_sd = 0.3, randomise = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_zones <- max(1L, as.integer(n_zones))
  if (is.null(covariates)) {
    covariates <- list(
      elevation = list(mean = 100, sd = 3, range = yield_range * 1.5),
      clay = list(mean = 30, sd = 6, range = yield_range * 0.8))
  }
  if (length(covariates) > 0L && is.null(names(covariates))) {
    stop("`covariates` must be a named list, one element per layer.",
         call. = FALSE)
  }

  grid <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  grid$x <- grid$col * cell_size
  grid$y <- grid$row * cell_size
  n <- nrow(grid)

  # Yield potential: the paddock the grower already has, before any treatment.
  pot <- .sim_surface(grid, range = yield_range / 3, psill = 1)
  pot <- (pot - mean(pot)) / stats::sd(pot)
  grid$yield_potential <- yield_mean + yield_sd * pot

  for (nm in names(covariates)) {
    spec <- covariates[[nm]]
    rho <- if (!is.null(spec$cor)) spec$cor else covariate_cor
    if (rho < -1 || rho > 1) {
      stop("`cor` for covariate `", nm, "` must be between -1 and 1.",
           call. = FALSE)
    }
    rng <- if (!is.null(spec$range)) spec$range else yield_range
    own <- .sim_surface(grid, range = rng / 3, psill = 1)
    own <- (own - mean(own)) / stats::sd(own)
    # Mix the shared surface with the layer's own, so the correlation with
    # yield potential is what was asked for rather than whatever falls out.
    mixed <- rho * pot + sqrt(max(0, 1 - rho^2)) * own
    mu <- if (!is.null(spec$mean)) spec$mean else 0
    sdv <- if (!is.null(spec$sd)) spec$sd else 1
    grid[[nm]] <- mu + sdv * as.numeric(scale(mixed))
  }

  if (is.null(zone_props)) zone_props <- rep(1 / n_zones, n_zones)
  if (length(zone_props) != n_zones) {
    stop("`zone_props` must have one element per zone.", call. = FALSE)
  }
  zone_props <- zone_props / sum(zone_props)
  len <- n_row * cell_size
  breaks <- c(0, cumsum(zone_props) * len)
  breaks[length(breaks)] <- len + cell_size
  grid$zone <- factor(cut(grid$y, breaks = breaks, labels = FALSE,
                          include.lowest = TRUE), levels = seq_len(n_zones))

  truth <- list(yield_mean = yield_mean, yield_range = yield_range,
                yield_sd = yield_sd, n_zones = n_zones,
                zone_breaks = breaks, zone_props = zone_props,
                covariate_cor = vapply(names(covariates), function(nm)
                  stats::cor(grid[[nm]], grid$yield_potential), numeric(1)),
                noise_sd = noise_sd, cell_size = cell_size)

  if (!is.null(treatments)) {
    n_t <- length(treatments)
    if (is.null(treat_effect)) {
      treat_effect <- seq(0, 0.8, length.out = n_t)
    }
    if (length(treat_effect) != n_t) {
      stop("`treat_effect` must have one element per treatment.",
           call. = FALSE)
    }
    if (is.null(zone_effect)) zone_effect <- rep(1, n_zones)
    if (length(zone_effect) != n_zones) {
      stop("`zone_effect` must have one element per zone.", call. = FALSE)
    }
    plot_width <- floor(n_col / (n_t * n_rep)) * cell_size
    if (plot_width < cell_size) {
      stop("The paddock is only ", n_col, " cells across, too narrow for ",
           n_t, " treatments in ", n_rep, " reps. Widen it, or use fewer ",
           "reps.", call. = FALSE)
    }
    des <- make_trial_design(treatments, n_rep = n_rep, plot_width = plot_width,
                             plot_length = n_row * cell_size,
                             cell_size = cell_size, randomise = randomise)
    des <- as.data.frame(des)[, c("row", "col", "treat", "rep", "plot")]
    grid <- merge(grid, des, by = c("row", "col"), all.x = TRUE,
                  sort = FALSE)
    grid <- grid[order(grid$row, grid$col), ]

    names(treat_effect) <- treatments
    eff <- treat_effect[as.character(grid$treat)] *
      zone_effect[as.integer(grid$zone)]
    eff[is.na(grid$treat)] <- NA_real_
    grid$yield <- grid$yield_potential + eff +
      stats::rnorm(nrow(grid), 0, noise_sd)

    truth$treatments <- treatments
    truth$treat_effect <- treat_effect
    truth$zone_effect <- zone_effect
    # The number that matters for interpreting any later analysis: how much of
    # the yield spread was never the treatment's doing.
    truth$var_potential <- stats::var(grid$yield_potential)
    truth$var_treatment <- stats::var(eff, na.rm = TRUE)
    truth$var_noise <- noise_sd^2
  }

  rownames(grid) <- NULL
  attr(grid, "truth") <- truth
  attr(grid, "partition") <- list(method = "simulated", n_zones = n_zones,
                                  covariates = "zone", breaks = breaks)
  class(grid) <- unique(c("ofe_zones", class(grid)))
  grid
}
