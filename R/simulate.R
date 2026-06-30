#' Simulate a synthetic OFE strip trial with a sparse point-source covariate
#'
#' Generates a rectangular grid trial with treatment strips, a spatially
#' correlated "dense" response (e.g. yield), and a spatially correlated
#' "point-source" covariate (e.g. soil, tissue, or disease measurements)
#' that is only observed at a sparse set of sampled locations. Intended for
#' testing and demonstrating [krige_point_samples()], [fit_integrated_kriged()],
#' and [fit_integrated_joint()].
#'
#' @param n_row,n_col Integer grid dimensions (rows = one axis, columns =
#'   the axis along which treatment strips run).
#' @param n_treat Integer number of treatment strips, applied as
#'   contiguous, approximately equal-width column blocks.
#' @param treat_effects Numeric vector of length `n_treat` giving the true
#'   treatment effect added to the dense response. Recycled with treatment
#'   labels `LETTERS[1:n_treat]`.
#' @param point_range,point_psill,point_nugget Variogram parameters
#'   (exponential model) for the spatially correlated point-source surface.
#' @param dense_var_weight Numeric; coefficient relating the point-source
#'   surface to the dense response (i.e. how much of the dense response's
#'   spatial pattern is explained by the point-source variable).
#' @param noise_sd Standard deviation of i.i.d. measurement noise added to
#'   the dense response.
#' @param n_point_samples Integer number of sparse point-source samples to
#'   draw from the grid.
#' @param point_obs_noise_sd Standard deviation of i.i.d. measurement
#'   (e.g. lab) noise added to the observed point-source samples.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with components:
#'   \describe{
#'     \item{grid}{Data frame, one row per grid cell, with `row`, `col`,
#'       `treat`, `point_true` (true point-source surface) and
#'       `dense_response` (e.g. simulated yield).}
#'     \item{point_samples}{Data frame of `n_point_samples` sparse
#'       observations with `row`, `col`, `point_true`, and `point_obs`
#'       (true value plus measurement noise).}
#'     \item{true_effects}{Named numeric vector of the true treatment
#'       effects used to simulate `dense_response`.}
#'   }
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 8, seed = 1)
#' head(sim$grid)
#' sim$point_samples
#' sim$true_effects
#'
#' @export
simulate_ofe_trial <- function(n_row = 40,
                                n_col = 20,
                                n_treat = 3,
                                treat_effects = c(0, 0.8, 1.6),
                                point_range = 8,
                                point_psill = 1.2,
                                point_nugget = 0.1,
                                dense_var_weight = 0.6,
                                noise_sd = 0.4,
                                n_point_samples = 10,
                                point_obs_noise_sd = 0.3,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(treat_effects) != n_treat) {
    stop("`treat_effects` must have length `n_treat`.", call. = FALSE)
  }

  grid <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))

  breaks <- seq(0, n_col, length.out = n_treat + 1)
  grid$treat <- cut(grid$col, breaks = breaks,
                     labels = LETTERS[seq_len(n_treat)],
                     include.lowest = TRUE)

  names(treat_effects) <- LETTERS[seq_len(n_treat)]

  point_vgm <- gstat::vgm(psill = point_psill, model = "Exp",
                           range = point_range, nugget = point_nugget)
  g <- gstat::gstat(formula = z ~ 1, locations = ~ col + row, dummy = TRUE,
                     beta = 0, model = point_vgm, nmax = 30)
  sim_out <- stats::predict(g, newdata = grid[, c("col", "row")], nsim = 1,
                             debug.level = 0)
  grid$point_true <- sim_out$sim1

  grid$dense_response <- treat_effects[as.character(grid$treat)] +
    dense_var_weight * grid$point_true +
    stats::rnorm(nrow(grid), 0, noise_sd)

  idx <- sample(nrow(grid), n_point_samples)
  point_samples <- grid[idx, c("row", "col", "point_true")]
  point_samples$point_obs <- point_samples$point_true +
    stats::rnorm(n_point_samples, 0, point_obs_noise_sd)

  list(
    grid = grid,
    point_samples = point_samples,
    true_effects = treat_effects
  )
}
