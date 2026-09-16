# Deriving pseudo-environments from the dense layer, and turning them into a
# residual structure the fitting engines will accept.

#' One-dimensional exponential variogram range
#'
#' A deliberately small, dependency-free variogram fit. The profile is a short,
#' regularly spaced series, so a grid search over the range parameter beats
#' `nls()`, which fails often enough on 20-40 points to be a nuisance inside an
#' automatic procedure. Two details matter and are easy to get wrong: the
#' nugget and partial sill are constrained to be non-negative (an unconstrained
#' fit will happily buy a tiny range with a negative nugget, collapsing the
#' estimate), and lags are weighted by their pair count, so the noisy long lags
#' do not dominate.
#' @keywords internal
#' @noRd
.profile_range <- function(v, spacing, extent) {
  m <- length(v)
  if (m < 6L) return(spacing)
  h_max <- max(3L, floor(m / 2))
  h <- seq_len(h_max)
  n_pair <- m - h
  gamma <- vapply(h, function(k) {
    d <- v[seq.int(1L, m - k)] - v[seq.int(1L + k, m)]
    0.5 * mean(d^2, na.rm = TRUE)
  }, numeric(1))
  if (!all(is.finite(gamma)) || stats::var(gamma) == 0) return(spacing)

  a <- .fit_exp_variogram(h, gamma, n_pair)
  if (is.na(a)) return(spacing)

  # Practical range of an exponential variogram: the lag at which the
  # correlation has decayed to 0.05. Capped at half the trial, because a range
  # wider than that is a statement that there is nothing to partition -- and
  # because a strong step in the profile inflates the variogram without ever
  # reaching a sill, which would otherwise read as an infinite range.
  min(3 * a * spacing, extent / 2)
}

#' Range parameter of an exponential variogram, by constrained grid search
#'
#' The nugget and partial sill enter linearly once the range is fixed, so the
#' search is one-dimensional. Both are constrained non-negative: unconstrained,
#' the fit will buy an arbitrarily small range with a negative nugget and the
#' estimate collapses. Lags are weighted by their pair count so the noisy long
#' lags do not dominate.
#'
#' @param h Numeric vector of lags.
#' @param gamma Numeric vector of empirical semivariances at those lags.
#' @param n_pair Numeric vector of pair counts.
#' @return The fitted range parameter, or `NA_real_` if nothing fits.
#' @keywords internal
#' @noRd
.fit_exp_variogram <- function(h, gamma, n_pair) {
  keep <- is.finite(h) & is.finite(gamma) & is.finite(n_pair) & n_pair > 0
  h <- h[keep]; gamma <- gamma[keep]; n_pair <- n_pair[keep]
  if (length(h) < 3L || stats::var(gamma) == 0) return(NA_real_)
  a_grid <- exp(seq(log(0.25 * min(h)), log(max(h)), length.out = 150))
  ss <- vapply(a_grid, function(a) {
    b <- 1 - exp(-h / a)
    X <- cbind(1, b)
    cf <- tryCatch(
      as.numeric(solve(crossprod(X * n_pair, X), crossprod(X * n_pair, gamma))),
      error = function(e) c(NA_real_, NA_real_))
    if (any(!is.finite(cf))) return(Inf)
    if (cf[2] < 0) return(Inf)
    if (cf[1] < 0) {
      cf[1] <- 0
      cf[2] <- sum(n_pair * b * gamma) / sum(n_pair * b * b)
    }
    sum(n_pair * (gamma - (cf[1] + cf[2] * b))^2)
  }, numeric(1))
  if (all(!is.finite(ss))) return(NA_real_)
  a_grid[which.min(ss)]
}

#' Optimal contiguous segmentation of a profile
#'
#' Dynamic programming over within-segment sums of squares. Exact, and cheap
#' at OFE scale (tens of positions), unlike the greedy splits that
#' quantile-based zoning amounts to.
#' @keywords internal
#' @noRd
.segment_profile <- function(v, w, k, min_pos) {
  m <- length(v)
  cw <- c(0, cumsum(w))
  cs <- c(0, cumsum(w * v))
  cq <- c(0, cumsum(w * v^2))
  cost <- function(i, j) {
    sw <- cw[j + 1L] - cw[i]
    if (sw <= 0) return(0)
    ss <- cq[j + 1L] - cq[i]
    s <- cs[j + 1L] - cs[i]
    max(ss - s^2 / sw, 0)
  }
  INF <- .Machine$double.xmax / 4
  F <- matrix(INF, nrow = k, ncol = m)
  B <- matrix(NA_integer_, nrow = k, ncol = m)
  for (j in seq_len(m)) if (j >= min_pos) F[1L, j] <- cost(1L, j)
  if (k >= 2L) {
    for (q in 2:k) {
      for (j in seq_len(m)) {
        if (j < q * min_pos) next
        lo <- q * min_pos - min_pos + 1L
        hi <- j - min_pos + 1L
        if (hi < lo) next
        for (i in lo:hi) {
          val <- F[q - 1L, i - 1L] + cost(i, j)
          if (val < F[q, j]) {
            F[q, j] <- val
            B[q, j] <- i
          }
        }
      }
    }
  }
  if (!is.finite(F[k, m]) || F[k, m] >= INF) return(NULL)
  ends <- integer(k)
  j <- m
  for (q in seq(k, 1L)) {
    ends[q] <- j
    if (q == 1L) break
    j <- B[q, j] - 1L
  }
  list(ends = ends, sse = F[k, m])
}

#' Derive pseudo-environments from the spatial pattern of the dense layer
#'
#' Splits a trial into contiguous pseudo-environments (PEs) along its length,
#' so that a treatment effect can be estimated separately within each. The cut
#' points come from the data rather than from a guess: the treatment signal is
#' removed, the residual field is collapsed to a profile along the chosen axis,
#' and the profile is segmented optimally by dynamic programming.
#'
#' The role of the spatial covariance is to stop the procedure inventing zones.
#' A smooth field will always look like it has "regions"; the question is
#' whether a region is wider than the correlation range, because anything
#' narrower is one realisation of the same correlated surface rather than a
#' distinct environment. So the practical range of an exponential variogram
#' fitted to the profile becomes the minimum zone width, and it also caps the
#' number of zones at `floor(trial length / range)`. Within that cap the number
#' of zones is chosen by BIC.
#'
#' Two things are worth knowing before trusting the output.
#'
#' First, a smooth field will be split even when nothing discrete is there. A
#' paddock whose yield varies continuously genuinely does have a better end and
#' a worse end, and the procedure will say so. On simulated fields with no step
#' at all and a practical range of about a third of the trial, the default
#' settings return one zone about a third of the time and two or three the
#' rest; with a real step of two-and-a-half field standard deviations they find
#' it, at the right place, essentially always. Read a zone as \dQuote{this part
#' of the paddock behaves differently}, not as evidence of a boundary.
#'
#' Second, PEs derived from the same yield data that are then used to estimate
#' zone-specific treatment effects will overstate those differences, because
#' the boundaries were placed where the residuals already differed. Use them
#' for the residual structure (via [adaptive_residual()]) with a clear
#' conscience; treat a zone-by-treatment interaction fitted on self-derived
#' zones as exploratory, and confirm it against an independent layer --
#' elevation, EM38, a prior season's yield -- by passing that layer as
#' `response` instead.
#'
#' @param data Data frame, one row per lattice cell (e.g. from
#'   [grid_dense_layer()]).
#' @param response Character; the column the zones are derived from. Passing an
#'   independent layer (elevation, EM38, last season's yield) rather than this
#'   season's yield avoids the selection problem described above.
#' @param along Character; the coordinate the trial is cut across -- normally
#'   the direction the strips run, so that every treatment appears in every
#'   zone. Typically `"row"` for strips running up the paddock.
#' @param treat Character or `NULL`; a treatment column to remove before
#'   zoning, so that the treatment pattern is not mistaken for a spatial one.
#' @param block Character or `NULL`; a further nuisance factor to remove
#'   (e.g. `"rep"`).
#' @param method `"segment"` (default) for the data-driven segmentation
#'   described above, or `"equal"` for equal-width slices, which is the
#'   convention this is meant to be compared against.
#' @param n_zones Integer or `NULL`; force a number of zones instead of
#'   choosing it by BIC. Required for `method = "equal"`.
#' @param max_zones Integer; upper limit on the number of zones considered.
#' @param min_zone_width Numeric or `NULL`; minimum zone extent in the units of
#'   `along`. Defaults to the fitted practical range.
#' @param name Character; name of the zone column added to `data`.
#'
#' @return `data` with an added factor column (named by `name`) giving the
#'   pseudo-environment of each cell. Details of the derivation are attached as
#'   `attr(, "partition")`: the fitted `range`, the `min_zone_width` used, the
#'   `breaks`, the per-zone summary, the BIC table, and the profile itself.
#'
#' @seealso [adaptive_residual()], which turns the zones into a `dsum()`
#'   residual formula, and [fit_ofe()], which fits it.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 40, n_col = 12, point_range = 6, seed = 3)
#' z <- partition_pseudo_env(sim$grid, response = "dense_response",
#'                           along = "row", treat = "treat")
#' table(z$zone)
#' attr(z, "partition")$range
#'
#' @export
partition_pseudo_env <- function(data, response, along = "row", treat = NULL,
                                 block = NULL, method = c("segment", "equal"),
                                 n_zones = NULL, max_zones = 8L,
                                 min_zone_width = NULL, name = "zone") {
  method <- match.arg(method)
  data <- as.data.frame(data)
  needed <- c(response, along, treat, block)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  pos_raw <- data[[along]]
  if (is.factor(pos_raw) || is.character(pos_raw)) {
    pos_raw <- suppressWarnings(as.numeric(as.character(pos_raw)))
  }
  if (!is.numeric(pos_raw) || all(is.na(pos_raw))) {
    stop("`along` must name a numeric (or numeric-valued factor) coordinate; ",
         "`", along, "` is not one.", call. = FALSE)
  }

  # Strip the treatment pattern first: a strip trial's treatment effect runs
  # across the strips, and left in place it would be read as spatial structure.
  resp <- data[[response]]
  rhs <- c(treat, block)
  if (length(rhs) > 0L) {
    f <- stats::as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
    m <- stats::lm(f, data = data, na.action = stats::na.exclude)
    resid_v <- stats::residuals(m)
  } else {
    resid_v <- resp - mean(resp, na.rm = TRUE)
  }

  ok <- !is.na(resid_v) & !is.na(pos_raw)
  agg <- stats::aggregate(list(value = resid_v[ok]),
                          by = list(pos = pos_raw[ok]),
                          FUN = mean, na.rm = TRUE)
  cnt <- stats::aggregate(list(n = resid_v[ok]), by = list(pos = pos_raw[ok]),
                          FUN = length)
  agg <- merge(agg, cnt, by = "pos")
  agg <- agg[order(agg$pos), ]
  m <- nrow(agg)
  if (m < 2L) stop("`along` has fewer than two distinct positions.",
                   call. = FALSE)
  spacing <- stats::median(diff(agg$pos))
  extent <- diff(range(agg$pos)) + spacing

  rng <- .profile_range(agg$value, spacing, extent)
  if (is.null(min_zone_width)) min_zone_width <- rng
  min_zone_width <- max(min_zone_width, spacing)
  min_pos <- max(1L, as.integer(ceiling(min_zone_width / spacing)))
  k_cap <- max(1L, min(as.integer(max_zones), floor(m / min_pos)))

  bic_tab <- NULL
  if (method == "equal") {
    if (is.null(n_zones)) {
      stop("`method = \"equal\"` needs `n_zones`.", call. = FALSE)
    }
    k <- as.integer(n_zones)
    brk <- seq(min(agg$pos) - spacing / 2, max(agg$pos) + spacing / 2,
               length.out = k + 1L)
  } else {
    k_try <- seq_len(k_cap)
    fits <- lapply(k_try, function(k) .segment_profile(agg$value, agg$n, k,
                                                       min_pos))
    keep <- !vapply(fits, is.null, logical(1))
    if (!any(keep)) {
      k <- 1L
      brk <- c(min(agg$pos) - spacing / 2, max(agg$pos) + spacing / 2)
    } else {
      k_try <- k_try[keep]
      fits <- fits[keep]
      n_obs <- sum(agg$n)
      # Put the weighted SSE back on the scale of the profile: the profile is
      # what is being segmented, and its m points are averages, so counting all
      # n_obs cells as independent evidence would penalise extra zones far too
      # lightly and split a smooth field into as many pieces as it is allowed.
      sse <- vapply(fits, function(f) f$sse, numeric(1)) / n_obs * m
      bic <- m * log(pmax(sse, .Machine$double.eps) / m) + (2 * k_try) * log(m)
      bic_tab <- data.frame(n_zones = k_try, sse = sse, bic = bic)
      if (!is.null(n_zones)) {
        sel <- match(as.integer(n_zones), k_try)
        if (is.na(sel)) {
          stop("`n_zones = ", n_zones, "` is not reachable: at most ", k_cap,
               " zone(s) fit while keeping each at least ",
               format(min_zone_width, digits = 3), " wide.", call. = FALSE)
        }
      } else {
        sel <- which.min(bic)
      }
      k <- k_try[sel]
      ends <- fits[[sel]]$ends
      cut_at <- agg$pos[ends[-length(ends)]]
      brk <- c(min(agg$pos) - spacing / 2, cut_at + spacing / 2,
               max(agg$pos) + spacing / 2)
    }
  }

  zone <- cut(pos_raw, breaks = brk, labels = FALSE, include.lowest = TRUE)
  zone[is.na(zone)] <- if (k == 1L) 1L else NA_integer_
  data[[name]] <- factor(zone, levels = seq_len(k))

  summ <- data.frame(
    zone = factor(seq_len(k)),
    start = brk[-length(brk)],
    end = brk[-1L],
    width = diff(brk),
    n = as.integer(table(factor(zone, levels = seq_len(k))))
  )

  attr(data, "partition") <- list(
    method = method, along = along, response = response,
    n_zones = k, range = rng, min_zone_width = min_zone_width,
    spacing = spacing, extent = extent, max_zones_possible = k_cap,
    breaks = brk, zones = summ, bic = bic_tab,
    profile = data.frame(pos = agg$pos, value = agg$value, n = agg$n)
  )
  data
}

#' Build a residual structure that matches the realised geometry
#'
#' Writes the `dsum()` residual formula for a pseudo-environment analysis, one
#' independent section per zone. The reason this needs a function rather than a
#' literal formula is that AR1 needs at least two levels in a dimension, and
#' real zones are not guaranteed to have them: a narrow zone may be one row
#' deep, and asking for `ar1()` there gives an unidentifiable parameter and a
#' failed fit. Each section therefore gets `ar1()` only in the dimensions where
#' it actually has extent, and `id()` elsewhere.
#'
#' The formula is plain text that both [fit_ofe()] and `asreml::asreml()`
#' accept, so the same analysis script runs with or without a licence.
#'
#' @param data Data frame containing the zone and position columns.
#' @param zone Character; the pseudo-environment column, e.g. from
#'   [partition_pseudo_env()].
#' @param row,col Character; the two lattice dimensions.
#'
#' @return A one-sided formula, with the per-zone geometry attached as
#'   `attr(, "geometry")` and the number of zones that had to be demoted to
#'   `id()` as `attr(, "n_degenerate")`.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 40, n_col = 12, point_range = 6, seed = 3)
#' z <- partition_pseudo_env(sim$grid, response = "dense_response",
#'                           along = "row", treat = "treat")
#' r <- adaptive_residual(z)
#' r
#' attr(r, "geometry")
#'
#' @export
adaptive_residual <- function(data, zone = "zone", row = "row", col = "col") {
  data <- as.data.frame(data)
  missing_cols <- setdiff(c(zone, row, col), names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  z <- factor(data[[zone]])
  lv <- levels(droplevels(z))
  if (length(lv) == 0L) stop("`", zone, "` has no levels.", call. = FALSE)

  geo <- do.call(rbind, lapply(lv, function(l) {
    i <- which(as.character(z) == l)
    nr <- length(unique(data[[row]][i]))
    nc <- length(unique(data[[col]][i]))
    data.frame(zone = l, n = length(i), n_row = nr, n_col = nc,
               struct = paste0(if (nr >= 2L) paste0("ar1(", row, ")") else
                 paste0("id(", row, ")"), ":",
                 if (nc >= 2L) paste0("ar1(", col, ")") else
                   paste0("id(", col, ")")),
               stringsAsFactors = FALSE)
  }))

  # One dsum() per distinct structure, listing the zones that share it, which
  # is the compact spelling asreml's refman uses.
  parts <- vapply(unique(geo$struct), function(s) {
    ls <- geo$zone[geo$struct == s]
    sprintf("dsum(~ %s | %s, levels = c(%s))", s, zone,
            paste0('"', ls, '"', collapse = ", "))
  }, character(1), USE.NAMES = FALSE)

  f <- stats::as.formula(paste("~", paste(parts, collapse = " + ")),
                         env = parent.frame())
  attr(f, "geometry") <- geo
  attr(f, "n_degenerate") <- sum(geo$n_row < 2L | geo$n_col < 2L)
  f
}
