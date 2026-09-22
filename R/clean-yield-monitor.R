# Cleaning a raw yield-monitor file before anything else touches it.
#
# grid_dense_layer() assumes the cloud it is given is sound. A file straight off
# a harvester is not: the opening metres of every pass read low while grain is
# still reaching the sensor, the closing metres read high, the machine slows
# into each turn, and a scatter of points are zero or several times the crop.
# Averaging those into cells does not cancel them -- the pass-end effects are
# systematic and sit at the same end of every pass, so they survive aggregation
# and land in the model as spatial structure.

#' Split a point cloud into harvester passes
#'
#' Uses the recording order where there is one -- a pass ends where the machine
#' turns or jumps -- and falls back on geometry otherwise: passes are parallel,
#' so projecting onto the axis across them turns the problem into finding gaps
#' in one dimension.
#' @keywords internal
#' @noRd
.derive_passes <- function(xv, yv, order = NULL, gap_mult = 4) {
  n <- length(xv)
  if (n < 3L) return(rep(1L, n))

  # Runs shorter than this are absorbed into their neighbour. A turn measured
  # over a window fires on every point whose window straddles the turn, so a
  # single pass boundary otherwise produces several breaks a few points apart
  # and leaves slivers behind.
  tidy_runs <- function(p, min_run = 5L) {
    r <- rle(p)
    while (any(r$lengths < min_run) && length(r$lengths) > 1L) {
      k <- which.min(r$lengths)
      r$values[k] <- if (k == 1L) r$values[2L] else r$values[k - 1L]
      p <- inverse.rle(r)
      p <- cumsum(c(TRUE, diff(p) != 0))
      r <- rle(p)
    }
    as.integer(p)
  }

  if (!is.null(order)) {
    o <- base::order(order)
    ox <- xv[o]
    oy <- yv[o]
    step <- sqrt(diff(ox)^2 + diff(oy)^2)
    typical <- stats::median(step[step > 0])

    # A pass ends where the machine jumps to the next swath, or turns.
    is_break <- c(FALSE, step > gap_mult * typical)

    # The turn has to be measured over a window: between two consecutive
    # points the bearing is mostly GPS noise. Even windowed it fires on every
    # point whose window spans the turn, which is what tidy_runs() cleans up.
    w <- 3L
    if (n > 2L * w + 1L) {
      i <- seq.int(w + 1L, n - w)
      v1x <- ox[i] - ox[i - w]; v1y <- oy[i] - oy[i - w]
      v2x <- ox[i + w] - ox[i]; v2y <- oy[i + w] - oy[i]
      denom <- sqrt(v1x^2 + v1y^2) * sqrt(v2x^2 + v2y^2)
      cosang <- ifelse(denom > 0, (v1x * v2x + v1y * v2y) / denom, 1)
      is_break[i] <- is_break[i] | cosang < 0
    }

    p <- tidy_runs(cumsum(is_break) + 1L)
    out <- integer(n)
    out[o] <- p
    return(out)
  }

  # No recording order: recover the direction of travel from the data. Not by
  # principal components of the positions -- an irregular boundary, a cut
  # corner, puts covariance between x and y and tilts the axes, and over a long
  # trial a tilt of a few degrees smears the across-pass coordinate far enough
  # to hide the swaths entirely. A point's nearest neighbour, on the other
  # hand, is almost always the next point along its own pass, so the bearings
  # to nearest neighbours point along the passes whatever shape the paddock is.
  idx <- if (n > 600L) sample.int(n, 600L) else seq_len(n)
  dm <- as.matrix(stats::dist(cbind(xv[idx], yv[idx])))
  diag(dm) <- Inf
  nn <- max.col(-dm, ties.method = "first")
  ang <- atan2(yv[idx][nn] - yv[idx], xv[idx][nn] - xv[idx])
  # Direction of travel is a line, not an arrow: alternate passes run opposite
  # ways, so bearings are averaged doubled and halved back.
  theta <- atan2(mean(sin(2 * ang)), mean(cos(2 * ang))) / 2
  across <- -sin(theta) * xv + cos(theta) * yv

  o <- base::order(across)
  d <- diff(across[o])
  pos <- sort(d[d > 0], decreasing = TRUE)
  if (length(pos) < 2L) return(rep(1L, n))
  # Between-pass gaps are a swath wide and within-pass gaps are millimetres,
  # so the sorted gaps fall off a cliff. Cut at the cliff rather than at a
  # multiple of the median, which is set by the within-pass spacing and would
  # split every pass into its individual points.
  top <- pos[seq_len(min(length(pos), 60L))]
  ratio <- top[-length(top)] / pmax(top[-1], .Machine$double.eps)
  k <- which.max(ratio)
  thresh <- (top[k] + top[k + 1L]) / 2
  p <- cumsum(c(TRUE, d > thresh))
  out <- integer(n)
  out[o] <- tidy_runs(p)
  out
}

#' Distance from each point to the boundary of the cloud's convex hull
#' @keywords internal
#' @noRd
.hull_distance <- function(xv, yv) {
  h <- grDevices::chull(xv, yv)
  hx <- xv[h]
  hy <- yv[h]
  m <- length(h)
  if (m < 3L) return(rep(0, length(xv)))
  best <- rep(Inf, length(xv))
  for (k in seq_len(m)) {
    j <- if (k == m) 1L else k + 1L
    ax <- hx[k]; ay <- hy[k]
    bx <- hx[j]; by <- hy[j]
    vx <- bx - ax; vy <- by - ay
    len2 <- vx^2 + vy^2
    t <- if (len2 == 0) rep(0, length(xv)) else
      pmin(1, pmax(0, ((xv - ax) * vx + (yv - ay) * vy) / len2))
    dx <- xv - (ax + t * vx)
    dy <- yv - (ay + t * vy)
    best <- pmin(best, sqrt(dx^2 + dy^2))
  }
  best
}

#' Clean a raw yield-monitor point cloud
#'
#' Removes the points a harvester records but did not really measure, before
#' they reach [grid_dense_layer()] and become part of the spatial model.
#'
#' # How well it works
#'
#' Scored against [simulate_yield_monitor()] with `defects = TRUE`, which
#' labels every point it damages, using `pass_trim = 12`, a speed column, and
#' `local_mad = 3`: 97% of the depressed pass starts and 90% of the inflated
#' pass ends are removed, for 1% of the sound points. On a cloud in physical
#' units, where `min_yield = 0` can be set, 94% of the injected outliers go as
#' well. Without that rule only 39% of them do, because a dropout reading of
#' zero is not far enough from the median to trip a dispersion threshold.
#'
#' # Why it cannot be left to averaging
#'
#' The tempting argument is that gridding averages several observations per
#' cell, so bad points wash out. Only some do. Two of the defects below are
#' *systematic*: the opening metres of every pass read low and the closing
#' metres read high, and every pass starts and ends at the same two edges of the
#' paddock. Averaging preserves that, and it arrives in the analysis as a
#' spatial trend along the direction of travel, which a spatial model will
#' happily fit. Only the random outliers are diluted, and they are the ones that
#' matter least.
#'
#' Whether the trend reaches the treatment estimate depends on how the passes
#' lie against the treatment strips. Run the passes along the strips and each
#' one sits inside a single treatment, so the fault falls on all of them
#' equally: on simulated trials of that shape, cleaning moves the recovered
#' contrast by about a percentage point. Run them across the strips and it does
#' not fall equally, and the trend becomes a treatment effect. The case for
#' cleaning regardless is that it costs little and the direction of harvest is
#' often not the analyst's to choose.
#'
#' # The rules
#'
#' Each is optional and each is off unless its argument is set, except the
#' outlier rules, which default to conservative values. Every removed point is
#' attributed to the first rule that caught it.
#'
#' \describe{
#'   \item{`missing`}{No yield, or no position.}
#'   \item{`range`}{Outside `yield_range`. There is no default: the plausible
#'     range depends on the crop and the units, and guessing it would be worse
#'     than not applying the rule.}
#'   \item{`zero`}{At or below `min_yield`. A zero reading while the machine is
#'     cutting is impossible rather than extreme, and it sits well inside any
#'     dispersion-based threshold -- which is why a statistical rule on its own
#'     leaves dropouts in. **Set this for real data**: `min_yield = 0` is right
#'     for a yield in physical units. It is off by default for the same reason
#'     `yield_range` is: the function cannot see the units, and a threshold in
#'     the wrong ones removes good data. The MAD rules need no such care
#'     because they are relative to the spread of the data itself.}
#'   \item{`outlier`}{Further than `outlier_mad` median absolute deviations
#'     from the median yield. The median and MAD are used rather than the mean
#'     and standard deviation because the outliers inflate the standard
#'     deviation, which is how a three-sigma rule ends up keeping the very
#'     points it was meant to remove. Even so this rule only catches gross
#'     spikes: the dispersion it measures includes the field's real variation,
#'     so the threshold has to stay loose.}
#'   \item{`local`}{Further than `local_mad` deviations from the median of the
#'     point's own neighbourhood, on a grid of `local_cell` metres. Taking the
#'     comparison local removes the field's own variation from it, so a much
#'     tighter threshold is safe -- this is the rule that catches a bad point
#'     sitting in a good part of the paddock. Off unless `local_mad` is set.
#'     Use 3: on labelled simulated clouds it lifts the share of injected
#'     outliers removed from 39% to 58% while costing 1% of the good points,
#'     and tightening it to 2 finds no more of them but removes 8%.}
#'   \item{`speed`}{Zero or negative speed, or further than `speed_mad` MADs
#'     from the median. Needs a `speed` column.}
#'   \item{`pass_end`}{Within `pass_trim` metres of either end of a pass, along
#'     the direction of travel. This is the grain-flow rule, and the one worth
#'     setting: a combine's fill time is a few seconds, which at working speed
#'     is of the order of ten metres.}
#'   \item{`edge`}{Within `edge_buffer` metres of the boundary of the cloud's
#'     convex hull -- the headland, where the machine is turning, part-full, or
#'     driving over ground it has already cut.}
#' }
#'
#' @param data Data frame of yield-monitor observations.
#' @param yield,x,y Character; the yield and projected coordinate columns.
#' @param pass Character or `NULL`; a pass or swath identifier. When absent it
#'   is derived -- from `order` if given, otherwise from the geometry, since
#'   passes are parallel.
#' @param order Character or `NULL`; the recording order (a time stamp or a
#'   sequence number). Used to derive passes and to find the ends of each.
#' @param speed Character or `NULL`; a ground-speed column.
#' @param yield_range Numeric of length 2, or `NULL`; plausible yield bounds.
#' @param min_yield Numeric; yields at or below this are treated as sensor
#'   dropouts. Off by default (`-Inf`) because the right value depends on the
#'   units; set it to 0 for a yield in physical units.
#' @param outlier_mad Numeric; MADs from the median beyond which a yield is
#'   dropped. `Inf` turns the rule off.
#' @param local_mad Numeric or `NULL`; deviations from the local median beyond
#'   which a yield is dropped. `NULL` (default) turns the rule off.
#' @param local_cell Numeric or `NULL`; the neighbourhood size in metres for
#'   `local_mad`. Defaults to a twentieth of the field's longer side.
#' @param speed_mad Numeric; the same for speed, when a speed column is given.
#' @param pass_trim Numeric; metres to drop from each end of every pass.
#' @param edge_buffer Numeric; metres to drop in from the field boundary.
#' @param drop Logical; return only the surviving rows (default), or every row
#'   with the verdict attached.
#'
#' @return With `drop = TRUE`, `data` with the removed rows gone. With
#'   `drop = FALSE`, every row, plus `.keep` and `.reason` columns -- feed that
#'   to [ofe_map()] with `value = ".reason"` to see what was taken and from
#'   where. Either way a per-rule tally is attached as `attr(, "cleaning")`.
#'
#' @seealso [grid_dense_layer()], which is what this feeds, and
#'   [simulate_yield_monitor()] with `defects = TRUE`, which generates a cloud
#'   whose bad points are labelled so a cleaning rule can be scored.
#'
#' @examples
#' sim <- simulate_yield_monitor(defects = TRUE, seed = 1)
#' clean <- clean_yield_monitor(sim$cloud, pass = "pass", order = "seq",
#'                              speed = "speed", pass_trim = 12, local_mad = 3)
#' attr(clean, "cleaning")
#'
#' # What was removed, and from where
#' audit <- clean_yield_monitor(sim$cloud, pass = "pass", order = "seq",
#'                              speed = "speed", pass_trim = 12, drop = FALSE)
#' table(audit$.reason)
#'
#' @export
clean_yield_monitor <- function(data, yield = "yield", x = "x", y = "y",
                                pass = NULL, order = NULL, speed = NULL,
                                yield_range = NULL, min_yield = -Inf,
                                outlier_mad = 4, local_mad = NULL,
                                local_cell = NULL, speed_mad = 3,
                                pass_trim = 0, edge_buffer = 0, drop = TRUE) {
  data <- as.data.frame(data)
  needed <- c(yield, x, y, pass, order, speed)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(yield_range) &&
      (length(yield_range) != 2L || !is.numeric(yield_range) ||
       yield_range[1] >= yield_range[2])) {
    stop("`yield_range` must be two increasing numbers.", call. = FALSE)
  }
  n <- nrow(data)
  yv <- as.numeric(data[[yield]])
  xc <- as.numeric(data[[x]])
  yc <- as.numeric(data[[y]])
  reason <- rep(NA_character_, n)

  mark <- function(reason, hit, label) {
    hit[is.na(hit)] <- FALSE
    reason[is.na(reason) & hit] <- label
    reason
  }

  reason <- mark(reason, is.na(yv) | is.na(xc) | is.na(yc), "missing")
  if (!is.null(yield_range)) {
    reason <- mark(reason, yv < yield_range[1] | yv > yield_range[2], "range")
  }
  if (is.finite(min_yield)) {
    # A dropout, not an outlier. A zero reading while the machine is cutting is
    # impossible rather than merely extreme, and it is well inside any
    # dispersion-based threshold -- which is exactly why a statistical rule
    # alone leaves them in.
    reason <- mark(reason, yv <= min_yield, "zero")
  }
  if (is.finite(outlier_mad)) {
    med <- stats::median(yv, na.rm = TRUE)
    disp <- stats::mad(yv, center = med, na.rm = TRUE)
    if (is.finite(disp) && disp > 0) {
      reason <- mark(reason, abs(yv - med) > outlier_mad * disp, "outlier")
    }
  }
  if (!is.null(local_mad)) {
    # A global threshold has to be loose enough to leave the field's own
    # variation alone, which makes it too loose to catch anything subtle.
    # Comparing each point with its own neighbourhood removes that variation
    # from the comparison, so a much tighter threshold is still safe.
    if (is.null(local_cell)) {
      span <- max(diff(range(xc, na.rm = TRUE)), diff(range(yc, na.rm = TRUE)))
      local_cell <- span / 20
    }
    cx <- floor((xc - min(xc, na.rm = TRUE)) / local_cell)
    cy <- floor((yc - min(yc, na.rm = TRUE)) / local_cell)
    cell <- paste(cx, cy, sep = "_")
    med <- stats::ave(yv, cell, FUN = function(z) stats::median(z, na.rm = TRUE))
    dev <- abs(yv - med)
    # One dispersion for the whole field, from the local deviations: a per-cell
    # MAD would be estimated from a handful of points and would itself be noise.
    scale_local <- stats::median(dev, na.rm = TRUE) * 1.4826
    if (is.finite(scale_local) && scale_local > 0) {
      reason <- mark(reason, dev > local_mad * scale_local, "local")
    }
  }

  if (!is.null(speed)) {
    sv <- as.numeric(data[[speed]])
    hit <- !is.na(sv) & sv <= 0
    if (is.finite(speed_mad)) {
      smed <- stats::median(sv, na.rm = TRUE)
      sdisp <- stats::mad(sv, center = smed, na.rm = TRUE)
      if (is.finite(sdisp) && sdisp > 0) {
        hit <- hit | abs(sv - smed) > speed_mad * sdisp
      }
    }
    reason <- mark(reason, hit, "speed")
  }

  pass_id <- NULL
  if (pass_trim > 0) {
    pass_id <- if (!is.null(pass)) as.character(data[[pass]]) else
      as.character(.derive_passes(xc, yc,
                                  order = if (!is.null(order)) data[[order]]))
    from_end <- rep(NA_real_, n)
    for (p in unique(pass_id)) {
      i <- which(pass_id == p)
      if (length(i) < 2L) {
        from_end[i] <- 0
        next
      }
      # Project onto the pass's own direction: the ends of a pass are the
      # extremes along the way it was driven, not along x or y.
      px <- xc[i] - mean(xc[i])
      py <- yc[i] - mean(yc[i])
      ev <- eigen(stats::cov(cbind(px, py)), symmetric = TRUE)
      along <- as.numeric(cbind(px, py) %*% ev$vectors[, 1])
      from_end[i] <- pmin(along - min(along), max(along) - along)
    }
    reason <- mark(reason, from_end < pass_trim, "pass_end")
  }

  if (edge_buffer > 0) {
    ok <- !is.na(xc) & !is.na(yc)
    dist_edge <- rep(NA_real_, n)
    dist_edge[ok] <- .hull_distance(xc[ok], yc[ok])
    reason <- mark(reason, dist_edge < edge_buffer, "edge")
  }

  keep <- is.na(reason)
  rules <- c("missing", "range", "zero", "outlier", "local", "speed",
             "pass_end", "edge")
  tally <- data.frame(
    rule = rules,
    removed = as.integer(vapply(rules, function(r) sum(reason == r,
                                                       na.rm = TRUE),
                                integer(1))),
    stringsAsFactors = FALSE)
  tally$percent <- round(100 * tally$removed / n, 2)
  tally <- rbind(tally,
                 data.frame(rule = "kept", removed = sum(keep),
                            percent = round(100 * sum(keep) / n, 2)))

  out <- if (drop) data[keep, , drop = FALSE] else {
    data$.keep <- keep
    data$.reason <- factor(ifelse(keep, "kept", reason),
                           levels = c("kept", rules))
    data
  }
  rownames(out) <- NULL
  attr(out, "cleaning") <- tally
  attr(out, "n_input") <- n
  out
}
