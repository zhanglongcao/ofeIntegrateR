# Zoning a paddock into contiguous pseudo-environments from environmental
# covariates -- elevation, EM38, soil tests -- rather than from yield.
#
# The contiguity requirement is what rules out k-means. Clustering cells on
# their covariates alone puts a cell in the zone its soil resembles, wherever
# it happens to sit, so a zone arrives as confetti: a machine cannot drive it,
# a sampler cannot stratify by it, and a strip trial cannot use it as a
# pseudo-environment. The method here (Assuncao et al. 2006, SKATER) builds a
# minimum spanning tree over the neighbourhood graph and prunes it, so every
# zone is a subtree and therefore connected by construction.

#' Neighbourhood graph over the cells
#' @keywords internal
#' @noRd
.cell_neighbours <- function(xy, row = NULL, col = NULL,
                             type = c("auto", "rook", "queen", "knn"),
                             n_neighbours = 8L) {
  type <- match.arg(type)
  n <- nrow(xy)
  have_grid <- !is.null(row) && !is.null(col)
  if (type == "auto") type <- if (have_grid) "rook" else "knn"
  if (type %in% c("rook", "queen") && !have_grid) {
    stop("`neighbours = \"", type, "\"` needs `row` and `col` columns. Use ",
         "neighbours = \"knn\" for irregularly placed cells.", call. = FALSE)
  }

  if (type %in% c("rook", "queen")) {
    key <- paste(row, col, sep = "_")
    pos <- stats::setNames(seq_len(n), key)
    shifts <- if (type == "rook") {
      list(c(1, 0), c(0, 1))
    } else {
      list(c(1, 0), c(0, 1), c(1, 1), c(1, -1))
    }
    ed <- lapply(shifts, function(s) {
      j <- pos[paste(row + s[1], col + s[2], sep = "_")]
      i <- which(!is.na(j))
      if (length(i) == 0L) return(NULL)
      cbind(i, unname(j[i]))
    })
    e <- do.call(rbind, ed)
  } else {
    if (n > 4000L) {
      stop("`neighbours = \"knn\"` forms an ", n, " x ", n, " distance ",
           "matrix. Aggregate to a coarser lattice with grid_dense_layer(), ",
           "which also gives the row/col columns the grid neighbourhoods use.",
           call. = FALSE)
    }
    d <- as.matrix(stats::dist(xy))
    diag(d) <- Inf
    kk <- min(as.integer(n_neighbours), n - 1L)
    e <- do.call(rbind, lapply(seq_len(n), function(i) {
      cbind(i, order(d[i, ])[seq_len(kk)])
    }))
  }
  if (is.null(e)) stop("The cells have no neighbours.", call. = FALSE)

  # Symmetrise and drop duplicates: an undirected edge list.
  e <- cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]))
  e <- e[e[, 1] != e[, 2], , drop = FALSE]
  e <- unique(e)
  e
}

#' Union-find over n elements
#' @keywords internal
#' @noRd
.uf_new <- function(n) list(parent = seq_len(n), rank = integer(n))

#' @keywords internal
#' @noRd
.uf_find <- function(uf, x) {
  while (uf$parent[x] != x) x <- uf$parent[x]
  x
}

#' Components of an edge list
#' @keywords internal
#' @noRd
.components <- function(edges, n) {
  uf <- .uf_new(n)
  for (k in seq_len(nrow(edges))) {
    a <- .uf_find(uf, edges[k, 1])
    b <- .uf_find(uf, edges[k, 2])
    if (a != b) uf$parent[b] <- a
  }
  roots <- vapply(seq_len(n), function(i) .uf_find(uf, i), integer(1))
  as.integer(factor(roots))
}

#' Join a disconnected neighbourhood graph by its shortest spatial links
#'
#' A lattice with a hole through it, or a paddock in two blocks, gives a graph
#' in pieces. Refusing to zone it would be unhelpful and zoning the pieces
#' separately would be wrong, so the pieces are linked by their closest cells
#' and the tree is built over the whole.
#' @keywords internal
#' @noRd
.connect_components <- function(edges, xy) {
  n <- nrow(xy)
  comp <- .components(edges, n)
  while (length(unique(comp)) > 1L) {
    sizes <- table(comp)
    small <- as.integer(names(sizes)[which.min(sizes)])
    a <- which(comp == small)
    b <- which(comp != small)
    d <- outer(xy[a, 1], xy[b, 1], "-")^2 + outer(xy[a, 2], xy[b, 2], "-")^2
    hit <- which.min(d)
    ia <- a[(hit - 1L) %% length(a) + 1L]
    ib <- b[(hit - 1L) %/% length(a) + 1L]
    edges <- rbind(edges, c(min(ia, ib), max(ia, ib)))
    comp <- .components(edges, n)
  }
  edges
}

#' Minimum spanning tree by Kruskal's algorithm
#' @keywords internal
#' @noRd
.mst <- function(edges, w, n) {
  o <- order(w)
  uf <- .uf_new(n)
  keep <- integer(0)
  for (k in o) {
    a <- .uf_find(uf, edges[k, 1])
    b <- .uf_find(uf, edges[k, 2])
    if (a != b) {
      if (uf$rank[a] < uf$rank[b]) {
        tmp <- a; a <- b; b <- tmp
      }
      uf$parent[b] <- a
      if (uf$rank[a] == uf$rank[b]) uf$rank[a] <- uf$rank[a] + 1L
      keep <- c(keep, k)
      if (length(keep) == n - 1L) break
    }
  }
  edges[keep, , drop = FALSE]
}

#' Within-cluster sum of squared deviations
#' @keywords internal
#' @noRd
.ssd <- function(idx, X) {
  if (length(idx) <= 1L) return(0)
  Xi <- X[idx, , drop = FALSE]
  sum((Xi - rep(colMeans(Xi), each = nrow(Xi)))^2)
}

#' Best edge to cut in one subtree, by the reduction in within-zone variance
#'
#' Subtree sums are accumulated once in a single reverse-breadth-first pass, so
#' every candidate edge is then evaluated in constant time. Testing each edge
#' by re-traversing the tree would be quadratic, which is what makes naive
#' SKATER implementations unusable on a full lattice.
#' @keywords internal
#' @noRd
.best_cut <- function(verts, tree_edges, X, min_cells) {
  nv <- length(verts)
  if (nv < 2L * min_cells) return(NULL)
  local <- match(seq_len(max(verts)), verts)          # global -> local index
  e <- cbind(local[tree_edges[, 1]], local[tree_edges[, 2]])
  adj <- vector("list", nv)
  for (k in seq_len(nrow(e))) {
    adj[[e[k, 1]]] <- c(adj[[e[k, 1]]], e[k, 2])
    adj[[e[k, 2]]] <- c(adj[[e[k, 2]]], e[k, 1])
  }

  # Breadth-first order from an arbitrary root, then accumulate up it.
  parent <- integer(nv)
  order_v <- integer(nv)
  seen <- logical(nv)
  head_i <- 1L; tail_i <- 1L
  order_v[1L] <- 1L; seen[1L] <- TRUE; parent[1L] <- 0L
  while (head_i <= tail_i) {
    v <- order_v[head_i]; head_i <- head_i + 1L
    for (u in adj[[v]]) {
      if (!seen[u]) {
        seen[u] <- TRUE
        parent[u] <- v
        tail_i <- tail_i + 1L
        order_v[tail_i] <- u
      }
    }
  }

  Xi <- X[verts, , drop = FALSE]
  p <- ncol(Xi)
  acc_n <- rep(1L, nv)
  acc_S <- Xi
  acc_Q <- Xi^2
  for (v in rev(order_v)) {
    pa <- parent[v]
    if (pa > 0L) {
      acc_n[pa] <- acc_n[pa] + acc_n[v]
      acc_S[pa, ] <- acc_S[pa, ] + acc_S[v, ]
      acc_Q[pa, ] <- acc_Q[pa, ] + acc_Q[v, ]
    }
  }

  N <- nv
  S <- acc_S[1L, ]
  Q <- acc_Q[1L, ]
  total <- sum(Q - S^2 / N)

  cand <- which(parent > 0L)
  n1 <- acc_n[cand]
  ok <- n1 >= min_cells & (N - n1) >= min_cells
  cand <- cand[ok]
  if (length(cand) == 0L) return(NULL)
  n1 <- acc_n[cand]
  S1 <- acc_S[cand, , drop = FALSE]
  Q1 <- acc_Q[cand, , drop = FALSE]
  S2 <- rep(S, each = length(cand)) - S1
  Q2 <- rep(Q, each = length(cand)) - Q1
  n2 <- N - n1
  ssd1 <- rowSums(Q1 - S1^2 / n1)
  ssd2 <- rowSums(Q2 - S2^2 / n2)
  gain <- total - ssd1 - ssd2

  best <- which.max(gain)
  v <- cand[best]
  # The cut edge is (v, parent[v]); the new zone is the subtree below v.
  sub <- integer(0)
  stack <- v
  while (length(stack) > 0L) {
    cur <- stack[length(stack)]
    stack <- stack[-length(stack)]
    sub <- c(sub, cur)
    kids <- adj[[cur]][parent[adj[[cur]]] == cur]
    stack <- c(stack, kids)
  }
  list(gain = gain[best],
       part1 = verts[sub],
       part2 = verts[setdiff(seq_len(nv), sub)],
       cut = c(verts[v], verts[parent[v]]))
}

#' Best axis-aligned cut of one rectangle
#'
#' Considers every horizontal and every vertical line that separates two
#' adjacent coordinate positions, and returns the one removing most within-zone
#' variance. Cutting a rectangle all the way across always leaves two
#' rectangles, so recursing on the result gives rectangular zones however many
#' times it is applied -- which is the property `method = "skater"` cannot
#' offer, since a pruned spanning tree is connected but any shape.
#'
#' Sums are accumulated once per coordinate line and then cumulated, so each
#' candidate cut costs O(p) rather than a pass over the cells.
#' @keywords internal
#' @noRd
.best_rect_cut <- function(v, ax, ay, X, min_cells) {
  Xv <- X[v, , drop = FALSE]
  N <- length(v)
  S <- colSums(Xv)
  Q <- colSums(Xv^2)
  total <- sum(Q - S^2 / N)
  best <- NULL

  for (axis in 1:2) {
    a <- if (axis == 1L) ax[v] else ay[v]
    ua <- sort(unique(a))
    if (length(ua) < 2L) next
    g <- match(a, ua)
    cn <- cumsum(tabulate(g, nbins = length(ua)))
    S_g <- rowsum(Xv, g, reorder = TRUE)
    Q_g <- rowsum(Xv^2, g, reorder = TRUE)
    cS <- matrix(apply(S_g, 2, cumsum), ncol = ncol(S_g))
    cQ <- matrix(apply(Q_g, 2, cumsum), ncol = ncol(Q_g))

    j <- seq_len(length(ua) - 1L)
    n1 <- cn[j]
    n2 <- N - n1
    ok <- n1 >= min_cells & n2 >= min_cells
    if (!any(ok)) next
    j <- j[ok]; n1 <- n1[ok]; n2 <- n2[ok]
    S1 <- cS[j, , drop = FALSE]
    Q1 <- cQ[j, , drop = FALSE]
    S2 <- rep(S, each = length(j)) - S1
    Q2 <- rep(Q, each = length(j)) - Q1
    gain <- total - rowSums(Q1 - S1^2 / n1) - rowSums(Q2 - S2^2 / n2)

    b <- which.max(gain)
    if (is.null(best) || gain[b] > best$gain) {
      at <- ua[j[b]]
      lower <- a <= at
      best <- list(gain = gain[b], axis = axis, at = at,
                   part1 = v[lower], part2 = v[!lower])
    }
  }
  best
}

#' Omnidirectional practical range of a multivariate covariate field
#' @keywords internal
#' @noRd
.omni_range <- function(xy, v, n_bin = 15L, max_pts = 1200L) {
  n <- nrow(xy)
  if (n > max_pts) {
    i <- sample.int(n, max_pts)
    xy <- xy[i, , drop = FALSE]
    v <- v[i]
  }
  d <- stats::dist(xy)
  g <- stats::dist(v)^2 / 2
  cutoff <- max(d) / 3
  keep <- d > 0 & d <= cutoff
  if (sum(keep) < 30L) return(NA_real_)
  br <- seq(0, cutoff, length.out = n_bin + 1L)
  bin <- cut(as.numeric(d)[keep], br, include.lowest = TRUE)
  gamma <- tapply(as.numeric(g)[keep], bin, mean)
  n_pair <- tapply(as.numeric(g)[keep], bin, length)
  h <- (br[-1L] + br[-length(br)]) / 2
  ok <- !is.na(gamma)
  a <- .fit_exp_variogram(h[ok], as.numeric(gamma[ok]),
                          as.numeric(n_pair[ok]))
  if (is.na(a)) return(NA_real_)
  3 * a
}

#' Partition a paddock into contiguous zones from environmental covariates
#'
#' Divides a paddock into pseudo-environments using covariates that are known
#' before harvest -- elevation, EM38 or gamma survey, soil test grids, a
#' previous season's imagery -- rather than from the yield being analysed. The
#' zones are guaranteed to be spatially contiguous.
#'
#' # Shape of the zones
#'
#' Three methods, differing only in what shape a zone is allowed to be. The
#' choice is a practical one, not a statistical one: a zone has to be something
#' a machine can drive and a sampling plan can stratify by.
#'
#' \describe{
#'   \item{`"rectangle"` (default)}{Rectangular blocks. The paddock is cut
#'     repeatedly by lines running the full width or full length of the region
#'     being split, each cut placed where it removes most within-zone variance.
#'     Cutting a rectangle across always leaves two rectangles, so every zone is
#'     a rectangle however many cuts are made. This is the shape a variable-rate
#'     prescription, a soil sampling grid or a set of management blocks wants,
#'     and it is the two-dimensional version of what [partition_pseudo_env()]
#'     does along a single axis. Where the paddock is not itself a rectangle --
#'     a clipped corner, a headland left out -- a zone is a rectangle
#'     intersected with the paddock, which is still a block you can drive.}
#'   \item{`"skater"`}{Contiguous zones of any shape. A minimum spanning tree
#'     is built over the neighbourhood graph -- edges weighted by distance
#'     between cells in standardised covariate space -- and pruned one edge at a
#'     time, always cutting the edge that removes most within-zone variance.
#'     Every zone is a subtree of a connected graph, so every zone is connected,
#'     but a zone may send a finger out between two others. Use it when the
#'     zones are there to describe the soil faithfully rather than to be
#'     worked, and when a boundary genuinely runs at an angle or along a
#'     contour. Being free to follow the boundary does not guarantee it fits
#'     better: both methods are greedy, and a greedy tree pruning can strand
#'     itself where a well-placed straight cut would not, so rectangles
#'     sometimes explain more of the covariate variance at the same `k`.
#'     Compare the `r2` column of the returned table on your own paddock rather
#'     than assuming. The method is SKATER (Assuncao et al. 2006).}
#'   \item{`"kmeans"`}{Plain k-means on the covariates, for comparison. It is
#'     not contiguous at all: cells are assigned to the zone their soil
#'     resembles wherever they happen to sit, so zones come back as confetti
#'     scattered across the paddock. The `patches` column of the zone summary
#'     counts the connected pieces of each zone -- 1 for every rectangle and
#'     every skater zone, more than 1 whenever k-means fragments.}
#' }
#'
#' What the choice costs depends entirely on the shape of the underlying
#' feature. On simulated paddocks with one clear feature in one covariate, the
#' variance explained at `k = 2` was:
#'
#' \tabular{lrr}{
#'   **feature**                  \tab **rectangle** \tab **skater** \cr
#'   block boundaries on the axes \tab 0.78 \tab 0.78 \cr
#'   a smooth gradient up the paddock \tab 0.69 \tab 0.67 \cr
#'   a boundary running diagonally \tab 0.25 \tab 0.98 \cr
#'   a round patch in the middle  \tab 0.13 \tab 0.97 \cr
#' }
#'
#' So rectangles cost nothing when the structure is already blocky or a
#' gradient, and cost almost everything when the boundary runs at an angle or
#' curves -- a creek line, a dune, a contour. If your covariate map has a
#' feature like that, either use `"skater"` and accept the awkward shapes, or
#' raise `k` so rectangles can approximate the boundary in steps.
#'
#' # Choosing k
#'
#' With `k = NULL`, zones are added while each new one earns its keep: the
#' default `criterion = "gain"` stops at the first `k` whose successor would
#' explain less than `min_gain` (5% by default) more of the covariate variance.
#' This is the elbow rule made explicit, and it is the one an agronomist can
#' argue with, because zones cost management effort and the threshold is where
#' that trade-off is stated.
#'
#' `criterion = "ch"` maximises the Calinski-Harabasz index instead. Be aware
#' that on a smoothly varying covariate -- a slope, a gradual texture change --
#' CH often has no interior maximum at all and simply picks `k_max`, because
#' every further split of a smooth surface really does separate it further.
#' That is a property of smooth fields rather than a failure of the index, and
#' it is why it is not the default here.
#'
#' Either way, if the partition explains less than `min_r2` of the covariate
#' variance, one zone is returned: a paddock uniform in its covariates should
#' be reported as uniform. The full criterion table is returned so `k` can be
#' overridden knowingly.
#'
#' The practical range of an exponential variogram fitted to the covariate
#' field is reported alongside, because it tells you the scale at which the
#' covariates vary and hence how far apart to sample. It deliberately does
#' **not** constrain `k`: a real zone boundary is a step, a step keeps the
#' variogram climbing without reaching a sill, and the fitted range is then
#' inflated by the very structure being looked for. Read it as context, not as
#' a limit.
#'
#' # Using the zones in a trial analysis
#'
#' Pass `treat` and the zone summary gains a `treatments` column counting the
#' treatment levels present in each zone. A zone missing a treatment cannot
#' support a `zone:treat` term, and you are warned rather than left to discover
#' it when the fit drops rank. Zones from this function are arbitrary shapes,
#' so a strip trial will often want [partition_pseudo_env()] instead, which
#' cuts across the strips and so keeps every treatment in every zone.
#'
#' @param data Data frame, one row per cell (e.g. from [grid_dense_layer()]).
#' @param covariates Character vector naming the columns of `data` to zone on,
#'   such as `c("elevation", "clay")`. Required, and there is deliberately no
#'   default: which layers define a zone is an agronomic decision. All must be
#'   numeric; rows missing any of them are dropped.
#' @param k Integer or `NULL`; the number of zones. `NULL` chooses it (see
#'   above).
#' @param row,col Character or `NULL`; integer lattice indices, used to build
#'   the grid neighbourhoods. Defaults to `"row"`/`"col"` when `data` has them.
#' @param x,y Character or `NULL`; projected coordinates in metres. Defaults to
#'   `"x_centre"`/`"y_centre"` when present, else `"x"`/`"y"`. Needed for the
#'   spatial range and for `neighbours = "knn"`.
#' @param method `"rectangle"` (default) for rectangular blocks,
#'   `"skater"` for contiguous zones of any shape, or `"kmeans"` for
#'   unconstrained clustering (not contiguous; for comparison). See above.
#' @param neighbours `"auto"` (grid neighbourhoods when `row`/`col` are
#'   available, otherwise nearest neighbours), `"rook"`, `"queen"` or `"knn"`.
#' @param n_neighbours Integer; neighbours per cell when `neighbours = "knn"`.
#' @param k_max Integer; largest number of zones considered when `k` is
#'   chosen automatically.
#' @param criterion `"gain"` (default) to stop when another zone would explain
#'   less than `min_gain` more of the covariate variance, or `"ch"` to maximise
#'   the Calinski-Harabasz index. See above.
#' @param min_gain Numeric; the share of covariate variance an extra zone must
#'   explain to be kept, under `criterion = "gain"`.
#' @param min_cells Integer or `NULL`; smallest zone allowed. Defaults to 2% of
#'   the cells, with a floor of 5.
#' @param min_r2 Numeric; if the chosen partition explains less than this share
#'   of covariate variance, one zone is returned instead.
#' @param scale Logical; standardise each covariate to mean 0 and standard
#'   deviation 1 before clustering (default `TRUE`). Leave it on unless the
#'   covariates are already on a common scale -- otherwise elevation in metres
#'   will silently outrank pH.
#' @param n_pc Integer or `NULL`; cluster on the first `n_pc` principal
#'   components instead of the covariates themselves. Useful when several
#'   layers measure much the same thing.
#' @param treat Character or `NULL`; a treatment column, used only to report
#'   and warn about treatment coverage within zones.
#' @param name Character; name of the zone column added to `data`.
#'
#' @return `data` with an added factor column (named by `name`) giving the zone
#'   of each cell; cells dropped for missing covariates get `NA`. The
#'   derivation is attached as `attr(, "partition")`: the `method`, `k`,
#'   `covariates`, the fitted `range`, the `table` of criteria over k (`ssd`,
#'   `r2`, the marginal `gain`, and `ch`), and a `zones` summary with each
#'   zone's cell count, area, number of connected `patches`, corner
#'   coordinates, covariate means and (with `treat`) treatment coverage.
#'
#' @references
#' Assuncao, R.M., Neves, M.C., Camara, G. and da Costa Freitas, C. (2006)
#' Efficient regionalization techniques for socio-economic geographical units
#' using minimum spanning trees. *International Journal of Geographical
#' Information Science* **20**, 797-811.
#'
#' @seealso [partition_pseudo_env()] for slicing a strip trial across its
#'   length, and [adaptive_residual()] to turn either set of zones into a
#'   residual structure.
#'
#' @examples
#' # A paddock with a ridge running through it and a sandy corner
#' g <- expand.grid(row = 1:24, col = 1:18)
#' g$x <- g$col * 10
#' g$y <- g$row * 10
#' g$elevation <- 100 + 6 * exp(-((g$y - 120)^2) / 2000)
#' g$clay <- ifelse(g$x > 120 & g$y > 140, 18, 32)
#' set.seed(1)
#' g$elevation <- g$elevation + stats::rnorm(nrow(g), 0, 0.2)
#' g$clay <- g$clay + stats::rnorm(nrow(g), 0, 1)
#'
#' z <- partition_paddock(g, covariates = c("elevation", "clay"))
#' attr(z, "partition")$zones
#' table(z$zone)
#'
#' @export
partition_paddock <- function(data, covariates, k = NULL,
                              row = NULL, col = NULL, x = NULL, y = NULL,
                              method = c("rectangle", "skater", "kmeans"),
                              neighbours = c("auto", "rook", "queen", "knn"),
                              n_neighbours = 8L, k_max = 8L,
                              criterion = c("gain", "ch"), min_gain = 0.05,
                              min_cells = NULL, min_r2 = 0.10, scale = TRUE,
                              n_pc = NULL, treat = NULL, name = "zone") {
  method <- match.arg(method)
  neighbours <- match.arg(neighbours)
  criterion <- match.arg(criterion)
  data <- as.data.frame(data)

  if (missing(covariates) || is.null(covariates) ||
      length(covariates) == 0L || all(is.na(covariates))) {
    stop("`covariates` is required: name the columns of `data` to zone on, ",
         "e.g. partition_paddock(dat, covariates = c(\"elevation\", ",
         "\"soil\")). There is no sensible default -- which layers define a ",
         "zone is the agronomy, not something this function can guess.",
         call. = FALSE)
  }
  if (!is.character(covariates)) {
    stop("`covariates` must be column names, as character, not ",
         class(covariates)[1], ".", call. = FALSE)
  }
  dup <- unique(covariates[duplicated(covariates)])
  if (length(dup) > 0L) {
    stop("`covariates` names the same column twice: ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(row) && "row" %in% names(data)) row <- "row"
  if (is.null(col) && "col" %in% names(data)) col <- "col"
  if (is.null(x)) x <- if ("x_centre" %in% names(data)) "x_centre" else
    if ("x" %in% names(data)) "x" else NULL
  if (is.null(y)) y <- if ("y_centre" %in% names(data)) "y_centre" else
    if ("y" %in% names(data)) "y" else NULL
  if (is.null(x) || is.null(y)) {
    if (is.null(row) || is.null(col)) {
      stop("Need coordinates: either `x`/`y`, or `row`/`col` lattice indices.",
           call. = FALSE)
    }
    x <- row
    y <- col
  }

  needed <- c(covariates, x, y, row, col, treat)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  num_bad <- covariates[!vapply(data[covariates], is.numeric, logical(1))]
  if (length(num_bad) > 0L) {
    stop("`covariates` must be numeric; these are not: ",
         paste(num_bad, collapse = ", "), ".", call. = FALSE)
  }

  keep <- stats::complete.cases(data[, c(covariates, x, y), drop = FALSE])
  n_drop <- sum(!keep)
  idx_keep <- which(keep)
  n <- length(idx_keep)
  if (n < 4L) stop("Fewer than four cells have complete covariates.",
                   call. = FALSE)

  X <- as.matrix(data[idx_keep, covariates, drop = FALSE])
  if (scale) {
    sdv <- apply(X, 2, stats::sd)
    flat <- sdv == 0 | !is.finite(sdv)
    if (all(flat)) {
      stop("Every covariate is constant; there is nothing to partition.",
           call. = FALSE)
    }
    if (any(flat)) {
      warning("Constant covariate(s) dropped: ",
              paste(covariates[flat], collapse = ", "), ".", call. = FALSE)
      X <- X[, !flat, drop = FALSE]
      sdv <- sdv[!flat]
    }
    X <- (X - rep(colMeans(X), each = n)) / rep(sdv, each = n)
  }
  if (!is.null(n_pc)) {
    n_pc <- min(as.integer(n_pc), ncol(X))
    X <- stats::prcomp(X, center = TRUE, scale. = FALSE)$x[, seq_len(n_pc),
                                                           drop = FALSE]
  }

  xy <- cbind(as.numeric(data[[x]][idx_keep]), as.numeric(data[[y]][idx_keep]))
  rowv <- if (is.null(row)) NULL else as.integer(data[[row]][idx_keep])
  colv <- if (is.null(col)) NULL else as.integer(data[[col]][idx_keep])

  edges <- .cell_neighbours(xy, rowv, colv, neighbours, n_neighbours)
  edges <- .connect_components(edges, xy)

  # Cell area from the spacing of the neighbourhood links, for the range guard.
  link_d <- sqrt((xy[edges[, 1], 1] - xy[edges[, 2], 1])^2 +
                   (xy[edges[, 1], 2] - xy[edges[, 2], 2])^2)
  spacing <- stats::median(link_d[link_d > 0])
  if (!is.finite(spacing) || spacing <= 0) spacing <- 1
  area_total <- n * spacing^2

  rng <- .omni_range(xy, X[, 1])

  if (is.null(min_cells)) min_cells <- max(5L, ceiling(0.02 * n))
  min_cells <- max(1L, as.integer(min_cells))
  k_allow <- max(1L, min(as.integer(k_max), floor(n / min_cells)))
  k_search <- if (is.null(k)) k_allow else max(as.integer(k), 1L)
  if (!is.null(k) && k > floor(n / min_cells)) {
    stop("`k = ", k, "` needs at least ", k * min_cells, " cells at ",
         "`min_cells = ", min_cells, "`, but only ", n, " are available.",
         call. = FALSE)
  }

  ssd_total <- .ssd(seq_len(n), X)
  if (method == "rectangle") {
    ax <- match(xy[, 1], sort(unique(xy[, 1])))
    ay <- match(xy[, 2], sort(unique(xy[, 2])))
    rects <- list(seq_len(n))
    seq_ssd <- ssd_total
    labels_at_k <- list(rep(1L, n))
    for (step in seq_len(max(k_search, 1L) - 1L)) {
      cuts <- lapply(rects, function(v) .best_rect_cut(v, ax, ay, X, min_cells))
      gains <- vapply(cuts, function(z) if (is.null(z)) -Inf else z$gain,
                      numeric(1))
      if (all(!is.finite(gains))) break
      ci <- which.max(gains)
      cut <- cuts[[ci]]
      rects[[ci]] <- cut$part1
      rects[[length(rects) + 1L]] <- cut$part2
      lab <- integer(n)
      for (j in seq_along(rects)) lab[rects[[j]]] <- j
      labels_at_k[[length(labels_at_k) + 1L]] <- lab
      seq_ssd <- c(seq_ssd,
                   sum(vapply(rects, function(v) .ssd(v, X), numeric(1))))
    }
  } else if (method == "skater") {
    w <- sqrt(rowSums((X[edges[, 1], , drop = FALSE] -
                         X[edges[, 2], , drop = FALSE])^2))
    tree <- .mst(edges, w, n)
    tree_w <- sqrt(rowSums((X[tree[, 1], , drop = FALSE] -
                              X[tree[, 2], , drop = FALSE])^2))

    clusters <- list(seq_len(n))
    cl_edges <- list(tree)
    seq_ssd <- ssd_total
    labels_at_k <- list(rep(1L, n))
    for (step in seq_len(max(k_search, 1L) - 1L)) {
      cuts <- lapply(seq_along(clusters), function(ci) {
        .best_cut(clusters[[ci]], cl_edges[[ci]], X, min_cells)
      })
      gains <- vapply(cuts, function(z) if (is.null(z)) -Inf else z$gain,
                      numeric(1))
      if (all(!is.finite(gains))) break
      ci <- which.max(gains)
      cut <- cuts[[ci]]
      e <- cl_edges[[ci]]
      in1 <- e[, 1] %in% cut$part1 & e[, 2] %in% cut$part1
      in2 <- e[, 1] %in% cut$part2 & e[, 2] %in% cut$part2
      clusters[[ci]] <- cut$part1
      cl_edges[[ci]] <- e[in1, , drop = FALSE]
      clusters[[length(clusters) + 1L]] <- cut$part2
      cl_edges[[length(cl_edges) + 1L]] <- e[in2, , drop = FALSE]

      lab <- integer(n)
      for (j in seq_along(clusters)) lab[clusters[[j]]] <- j
      labels_at_k[[length(labels_at_k) + 1L]] <- lab
      # NB anonymous wrapper: passing `X = X` through vapply's dots would be
      # matched to vapply's own `X` argument instead of .ssd's.
      seq_ssd <- c(seq_ssd,
                   sum(vapply(clusters, function(v) .ssd(v, X), numeric(1))))
    }
  } else {
    labels_at_k <- list(rep(1L, n))
    seq_ssd <- ssd_total
    for (kk in seq_len(max(k_search, 1L))[-1L]) {
      km <- stats::kmeans(X, centers = kk, nstart = 25)
      labels_at_k[[kk]] <- km$cluster
      seq_ssd <- c(seq_ssd, sum(km$withinss))
    }
  }

  k_reached <- length(labels_at_k)
  r2 <- 1 - seq_ssd / ssd_total
  crit <- data.frame(
    k = seq_len(k_reached),
    ssd = seq_ssd,
    r2 = r2,
    gain = c(NA_real_, diff(r2)),
    ch = c(NA_real_, vapply(seq_len(k_reached)[-1L], function(kk) {
      ((ssd_total - seq_ssd[kk]) / (kk - 1)) / (seq_ssd[kk] / (n - kk))
    }, numeric(1)))
  )

  if (is.null(k)) {
    k_use <- if (k_reached < 2L) 1L else if (criterion == "ch") {
      which.max(crit$ch[-1L]) + 1L
    } else {
      # Keep adding zones while the next one explains at least `min_gain` more
      # of the covariate variance; stop at the first that does not.
      worth <- crit$gain[-1L] >= min_gain
      if (!worth[1L]) 1L else {
        stop_at <- which(!worth)
        if (length(stop_at) == 0L) k_reached else stop_at[1L]
      }
    }
    if (crit$r2[k_use] < min_r2) k_use <- 1L
  } else {
    k_use <- min(as.integer(k), k_reached)
    if (k_use < k) {
      warning("Only ", k_use, " zone(s) could be separated while keeping ",
              "each at least ", min_cells, " cells.", call. = FALSE)
    }
  }
  lab <- labels_at_k[[k_use]]

  # Order zones by their first covariate's mean, so zone 1 is the low end and
  # the numbering means something across re-runs.
  ord <- order(tapply(X[, 1], lab, mean))
  lab <- match(lab, ord)

  out_zone <- rep(NA_integer_, nrow(data))
  out_zone[idx_keep] <- lab
  data[[name]] <- factor(out_zone, levels = seq_len(k_use))

  patches <- vapply(seq_len(k_use), function(j) {
    v <- which(lab == j)
    e <- edges[edges[, 1] %in% v & edges[, 2] %in% v, , drop = FALSE]
    if (length(v) == 1L) return(1L)
    length(unique(.components(cbind(match(e[, 1], v), match(e[, 2], v)),
                              length(v))))
  }, integer(1))

  summ <- data.frame(zone = factor(seq_len(k_use)),
                     n = as.integer(table(factor(lab, levels = seq_len(k_use)))),
                     stringsAsFactors = FALSE)
  summ$area <- summ$n * spacing^2
  summ$patches <- patches
  # The corners of each zone, so a block can actually be marked out. For
  # `method = "rectangle"` this is the zone itself; for the others it is the
  # bounding box of a shape that may not fill it.
  for (nm in c("x_min", "x_max", "y_min", "y_max")) summ[[nm]] <- NA_real_
  for (j in seq_len(k_use)) {
    v <- which(lab == j)
    summ$x_min[j] <- min(xy[v, 1]); summ$x_max[j] <- max(xy[v, 1])
    summ$y_min[j] <- min(xy[v, 2]); summ$y_max[j] <- max(xy[v, 2])
  }
  cov_used <- colnames(X)
  for (cv in covariates) {
    summ[[cv]] <- as.numeric(tapply(data[[cv]][idx_keep], lab, mean))
  }
  if (!is.null(treat)) {
    tv <- factor(data[[treat]][idx_keep])
    summ$treatments <- as.integer(tapply(tv, lab, function(z)
      length(unique(z[!is.na(z)]))))
    short <- summ$zone[summ$treatments < nlevels(tv)]
    if (length(short) > 0L) {
      warning("Zone(s) ", paste(short, collapse = ", "), " do not contain ",
              "every level of `", treat, "`, so a ", name, ":", treat,
              " term is not estimable there. partition_pseudo_env() cuts ",
              "across the strips and keeps every treatment in every zone.",
              call. = FALSE)
    }
  }

  attr(data, "partition") <- list(
    method = method, k = k_use, covariates = covariates,
    clustered_on = cov_used, scaled = scale, n_pc = n_pc,
    range = rng, criterion = criterion, min_gain = min_gain,
    min_cells = min_cells, min_r2 = min_r2,
    spacing = spacing, area = area_total, n_cells = n, n_dropped = n_drop,
    neighbours = neighbours, table = crit, zones = summ
  )
  data
}
