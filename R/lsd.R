# Compact letter displays: turning a set of pairwise comparisons into the
# a/b/c annotation a trial report prints next to its treatment means.

#' A supply of display letters
#'
#' a-z, then A-Z, then aa, bb, ... A display needing more than 52 letters is
#' already unreadable, but failing there would be worse than being ugly.
#' @keywords internal
#' @noRd
.letter_pool <- function(n) {
  base <- c(letters, LETTERS)
  if (n <= length(base)) return(base[seq_len(n)])
  reps <- ceiling(n / length(base))
  unlist(lapply(seq_len(reps), function(k) strrep(base, k)))[seq_len(n)]
}

#' All maximal cliques of an undirected graph (Bron-Kerbosch with pivoting)
#'
#' A compact letter display is exactly the set of maximal cliques of the
#' "not significantly different" graph: one letter per clique, shared by every
#' level in it. Enumerating them directly gives the minimal correct display,
#' where the sweep-and-absorb algorithms used elsewhere can leave redundant
#' letters behind.
#' @keywords internal
#' @noRd
.maximal_cliques <- function(adj, max_cliques = 10000L) {
  n <- nrow(adj)
  if (n == 0L) return(list())
  diag(adj) <- FALSE
  nb <- lapply(seq_len(n), function(i) which(adj[i, ]))
  res <- list()
  overflow <- FALSE

  bk <- function(R, P, X) {
    if (overflow) return(invisible(NULL))
    if (length(P) == 0L && length(X) == 0L) {
      res[[length(res) + 1L]] <<- R
      if (length(res) > max_cliques) overflow <<- TRUE
      return(invisible(NULL))
    }
    PX <- c(P, X)
    piv <- PX[which.max(vapply(PX, function(u) length(intersect(nb[[u]], P)),
                               integer(1)))]
    for (v in setdiff(P, nb[[piv]])) {
      bk(c(R, v), intersect(P, nb[[v]]), intersect(X, nb[[v]]))
      P <- setdiff(P, v)
      X <- c(X, v)
    }
    invisible(NULL)
  }

  bk(integer(0), seq_len(n), integer(0))
  if (overflow) return(NULL)
  res
}

#' Compact letter display from any set of pairwise comparisons
#'
#' Turns pairwise p-values into the a/b/c annotation: two levels share a letter
#' when they are not significantly different. Levels are ordered by their means
#' where these are supplied, so `"a"` marks the best treatment.
#'
#' This is deliberately separate from [ofe_lsd()] and takes only the
#' comparisons, so it can letter the output of anything -- `asreml::predict()`,
#' `emmeans`, a table typed in by hand -- and not only a fit from this package.
#'
#' @param comparisons Either a symmetric matrix of p-values whose dimnames are
#'   the level names, or a data frame with a p-value column and the two levels
#'   of each comparison. The level columns are taken from `level1` and `level2`
#'   if present, otherwise from a `contrast` column of the form `"B - A"`,
#'   otherwise from the first two non-numeric columns.
#' @param alpha Significance level; comparisons with `p.value > alpha` are
#'   treated as "not different" and so share a letter.
#' @param means Optional named numeric vector of the means, used to order the
#'   levels. Without it the levels keep the order they appear in.
#' @param decreasing Logical; with `means` supplied, order from the largest
#'   mean down, so `"a"` goes to the highest. Set `FALSE` if smaller is better.
#' @param p_col Character; name of the p-value column when `comparisons` is a
#'   data frame.
#'
#' @return A named character vector of letter groups, one per level, in the
#'   order the letters were assigned.
#'
#' @examples
#' cmp <- data.frame(
#'   level1 = c("A", "A", "B"),
#'   level2 = c("B", "C", "C"),
#'   p.value = c(0.40, 0.001, 0.002)
#' )
#' compact_letters(cmp, means = c(A = 3.0, B = 3.1, C = 5.0))
#'
#' @seealso [ofe_lsd()] to go straight from a fitted model to a lettered table.
#' @export
compact_letters <- function(comparisons, alpha = 0.05, means = NULL,
                            decreasing = TRUE, p_col = "p.value") {
  if (is.matrix(comparisons)) {
    lv <- rownames(comparisons)
    if (is.null(lv)) lv <- colnames(comparisons)
    if (is.null(lv)) {
      stop("A p-value matrix needs dimnames giving the level names.",
           call. = FALSE)
    }
    pm <- comparisons
    dimnames(pm) <- list(lv, lv)
  } else {
    cmp <- as.data.frame(comparisons)
    if (!p_col %in% names(cmp)) {
      stop("`comparisons` has no column `", p_col, "`.", call. = FALSE)
    }
    if (all(c("level1", "level2") %in% names(cmp))) {
      l1 <- as.character(cmp$level1)
      l2 <- as.character(cmp$level2)
    } else if ("contrast" %in% names(cmp)) {
      parts <- strsplit(as.character(cmp$contrast), " - ", fixed = TRUE)
      if (any(lengths(parts) != 2L)) {
        stop("`contrast` must be of the form \"B - A\"; supply `level1` and ",
             "`level2` columns instead.", call. = FALSE)
      }
      l1 <- vapply(parts, `[`, character(1), 2L)
      l2 <- vapply(parts, `[`, character(1), 1L)
    } else {
      chr <- names(cmp)[!vapply(cmp, is.numeric, logical(1))]
      if (length(chr) < 2L) {
        stop("Cannot find the two level columns in `comparisons`. Supply ",
             "`level1` and `level2`.", call. = FALSE)
      }
      l1 <- as.character(cmp[[chr[1]]])
      l2 <- as.character(cmp[[chr[2]]])
    }
    lv <- unique(c(l1, l2))
    pm <- matrix(NA_real_, length(lv), length(lv), dimnames = list(lv, lv))
    idx <- cbind(match(l1, lv), match(l2, lv))
    pm[idx] <- cmp[[p_col]]
    pm[idx[, 2:1, drop = FALSE]] <- cmp[[p_col]]
  }

  if (!is.null(means)) {
    if (is.null(names(means))) {
      stop("`means` must be named with the level names.", call. = FALSE)
    }
    missing_lv <- setdiff(lv, names(means))
    if (length(missing_lv) > 0L) {
      stop("`means` is missing level(s): ",
           paste(missing_lv, collapse = ", "), ".", call. = FALSE)
    }
    ord <- order(means[lv], decreasing = decreasing)
    lv <- lv[ord]
    pm <- pm[lv, lv, drop = FALSE]
  }

  # Missing comparisons would silently become "different"; a level that was
  # never compared belongs in no-one's group, which is not what a reader takes
  # from a blank. Treat an absent p-value as not-different and say so.
  if (anyNA(pm[upper.tri(pm)])) {
    warning("Some pairs have no p-value; they are treated as not different.",
            call. = FALSE)
  }
  ns <- is.na(pm) | pm > alpha
  diag(ns) <- TRUE

  cliques <- .maximal_cliques(ns)
  if (is.null(cliques)) {
    stop("This comparison structure needs an unreasonable number of letters. ",
         "Use a multiplicity adjustment, or fewer levels.", call. = FALSE)
  }
  # Order the letters so the display reads as a staircase down the means.
  first <- vapply(cliques, min, integer(1))
  size <- vapply(cliques, length, integer(1))
  cliques <- cliques[order(first, -size)]

  pool <- .letter_pool(length(cliques))
  out <- vapply(seq_along(lv), function(i) {
    paste0(pool[vapply(cliques, function(cl) i %in% cl, logical(1))],
           collapse = "")
  }, character(1))
  stats::setNames(out, lv)
}

#' Predicted means with LSD letters
#'
#' The table a trial report prints: each level's predicted mean, its standard
#' error, and the a/b/c letters saying which means are separable. Levels
#' sharing a letter are not significantly different at `alpha`.
#'
#' Two ways of deciding "different" are offered, and they part company on
#' exactly the data this package is for.
#'
#' \describe{
#'   \item{`use = "pairwise"` (default)}{Each pair is judged on its own
#'     standard error of difference. This is the right choice for a spatial
#'     model of a strip trial: neighbouring strips are compared more precisely
#'     than distant ones, and the SEDs genuinely differ.}
#'   \item{`use = "lsd"`}{One least significant difference,
#'     \eqn{t_{1-\alpha/2, \nu} \times \overline{SED}}, applied to every
#'     comparison. This is the textbook LSD and what a balanced design implies;
#'     it is also what most published tables mean by "LSD (5%)". Quote the
#'     value with the table.}
#' }
#'
#' The default `adjust = "none"` gives unprotected LSD comparisons, which is
#' the convention these letters normally carry -- and which does not control
#' the family-wise error rate. With more than three or four treatments, say so
#' in the caption or pass `adjust = "tukey"`.
#'
#' @param object An `ofe_fit` from [fit_ofe()].
#' @param term Character; the factor whose means are to be lettered.
#' @param by Character or `NULL`; letter the means within each level of a
#'   second factor, as in `ofe_lsd(fit, "treat", by = "zone")` for a
#'   pseudo-environment analysis. Comparisons and letters stay inside a level;
#'   the same letter in two zones means nothing.
#' @param alpha Significance level for the letters (default 0.05).
#' @param adjust Multiplicity adjustment for the pairwise p-values: `"none"`
#'   (default), `"tukey"`, `"sidak"`, or any method taken by
#'   [stats::p.adjust()]. Ignored when `use = "lsd"`.
#' @param use `"pairwise"` to judge each comparison on its own SED, or
#'   `"lsd"` to apply a single average LSD to all of them (see above).
#' @param decreasing Logical; assign `"a"` to the largest mean (default). Set
#'   `FALSE` where smaller is better, such as a disease score.
#' @param sort Logical; order the rows by mean rather than by factor level.
#' @param level Numeric; confidence level for the interval on each mean.
#'
#' @return A data frame with the `by` level where given, the `term` level,
#'   `estimate`, `std.error`, `lower`, `upper` and `group` (the letters).
#'   Attached as attributes: `lsd`, a one-row-per-group data frame of the
#'   average SED, the LSD, the degrees of freedom and the settings used; and
#'   `comparisons`, the full pairwise table behind the letters.
#'
#' @seealso [ofe_means()] for the means and comparisons on their own, and
#'   [compact_letters()] to letter comparisons from any other source.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 1)
#' fit <- fit_ofe(dense_response ~ treat, data = sim$grid)
#'
#' tab <- ofe_lsd(fit, "treat")
#' tab
#' attr(tab, "lsd")
#'
#' # Protected against multiplicity instead
#' ofe_lsd(fit, "treat", adjust = "tukey")
#'
#' @export
ofe_lsd <- function(object, term, by = NULL, alpha = 0.05, adjust = "none",
                    use = c("pairwise", "lsd"), decreasing = TRUE,
                    sort = TRUE, level = 0.95) {
  use <- match.arg(use)
  adjust <- match.arg(adjust, .ofe_adjust_methods)
  if (use == "lsd" && adjust != "none") {
    stop("`adjust` applies to the pairwise p-values and has no meaning for a ",
         "single average LSD. Use `use = \"pairwise\"`, or leave ",
         "`adjust = \"none\"`.", call. = FALSE)
  }
  cells <- .ofe_cells(object, term, by)
  den_df <- cells$den_df
  tcrit <- stats::qt(1 - alpha / 2, df = den_df)
  qlev <- stats::qt(1 - (1 - level) / 2, df = den_df)

  grp <- if (is.null(by)) rep(1L, length(cells$level)) else
    match(cells$by, cells$by_levels)

  means_out <- list()
  cmp_out <- list()
  lsd_out <- list()

  for (g in sort.int(unique(grp))) {
    i <- which(grp == g)
    est <- cells$estimate[i]
    V <- cells$V[i, i, drop = FALSE]
    lv <- cells$level[i]
    by_lab <- if (is.null(by)) NA_character_ else cells$by[i][1]

    cmp <- .ofe_pairs(est, V, lv, den_df, adjust)
    avg_sed <- mean(cmp$std.error)
    lsd_value <- tcrit * avg_sed

    if (use == "lsd") {
      # A single LSD, as a published table means it: the comparison is made
      # against the average SED rather than each pair's own.
      cmp$p.value <- 2 * stats::pt(abs(cmp$estimate) / avg_sed, den_df,
                                   lower.tail = FALSE)
      sig <- abs(cmp$estimate) > lsd_value
      cmp_for_letters <- cmp
      cmp_for_letters$p.value <- ifelse(sig, 0, 1)
    } else {
      cmp_for_letters <- cmp
    }

    mu <- stats::setNames(est, lv)
    letters_g <- compact_letters(cmp_for_letters, alpha = alpha, means = mu,
                                 decreasing = decreasing)

    se <- sqrt(diag(V))
    tab <- data.frame(level = lv, estimate = unname(est),
                      std.error = unname(se),
                      lower = unname(est - qlev * se),
                      upper = unname(est + qlev * se),
                      group = unname(letters_g[lv]),
                      stringsAsFactors = FALSE)
    if (sort) tab <- tab[match(names(letters_g), tab$level), , drop = FALSE]
    names(tab)[1] <- term

    lsd_row <- data.frame(average_sed = avg_sed, lsd = lsd_value, df = den_df,
                          alpha = alpha, adjust = adjust, use = use,
                          stringsAsFactors = FALSE)
    if (!is.null(by)) {
      addb <- function(x) cbind(stats::setNames(
        data.frame(rep(by_lab, nrow(x)), stringsAsFactors = FALSE), by), x)
      tab <- addb(tab)
      cmp <- addb(cmp)
      lsd_row <- addb(lsd_row)
    }
    means_out[[length(means_out) + 1L]] <- tab
    cmp_out[[length(cmp_out) + 1L]] <- cmp
    lsd_out[[length(lsd_out) + 1L]] <- lsd_row
  }

  out <- do.call(rbind, means_out)
  rownames(out) <- NULL
  attr(out, "lsd") <- do.call(rbind, lsd_out)
  cmp_all <- do.call(rbind, cmp_out)
  rownames(cmp_all) <- NULL
  attr(out, "comparisons") <- cmp_all
  out
}
