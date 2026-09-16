# Parsing and evaluation of ASReml-style residual structure formulae.
#
# The point of this file is that `~ ar1(row):ar1(col)` and
# `~ dsum(~ ar1(row):ar1(col) | zone)` should mean the same thing to
# [fit_ofe()] as they do to asreml, so that a workflow written against one
# runs against the other. Nothing here is exported: the user-facing surface is
# the `residual` argument of [fit_ofe()] and the formula built by
# [adaptive_residual()].

#' @keywords internal
#' @noRd
.ofe_struct_fns <- c("id", "ar1", "exp", "diag")

#' Split a structure expression on `:` into its separable components
#' @keywords internal
#' @noRd
.parse_struct <- function(expr) {
  comps <- list()
  walk <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), ":")) {
      walk(e[[2]])
      walk(e[[3]])
    } else {
      comps[[length(comps) + 1L]] <<- .parse_comp(e)
    }
  }
  walk(expr)
  comps
}

#' @keywords internal
#' @noRd
.parse_comp <- function(e) {
  # A bare variable name means id(), matching asreml's shorthand.
  if (is.name(e)) return(list(fn = "id", var = as.character(e)))
  if (!is.call(e)) {
    stop("Cannot interpret `", deparse(e), "` as a residual structure.",
         call. = FALSE)
  }
  fn <- as.character(e[[1]])
  if (!fn %in% .ofe_struct_fns) {
    stop("Unsupported residual structure `", fn, "()`. Supported: ",
         paste0(.ofe_struct_fns, "()", collapse = ", "), ".", call. = FALSE)
  }
  if (length(e) < 2L || !is.name(e[[2]])) {
    stop("`", fn, "()` needs a single variable name, e.g. ", fn, "(row).",
         call. = FALSE)
  }
  list(fn = fn, var = as.character(e[[2]]))
}

#' Split a residual formula into independent sections
#'
#' Returns a list of sections, each `list(struct = <list of components>,
#' level = <section level or NA>, var = <section variable or NA>)`. A formula
#' with no `dsum()` gives one section covering every observation; `dsum()`
#' gives one section per named level, each with its own parameters, exactly as
#' asreml does.
#' @keywords internal
#' @noRd
.parse_residual <- function(residual, data) {
  if (!inherits(residual, "formula") || length(residual) != 2L) {
    stop("`residual` must be a one-sided formula, e.g. ~ ar1(row):ar1(col).",
         call. = FALSE)
  }
  terms_list <- list()
  walk <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), "+")) {
      walk(e[[2]])
      walk(e[[3]])
    } else {
      terms_list[[length(terms_list) + 1L]] <<- e
    }
  }
  walk(residual[[2]])

  sections <- list()
  for (e in terms_list) {
    is_dsum <- is.call(e) && identical(as.character(e[[1]]), "dsum")
    if (!is_dsum) {
      if (length(terms_list) > 1L) {
        stop("A residual formula with `+` must be a sum of `dsum()` terms.",
             call. = FALSE)
      }
      sections[[1L]] <- list(struct = .parse_struct(e), var = NA_character_,
                             level = NA_character_)
      next
    }
    args <- as.list(e)[-1L]
    nms <- names(args)
    if (is.null(nms)) nms <- rep("", length(args))
    f_idx <- which(nms == "" | nms == "formula")[1L]
    if (is.na(f_idx)) {
      stop("`dsum()` needs a formula, e.g. dsum(~ ar1(row):ar1(col) | zone).",
           call. = FALSE)
    }
    inner <- args[[f_idx]]
    if (!(is.call(inner) && identical(as.character(inner[[1]]), "~"))) {
      stop("The first argument of `dsum()` must be a formula.", call. = FALSE)
    }
    body <- inner[[length(inner)]]
    if (!(is.call(body) && identical(as.character(body[[1]]), "|"))) {
      stop("`dsum()` needs a section variable after `|`, e.g. ",
           "dsum(~ ar1(row):ar1(col) | zone).", call. = FALSE)
    }
    struct <- .parse_struct(body[[2]])
    svar <- as.character(body[[3]])
    if (is.null(data[[svar]])) {
      stop("Section variable `", svar, "` not found in `data`.", call. = FALSE)
    }
    all_levels <- levels(factor(data[[svar]]))
    lv_idx <- which(nms == "levels")[1L]
    lv <- if (is.na(lv_idx)) all_levels else {
      as.character(eval(args[[lv_idx]], envir = data, enclos = parent.frame()))
    }
    unknown <- setdiff(lv, all_levels)
    if (length(unknown) > 0L) {
      stop("`dsum()` names level(s) not present in `", svar, "`: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    for (l in lv) {
      sections[[length(sections) + 1L]] <-
        list(struct = struct, var = svar, level = l)
    }
  }

  if (length(sections) == 0L) {
    stop("`residual` defines no structure.", call. = FALSE)
  }
  covered <- sections[[1L]]$var
  if (!is.na(covered)) {
    lv_used <- vapply(sections, function(s) s$level, character(1))
    if (anyDuplicated(lv_used)) {
      stop("A level of `", covered, "` appears in more than one `dsum()`: ",
           paste(unique(lv_used[duplicated(lv_used)]), collapse = ", "), ".",
           call. = FALSE)
    }
  }
  sections
}

#' Ordinal positions for an AR1 dimension
#'
#' AR1 is defined on equally spaced, ordered levels. A factor uses its level
#' order, which is asreml's convention and lets a caller control the lattice by
#' setting levels. A numeric column is read as a position on a lattice whose
#' spacing is the commonest gap between its distinct values, so that a row
#' absent from the data -- dropped for a missing yield, or never harvested --
#' still contributes its lag instead of being closed up. Ranking the distinct
#' values instead would make rows 5 and 7 adjacent, biasing every correlation
#' upwards on exactly the incomplete lattices that OFE data arrive in.
#' @keywords internal
#' @noRd
.ordinal_position <- function(v) {
  if (is.factor(v)) return(as.integer(v))
  if (is.character(v)) return(as.integer(factor(v)))
  u <- sort(unique(v))
  if (length(u) < 2L) return(rep(1, length(v)))
  gaps <- diff(u)
  tab <- table(signif(gaps, 8))
  spacing <- as.numeric(names(tab)[which.max(tab)])
  if (!is.finite(spacing) || spacing <= 0) return(match(v, u))
  (v - u[1]) / spacing + 1
}

#' Precompute the fixed pieces of every section
#'
#' Lag and equality matrices do not depend on the variance parameters, so they
#' are built once here and reused at every likelihood evaluation. This is what
#' keeps a few hundred REML iterations affordable.
#' @keywords internal
#' @noRd
.build_residual_spec <- function(residual, data) {
  sections <- .parse_residual(residual, data)
  n <- nrow(data)
  p_name <- character(0)
  p_type <- character(0)

  add_par <- function(name, type) {
    p_name[length(p_name) + 1L] <<- name
    p_type[length(p_type) + 1L] <<- type
    length(p_name)
  }

  multi <- length(sections) > 1L
  out <- vector("list", length(sections))
  used <- logical(n)

  for (s in seq_along(sections)) {
    sec <- sections[[s]]
    idx <- if (is.na(sec$var)) seq_len(n) else which(as.character(data[[sec$var]]) == sec$level)
    if (length(idx) == 0L) {
      stop("Section `", sec$level, "` of `", sec$var, "` has no observations.",
           call. = FALSE)
    }
    used[idx] <- TRUE
    tag <- if (is.na(sec$level)) "R" else paste0(sec$var, "!", sec$level)

    comps <- vector("list", length(sec$struct))
    for (k in seq_along(sec$struct)) {
      cp <- sec$struct[[k]]
      v <- data[[cp$var]]
      if (is.null(v)) {
        stop("Residual structure variable `", cp$var, "` not found in `data`.",
             call. = FALSE)
      }
      fn <- cp$fn
      if (fn == "ar1") {
        pos <- .ordinal_position(v)[idx]
        if (length(unique(pos)) < 2L) {
          warning("`ar1(", cp$var, ")` has one level in section ", tag,
                  "; fitting `id(", cp$var, ")` there instead. ",
                  "adaptive_residual() makes this choice for you.",
                  call. = FALSE)
          fn <- "id"
        } else {
          cp$D <- abs(outer(pos, pos, "-"))
          cp$par <- add_par(paste0(tag, "!", cp$var, "!cor"), "cor")
        }
      }
      if (fn == "exp") {
        if (!is.numeric(v)) {
          stop("`exp(", cp$var, ")` needs a numeric coordinate.", call. = FALSE)
        }
        d <- abs(outer(v[idx], v[idx], "-"))
        if (max(d) <= 0) {
          stop("`exp(", cp$var, ")` has no spread in section ", tag, ".",
               call. = FALSE)
        }
        cp$D <- d
        cp$scale <- stats::median(d[d > 0])
        cp$par <- add_par(paste0(tag, "!", cp$var, "!range"), "range")
      }
      if (fn %in% c("id", "diag")) {
        f <- factor(v[idx])
        l <- as.integer(f)
        cp$Ieq <- outer(l, l, "==") + 0
        cp$lev <- l
        cp$labels <- levels(f)
        if (fn == "diag") {
          if (nlevels(f) < 2L) {
            stop("`diag(", cp$var, ")` needs at least two levels in section ",
                 tag, ".", call. = FALSE)
          }
          # First level is the reference at 1; the rest are relative variances.
          cp$par <- vapply(levels(f)[-1L], function(lb)
            add_par(paste0(tag, "!", cp$var, "!", lb, "!var"), "var"),
            integer(1), USE.NAMES = FALSE)
        }
      }
      cp$fn <- fn
      comps[[k]] <- cp
    }

    nu_par <- NA_integer_
    if (multi && s > 1L) nu_par <- add_par(paste0(tag, "!var"), "var")

    out[[s]] <- list(idx = idx, comps = comps, nu_par = nu_par, tag = tag,
                     n = length(idx))
  }

  if (!all(used)) {
    stop(sum(!used), " observation(s) are in no residual section. Every row ",
         "must belong to a level named by `dsum()`.", call. = FALSE)
  }

  ptab <- data.frame(name = p_name, type = p_type, stringsAsFactors = FALSE)
  list(sections = out, ptab = ptab, n_par = nrow(ptab))
}

#' Evaluate one section's covariance at the current parameters
#' @keywords internal
#' @noRd
.eval_section <- function(sec, theta) {
  m <- sec$n
  M <- matrix(1, m, m)
  for (cp in sec$comps) {
    Mi <- switch(
      cp$fn,
      id   = cp$Ieq,
      ar1  = theta[cp$par]^cp$D,
      exp  = exp(-cp$D / theta[cp$par]),
      diag = {
        v <- c(1, theta[cp$par])
        cp$Ieq * matrix(v[cp$lev], m, m)
      }
    )
    M <- M * Mi
  }
  if (!is.na(sec$nu_par)) M <- M * theta[sec$nu_par]
  M
}

#' Assemble the residual covariance (relative to sigma^2) for the whole trial
#' @keywords internal
#' @noRd
.residual_matrix <- function(spec, theta, n) {
  H <- matrix(0, n, n)
  for (sec in spec$sections) {
    H[sec$idx, sec$idx] <- .eval_section(sec, theta)
  }
  H
}
