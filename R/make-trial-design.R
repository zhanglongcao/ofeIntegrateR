#' Lay out an on-farm strip trial
#'
#' Builds the trial map: treatment plots of a given width and length, arranged
#' in replicate blocks, snapped to a lattice at the resolution the dense layer
#' will be gridded to. The result is the plan you take to the field, and it has
#' the same columns as [grid_dense_layer()], so sampling locations can be chosen
#' from it with [place_point_samples()] before any yield data exist.
#'
#' Two layouts are supported, following the geometry compared in the GRDC
#' AAGI-CU-RD-OFE project:
#'
#' \describe{
#'   \item{`"strip"`}{Every plot runs the full length of the trial, side by
#'     side across its width. The simplest layout to drive and the usual
#'     default.}
#'   \item{`"stack"`}{The same plots split into two tiers separated by a buffer,
#'     halving the width and roughly doubling the length. Useful where the
#'     paddock is too narrow for a single row of strips, and it gives a squarer
#'     footprint, which changes how much spatial variation the trial spans.}
#' }
#'
#' `randomise` is not a formality. Grower strip trials have conventionally been
#' randomised by analogy with small-plot work, but under strong spatial
#' correlation a systematic arrangement can estimate treatment contrasts more
#' precisely, because it spreads each treatment evenly across the field's
#' spatial gradient. Both are produced here so the choice can be made
#' deliberately rather than by habit.
#'
#' @param treatments Character vector of treatment labels, e.g.
#'   `c("N0", "N60", "N120")`.
#' @param n_rep Number of replicate blocks. Each block contains every treatment
#'   once.
#' @param layout `"strip"` or `"stack"`, as described above.
#' @param plot_width Width of one plot in metres, normally the working width of
#'   the applicator or seeder.
#' @param plot_length Length of one plot in metres, along the direction of
#'   travel. For `"stack"` this is the length of each tier, not of the trial.
#' @param cell_size Lattice resolution in metres, normally the harvester swath,
#'   so the design and the harvested data share a grid.
#' @param gap Buffer between the two tiers of a `"stack"` layout, in metres.
#'   Ignored for `"strip"`.
#' @param randomise Randomise treatment order within each replicate block
#'   (a randomised complete block design). `FALSE` repeats the same order in
#'   every block, giving a systematic arrangement.
#' @param origin Length-2 numeric giving the south-west corner of the trial in
#'   field coordinates, for placing the trial within a paddock.
#' @param seed Optional integer seed, so a randomised layout is reproducible.
#'
#' @return A data frame with one row per lattice cell: `row`, `col`, `x`, `y`
#'   (cell centres in metres), `treat`, `rep`, `plot` (plot identifier), and
#'   `tier` for stacked layouts. Cells falling in a buffer carry `NA`
#'   treatment. The design is recorded in `attr(, "design")`.
#'
#' @examples
#' # Three nitrogen rates, four replicates, 18 m applicator, 200 m runs
#' d <- make_trial_design(c("N0", "N60", "N120"), n_rep = 4,
#'                        plot_width = 18, plot_length = 200, seed = 1)
#' attr(d, "design")$trial_width
#' table(d$treat)
#'
#' # Treatment order differs between blocks when randomised
#' unique(d[, c("rep", "plot", "treat")])[1:6, ]
#'
#' # The same plots stacked into two tiers, for a narrower paddock
#' s <- make_trial_design(c("N0", "N60", "N120"), n_rep = 4, layout = "stack",
#'                        plot_width = 18, plot_length = 200, gap = 10, seed = 1)
#' attr(s, "design")$trial_width
#'
#' # Feeds the sampling design directly
#' set.seed(1)
#' cores <- place_point_samples(d[!is.na(d$treat), ], n = 24,
#'                              design = "stratified",
#'                              coords = c("x", "y"), strata = "treat")
#' table(cores$treat)
#'
#' @export
make_trial_design <- function(treatments,
                               n_rep = 4,
                               layout = c("strip", "stack"),
                               plot_width = 18,
                               plot_length = 200,
                               cell_size = 9,
                               gap = 0,
                               randomise = TRUE,
                               origin = c(0, 0),
                               seed = NULL) {
  layout <- match.arg(layout)
  treatments <- as.character(treatments)
  n_treat <- length(treatments)

  if (n_treat < 2L) {
    stop("`treatments` must name at least two treatments.", call. = FALSE)
  }
  if (n_rep < 1L) stop("`n_rep` must be at least 1.", call. = FALSE)
  if (plot_width <= 0 || plot_length <= 0 || cell_size <= 0) {
    stop("`plot_width`, `plot_length` and `cell_size` must be positive.",
         call. = FALSE)
  }
  if (cell_size > plot_width) {
    stop("`cell_size` (", cell_size, " m) is wider than a plot (", plot_width,
         " m), so no cell could lie inside a single treatment. Use a finer ",
         "lattice or wider plots.", call. = FALSE)
  }
  n_plots <- n_treat * n_rep
  if (layout == "stack" && n_plots %% 2L != 0L) {
    stop("A stacked layout splits the plots into two tiers, so ",
         "`length(treatments) * n_rep` must be even; here it is ", n_plots, ".",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  # Randomised complete blocks: every block holds each treatment once, in a
  # fresh order. Systematic keeps one order throughout.
  block_order <- function() if (randomise) sample(treatments) else treatments
  book <- do.call(rbind, lapply(seq_len(n_rep), function(b) {
    data.frame(rep = b, position = seq_len(n_treat), treat = block_order(),
               stringsAsFactors = FALSE)
  }))
  book$plot <- seq_len(nrow(book))

  if (layout == "strip") {
    book$tier <- 1L
    book$col_index <- book$plot
    tier_len <- plot_length
    trial_width <- n_plots * plot_width
    trial_length <- plot_length
  } else {
    per_tier <- n_plots / 2L
    book$tier <- ifelse(book$plot <= per_tier, 1L, 2L)
    book$col_index <- ifelse(book$tier == 1L, book$plot, book$plot - per_tier)
    tier_len <- plot_length
    trial_width <- per_tier * plot_width
    trial_length <- 2 * plot_length + gap
  }

  n_col <- floor(trial_width / cell_size)
  n_row <- floor(trial_length / cell_size)
  if (n_col < 1L || n_row < 1L) {
    stop("`cell_size` leaves no whole cells in the trial.", call. = FALSE)
  }

  g <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  g$x <- origin[1] + (g$col - 0.5) * cell_size
  g$y <- origin[2] + (g$row - 0.5) * cell_size

  # Which plot does each cell centre fall in?
  rel_x <- g$x - origin[1]
  rel_y <- g$y - origin[2]
  col_index <- floor(rel_x / plot_width) + 1L

  if (layout == "strip") {
    tier <- rep(1L, nrow(g))
    in_buffer <- rep(FALSE, nrow(g))
  } else {
    tier <- ifelse(rel_y <= tier_len, 1L,
                   ifelse(rel_y > tier_len + gap, 2L, NA_integer_))
    in_buffer <- is.na(tier)
  }

  key <- paste(tier, col_index)
  book_key <- paste(book$tier, book$col_index)
  idx <- match(key, book_key)

  g$treat <- factor(book$treat[idx], levels = treatments)
  g$rep <- book$rep[idx]
  g$plot <- book$plot[idx]
  g$tier <- tier
  g$treat[in_buffer | col_index > max(book$col_index)] <- NA
  g$rep[in_buffer] <- NA_integer_
  g$plot[in_buffer] <- NA_integer_

  g <- g[order(g$row, g$col), c("row", "col", "x", "y", "treat", "rep",
                                 "plot", "tier")]
  rownames(g) <- NULL
  attr(g, "design") <- list(
    layout = layout, treatments = treatments, n_rep = n_rep,
    n_plots = n_plots, plot_width = plot_width, plot_length = plot_length,
    cell_size = cell_size, gap = if (layout == "stack") gap else 0,
    randomise = randomise, trial_width = trial_width,
    trial_length = trial_length, area_ha = trial_width * trial_length / 10000)
  g
}
