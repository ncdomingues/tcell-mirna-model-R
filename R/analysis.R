# analysis.R -------------------------------------------------------------------
#
# Summaries and plots on top of a `ginsim_run()` result. Base R graphics only --
# the point of this repo is that it runs on a bare R install with no packages.
#
# The 2017 plotting scripts were ~200 lines of copy-pasted `lines(...)` calls
# with a hand-maintained colour vector and legend, which is how several of the
# original figures ended up with two series sharing a colour and one label
# pointing at the wrong node. These helpers take node names and derive the rest.
# ------------------------------------------------------------------------------


# The master transcription factors that define each CD4+ T-helper fate. This is
# the readout the whole model exists to produce.
MASTER_TFS <- c(Th1 = "TBX21", Th2 = "GATA3", Th17 = "RORC", iTreg = "FOXP3",
                Tfh = "BCL6", Th22 = "AHR", Th9 = "SPI1-PU1")

CYTOKINES <- c("IL2", "IL4", "IL6", "IL9", "IL10", "IL17", "IL21", "IL22",
               "IFNG", "TNF", "TGFB")

MIRNAS <- c("miR-34c-5p", "miR-155-5p")


#' A colourblind-safe palette, recycled with distinct line types
#'
#' Okabe-Ito, which stays legible in print and to the ~8% of men with red-green
#' colour vision deficiency -- unlike the 2017 figures' red/green pairs.
series_style <- function(n) {
  okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#000000")
  list(col = okabe[(seq_len(n) - 1L) %% length(okabe) + 1L],
       lty = ((seq_len(n) - 1L) %/% length(okabe)) + 1L)
}


#' Plot mean trajectories for a set of nodes
#'
#' @param run    a `ginsim_run()` result
#' @param names  node names to draw
#' @param title  plot title
plot_trajectories <- function(run, names, title) {
  style <- series_style(length(names))
  time  <- seq_len(ncol(run$evol)) - 1L
  ymax  <- max(1, ceiling(max(run$evol[node_id(run$nodes, names), ])))

  op <- par(mar = c(4.5, 4.5, 3, 9), xpd = FALSE)
  on.exit(par(op), add = TRUE)

  plot(NA, xlim = range(time), ylim = c(0, ymax),
       xlab = "Update step", ylab = "Mean level across the ensemble",
       main = title, bty = "l", las = 1)
  grid(col = "grey92", lty = 1)

  for (i in seq_along(names)) {
    lines(time, run$evol[node_id(run$nodes, names[i]), ],
          col = style$col[i], lty = style$lty[i], lwd = 2)
  }

  par(xpd = NA)
  legend(par("usr")[2] + 0.02 * diff(par("usr")[1:2]), par("usr")[4],
         legend = names, col = style$col, lty = style$lty, lwd = 2,
         bty = "n", cex = 0.75, seg.len = 1.6)
}


#' Grouped bar chart of final activity across conditions
#'
#' @param activity  named-list of `final_activity()` vectors, one per condition
#' @param names     node names to show
plot_condition_bars <- function(activity, names, title,
                                ylab = "Mean occupancy over the second half of the run") {
  mat <- vapply(activity, function(a) a[names], numeric(length(names)))
  mat <- matrix(mat, nrow = length(names),
                dimnames = list(names, names(activity)))
  style <- series_style(ncol(mat))

  op <- par(mar = c(6.5, 4.5, 3, 9), xpd = FALSE)
  on.exit(par(op), add = TRUE)

  bp <- barplot(t(mat), beside = TRUE, col = style$col, border = NA,
                ylim = c(0, max(1, max(mat, na.rm = TRUE) * 1.1)),
                ylab = ylab, main = title, las = 2, cex.names = 0.8)
  abline(h = 0, col = "grey40")

  par(xpd = NA)
  legend(max(bp) + 0.5, par("usr")[4], legend = colnames(mat),
         fill = style$col, border = NA, bty = "n", cex = 0.8)
  invisible(mat)
}


#' Open a PNG device with consistent sizing, run `expr`, close it
to_png <- function(path, expr, width = 1400, height = 850, res = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  message("wrote ", path)
}


#' Long-format data frame of every node's trajectory, for writing to CSV
trajectory_table <- function(run, condition) {
  n_time <- ncol(run$evol)
  data.frame(
    condition = condition,
    node      = rep(run$nodes[, 2], each = n_time),
    step      = rep(seq_len(n_time) - 1L, times = nrow(run$evol)),
    level     = as.vector(t(run$evol)),
    stringsAsFactors = FALSE
  )
}


#' Distinct attractors reached, with how many trajectories reached each
attractor_summary <- function(run) {
  if (run$n_steady == 0) {
    return(data.frame(attractor = integer(0), trajectories = integer(0),
                      stringsAsFactors = FALSE))
  }
  key   <- apply(run$ss, 1, paste, collapse = "")
  tab   <- sort(table(key), decreasing = TRUE)
  data.frame(attractor = seq_along(tab), trajectories = as.integer(tab),
             signature = names(tab), stringsAsFactors = FALSE)
}
