# 05_compare_tables.R ----------------------------------------------------------
#
# Two rule tables for the same biology, run through the same engine under the
# same protocol, so that any difference is attributable to the tables.
#
#   A. model/rules_corrected.txt -- the 2017 working table, 91 nodes, 86 rules.
#      What was actually simulated for the thesis.
#   B. model/thesis_table5/     -- Table 5 as printed in the thesis, 81 nodes,
#      89 rules. Ported from the digitization in the Python rebuild
#      (github.com/ncdomingues/tcell-mirna-model), which transcribed it from
#      Results_PhDthesis.docx. See model/thesis_table5/README.md.
#
# They are not the same network. The working table is the February 2017 file;
# Table 5 is what went into the write-up two years later. This script is the
# first time the two have been run against each other.
#
# The stimulus is derived the same way for both: every node no rule targets is
# switched on, except exogenous polarising cytokines. That is a naive CD4+
# T cell with antigen, CD3/CD28 engagement and its constitutive co-factors, and
# nothing else -- the TCR + IL2 condition, where IL2 is autocrine.
#
#   Rscript scripts/05_compare_tables.R
# ------------------------------------------------------------------------------

source("R/ginsim_engine.R")
source("R/model_io.R")
source("R/analysis.R")
source("R/table_diff.R")

NTIME <- 200
NREP  <- 200
SEED  <- 2019

TABLES <- list(
  "2017 working table" = list(
    model = load_model(rules_file = "rules_corrected.txt"),
    # IL12, IL23 and IL1 are exogenous inputs in this table; held off.
    exogenous = c("IL12", "IL23", "IL1")
  ),
  "thesis Table 5" = list(
    model = load_model("model/thesis_table5"),
    exogenous = character(0)   # no exogenous cytokine inputs in this table
  )
)

CONDITIONS <- list(
  "no miRNA"   = MIRNAS,
  "miR-34c-5p" = "miR-155-5p",
  "miR-155-5p" = "miR-34c-5p",
  "both"       = character(0)
)

activity <- list()
for (tbl in names(TABLES)) {
  nodes <- TABLES[[tbl]]$model$nodes
  rules <- TABLES[[tbl]]$model$rules
  stim  <- setdiff(unregulated_nodes(nodes, rules), TABLES[[tbl]]$exogenous)
  message(tbl, ": ", nrow(nodes), " nodes, ", nrow(rules), " rules")
  message("  stimulus: ", paste(stim, collapse = ", "))

  for (cond in names(CONDITIONS)) {
    message("  running ", cond, " ...")
    r <- suppressWarnings(ginsim_run(make_state(nodes, stim), nodes,
                                     knock_out(rules, CONDITIONS[[cond]]),
                                     ntime = NTIME, nrep = NREP, seed = SEED))
    activity[[tbl]][[cond]] <- mean_activity(r)
  }
}

dir.create("output", showWarnings = FALSE)

# --- where the two tables actually differ -------------------------------------
d <- diff_tables(TABLES[["2017 working table"]]$model$nodes,
                 TABLES[["2017 working table"]]$model$rules,
                 TABLES[["thesis Table 5"]]$model$nodes,
                 TABLES[["thesis Table 5"]]$model$rules,
                 label_a = "working_table", label_b = "table5")
write.csv(d, "output/table_diff.csv", row.names = FALSE)
cat("\nRule-level diff (", nrow(d), " node/value pairs):\n", sep = "")
print(table(d$verdict))

# --- side by side, master transcription factors -------------------------------
rows <- list()
for (tbl in names(TABLES)) {
  for (i in seq_along(MASTER_TFS)) {
    node <- unname(MASTER_TFS[i])
    if (!(node %in% names(activity[[tbl]][["no miRNA"]]))) next
    vals <- vapply(names(CONDITIONS),
                   function(c) unname(activity[[tbl]][[c]][node]), 0)
    rows[[length(rows) + 1]] <- data.frame(
      table = tbl, subset = names(MASTER_TFS)[i], node = node,
      `no miRNA` = round(vals[["no miRNA"]], 4),
      `miR-34c-5p` = round(vals[["miR-34c-5p"]], 4),
      `miR-155-5p` = round(vals[["miR-155-5p"]], 4),
      both = round(vals[["both"]], 4),
      effect_of_both = round(vals[["both"]] - vals[["no miRNA"]], 4),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
side <- do.call(rbind, rows)
write.csv(side, "output/table_comparison_master_tfs.csv", row.names = FALSE)
print(side, row.names = FALSE)

# --- do the two tables agree on the direction of the combined effect? ---------
direction <- function(x, tol = 0.05) {
  if (is.na(x)) NA_character_ else if (x > tol) "up" else if (x < -tol) "down" else "flat"
}
wide <- reshape(side[, c("table", "subset", "effect_of_both")],
                idvar = "subset", timevar = "table", direction = "wide")
names(wide) <- c("subset", "working_table", "table5")
wide$working_dir <- vapply(wide$working_table, direction, "")
wide$table5_dir  <- vapply(wide$table5, direction, "")
wide$agree <- ifelse(wide$working_dir == wide$table5_dir, "yes", "no")
write.csv(wide, "output/table_comparison_agreement.csv", row.names = FALSE)

cat("\nDirection of the combined-miRNA effect, both tables:\n")
print(wide[, c("subset", "working_table", "working_dir", "table5", "table5_dir", "agree")],
      row.names = FALSE)
cat("\nagreement:", sum(wide$agree == "yes"), "of", nrow(wide), "subsets\n")

# --- which miRNA carries the combined effect? ---------------------------------
driver <- data.frame(subset = names(MASTER_TFS), stringsAsFactors = FALSE)
for (tbl in names(TABLES)) {
  driver[[tbl]] <- vapply(unname(MASTER_TFS), function(node) {
    a <- activity[[tbl]]
    if (!(node %in% names(a[["no miRNA"]]))) return(NA_character_)
    base <- a[["no miRNA"]][node]
    d34  <- abs(a[["miR-34c-5p"]][node] - base)
    d155 <- abs(a[["miR-155-5p"]][node] - base)
    if (max(d34, d155) < 0.05) "neither" else
      if (abs(d34 - d155) < 0.05) "both equally" else
        if (d34 > d155) "miR-34c-5p" else "miR-155-5p"
  }, "")
}
write.csv(driver, "output/table_comparison_driver.csv", row.names = FALSE)
cat("\nWhich miRNA moves each subset furthest from baseline:\n")
print(driver, row.names = FALSE)

# --- epistasis: in combination, does one miRNA override the other? ------------
#
# Different question from the one above. A miRNA can have the larger solo effect
# and still lose when both are present. Compare the `both` value against each
# solo value: if it lands on one of them, that miRNA is epistatic.
epi <- data.frame(subset = names(MASTER_TFS), stringsAsFactors = FALSE)
for (tbl in names(TABLES)) {
  epi[[tbl]] <- vapply(unname(MASTER_TFS), function(node) {
    a <- activity[[tbl]]
    if (!(node %in% names(a[["no miRNA"]]))) return(NA_character_)
    solo <- c("miR-34c-5p" = unname(a[["miR-34c-5p"]][node]),
              "miR-155-5p" = unname(a[["miR-155-5p"]][node]))
    if (abs(diff(solo)) < 0.05) return("solo effects alike")
    d <- abs(solo - unname(a[["both"]][node]))
    if (abs(d[1] - d[2]) < 0.05) "additive" else names(which.min(d))
  }, "")
}
write.csv(epi, "output/table_comparison_epistasis.csv", row.names = FALSE)
cat("\nWith both present, which solo outcome does the pair land on:\n")
print(epi, row.names = FALSE)

# --- how much of the time is each miRNA actually on? --------------------------
#
# The obvious explanation for a miRNA mattering less would be fewer target
# edges. Here it is the opposite -- Table 5 carries more than twice as many
# miR-34c-5p repression terms -- so the thing to measure is whether the miRNA
# is ever expressed in the first place.
on <- data.frame(node = c(MIRNAS, "STAT3", "RORC"), stringsAsFactors = FALSE)
for (tbl in names(TABLES)) on[[tbl]] <- round(activity[[tbl]][["both"]][on$node], 3)
cat("\nOccupancy with both miRNAs present:\n")
print(on, row.names = FALSE)

cat("\nThe miR-34c-5p promoter rule in each table:\n")
for (tbl in names(TABLES)) {
  r <- TABLES[[tbl]]$model$rules
  r <- r[trimws(r[[2]]) == "miR-34c-5p", ]
  cat("  ", tbl, "\n")
  for (i in seq_len(nrow(r))) cat("    value ", r[i, 3], ": ", r[i, 4], "\n", sep = "")
}
write.csv(on, "output/table_comparison_mirna_activity.csv", row.names = FALSE)

# --- cross-check against the Python rebuild's published numbers ---------------
#
# Same table (Table 5), different engine, different protocol, different units
# (activation frequency at the final state, as a percentage). The magnitudes are
# not comparable and are not compared; the direction of the combined effect is.
py_file <- "comparison/python_rebuild_phenotypes.csv"
if (file.exists(py_file)) {
  py   <- read.csv(py_file, row.names = 1, check.names = FALSE)
  mine <- side[side$table == "thesis Table 5", ]
  sub  <- intersect(rownames(py), mine$subset)
  cross <- data.frame(
    subset      = sub,
    python      = round(py[sub, "both"] - py[sub, "no miRNA"], 1),
    this_repo   = round(mine$both[match(sub, mine$subset)] -
                          mine$`no miRNA`[match(sub, mine$subset)], 3),
    stringsAsFactors = FALSE
  )
  cross$sign_agrees <- sign(cross$python) == sign(cross$this_repo)
  write.csv(cross, "output/table5_vs_python_rebuild.csv", row.names = FALSE)
  cat("\nTable 5 under two independent implementations",
      "(direction of the combined effect):\n")
  print(cross, row.names = FALSE)
  cat("signs agreeing:", sum(cross$sign_agrees), "of", nrow(cross), "\n")
}

# --- figure -------------------------------------------------------------------
to_png("figures/table_comparison.png",
       {
         op <- par(mfrow = c(1, 2), mar = c(6.5, 4.5, 3.5, 1), oma = c(0, 0, 0, 8))
         on.exit(par(op), add = TRUE)
         style <- series_style(length(CONDITIONS))
         for (tbl in names(TABLES)) {
           nodes_present <- unname(MASTER_TFS)[
             unname(MASTER_TFS) %in% names(activity[[tbl]][["no miRNA"]])]
           mat <- vapply(names(CONDITIONS),
                         function(c) activity[[tbl]][[c]][nodes_present],
                         numeric(length(nodes_present)))
           mat <- matrix(mat, nrow = length(nodes_present),
                         dimnames = list(nodes_present, names(CONDITIONS)))
           barplot(t(mat), beside = TRUE, col = style$col, border = NA,
                   ylim = c(0, 1), las = 2, cex.names = 0.8,
                   ylab = if (tbl == names(TABLES)[1])
                     "Mean occupancy over the second half of the run" else "",
                   main = tbl)
           abline(h = 0, col = "grey40")
         }
         par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
         plot(0, 0, type = "n", bty = "n", axes = FALSE)
         legend("right", legend = names(CONDITIONS), fill = style$col,
                border = NA, bty = "n", cex = 0.85)
       },
       width = 1700, height = 800)

message("done")
