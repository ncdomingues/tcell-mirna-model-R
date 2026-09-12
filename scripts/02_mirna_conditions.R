# 02_mirna_conditions.R --------------------------------------------------------
#
# The question the model was built to answer: what do miR-34c-5p and miR-155-5p
# do to CD4+ T-helper fate, separately and together?
#
# Four versions of the same network under identical TCR + CD28 stimulation --
# no miRNA, miR-34c-5p only, miR-155-5p only, both. A miRNA is "absent" when its
# own transcriptional rules are dropped: it then stays at 0 for the whole run, so
# every `!miR-...` term elsewhere in the table evaluates as if it were not there.
#
#   Rscript scripts/02_mirna_conditions.R
# ------------------------------------------------------------------------------

source("R/ginsim_engine.R")
source("R/model_io.R")
source("R/analysis.R")

NTIME <- 200
NREP  <- 200
SEED  <- 2019

model <- load_model(rules_file = "rules_corrected.txt")
state <- stimulation_state(model$nodes, "IL2")

conditions <- list(
  "no miRNA"     = MIRNAS,
  "miR-34c-5p"   = "miR-155-5p",
  "miR-155-5p"   = "miR-34c-5p",
  "both"         = character(0)
)

runs <- list()
for (cond in names(conditions)) {
  message("running ", cond, " ...")
  rules <- knock_out(model$rules, conditions[[cond]])
  runs[[cond]] <- ginsim_run(state, model$nodes, rules,
                             ntime = NTIME, nrep = NREP, seed = SEED)
}

activity <- lapply(runs, mean_activity)

dir.create("output", showWarnings = FALSE)

# --- sanity check: a knocked-out miRNA must actually stay off -----------------
for (cond in names(conditions)) {
  off <- conditions[[cond]]
  if (length(off) && any(activity[[cond]][off] > 0)) {
    stop("knock-out failed: ", paste(off, collapse = ", "), " active in '", cond, "'")
  }
}

# --- master transcription factors --------------------------------------------
tf <- data.frame(subset = names(MASTER_TFS), node = unname(MASTER_TFS),
                 stringsAsFactors = FALSE)
for (cond in names(conditions)) {
  tf[[cond]] <- round(activity[[cond]][tf$node], 4)
}
tf$`effect of both` <- round(tf$both - tf$`no miRNA`, 4)
write.csv(tf, "output/mirna_master_tfs.csv", row.names = FALSE)
print(tf, row.names = FALSE)

# --- every node, for anyone who wants to dig ---------------------------------
all_nodes <- data.frame(node = model$nodes[, 2], stringsAsFactors = FALSE)
for (cond in names(conditions)) all_nodes[[cond]] <- round(activity[[cond]], 4)
write.csv(all_nodes, "output/mirna_all_nodes.csv", row.names = FALSE)

# --- which nodes move most when both miRNAs are present? ---------------------
delta  <- activity$both - activity$`no miRNA`
moved  <- delta[order(abs(delta), decreasing = TRUE)]
moved  <- moved[abs(moved) > 0.01]
write.csv(data.frame(node = names(moved), change = round(unname(moved), 4)),
          "output/mirna_largest_changes.csv", row.names = FALSE)
cat("\nLargest shifts, both miRNAs vs none:\n")
print(round(head(moved, 15), 3))

# --- figures ------------------------------------------------------------------
to_png("figures/mirna_master_tfs.png",
       plot_condition_bars(activity, unname(MASTER_TFS),
                           "T-helper master regulators by miRNA condition"))

to_png("figures/mirna_cytokines.png",
       plot_condition_bars(activity, CYTOKINES,
                           "Cytokine output by miRNA condition"))

to_png("figures/mirna_trajectories_both.png",
       plot_trajectories(runs$both, c(MIRNAS, unname(MASTER_TFS)),
                         "Both miRNAs present: master regulators over time"))

message("done")
