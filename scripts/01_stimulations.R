# 01_stimulations.R ------------------------------------------------------------
#
# The 2017 experiment: stimulate a naive CD4+ T cell and watch the network
# resolve. One run per polarising condition, everything written to output/ and
# figures/.
#
#   Rscript scripts/01_stimulations.R
# ------------------------------------------------------------------------------

source("R/ginsim_engine.R")
source("R/model_io.R")
source("R/analysis.R")

NTIME <- 200   # as in original/2017_model/Procedures.txt ("200 time points")
NREP  <- 200   # "200 actualizations - GOOD!"
SEED  <- 2019

model      <- load_model(rules_file = "rules_corrected.txt")
conditions <- names(POLARISING_INPUTS)

runs <- list()
for (cond in conditions) {
  message("running ", cond, " ...")
  state <- stimulation_state(model$nodes, cond)
  runs[[cond]] <- ginsim_run(state, model$nodes, model$rules,
                             ntime = NTIME, nrep = NREP, seed = SEED)
}

dir.create("output", showWarnings = FALSE)

# --- trajectories, every node, every step ------------------------------------
traj <- do.call(rbind, Map(trajectory_table, runs, names(runs)))
write.csv(traj, "output/stimulation_trajectories.csv", row.names = FALSE)

# --- time-averaged activity per condition ------------------------------------
activity <- lapply(runs, mean_activity)
summary_tbl <- data.frame(node = model$nodes[, 2],
                          max_level = model$nodes[, 3],
                          stringsAsFactors = FALSE)
for (cond in conditions) summary_tbl[[cond]] <- round(activity[[cond]], 4)
write.csv(summary_tbl, "output/stimulation_activity.csv", row.names = FALSE)

# --- did anything settle? -----------------------------------------------------
settle <- data.frame(
  condition       = conditions,
  trajectories    = NREP,
  reached_steady  = vapply(runs, function(r) r$n_steady, integer(1)),
  distinct_states = vapply(runs, function(r) nrow(attractor_summary(r)), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(settle, "output/stimulation_attractors.csv", row.names = FALSE)
print(settle)

# --- figures ------------------------------------------------------------------
to_png("figures/trajectories_IL2.png",
       plot_trajectories(runs$IL2, c(MIRNAS, MASTER_TFS, "IL2", "NFKB1", "NFATC1"),
                         "TCR + CD28 stimulation: miRNAs and master regulators"))

to_png("figures/trajectories_IL2_cytokines.png",
       plot_trajectories(runs$IL2, CYTOKINES,
                         "TCR + CD28 stimulation: cytokine output"))

to_png("figures/master_tfs_by_condition.png",
       plot_condition_bars(activity, unname(MASTER_TFS),
                           "Master transcription factors across stimulations"))

to_png("figures/mirnas_by_condition.png",
       plot_condition_bars(activity, MIRNAS,
                           "miR-34c-5p and miR-155-5p across stimulations"),
       width = 1100)

message("done")
