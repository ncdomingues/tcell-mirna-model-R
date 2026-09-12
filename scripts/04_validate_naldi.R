# 04_validate_naldi.R ----------------------------------------------------------
#
# Does the network reproduce textbook T-helper biology at all? That question is
# separate from anything about miRNAs, and it has to be answered first.
#
# The Naldi et al. 2010 model -- one of the two published networks my model was
# built on -- was validated by its authors by showing that polarising cytokines
# drive the matching master transcription factor. Unlike my own table, it has
# real exogenous cytokine inputs (IL4_e, IL12_e, TGFB_e, ...), so the classic
# test can be run properly. This script runs it twice: on the base network, and
# with my two miRNAs added, to see whether the extension breaks the biology.
#
#   Rscript scripts/04_validate_naldi.R
# ------------------------------------------------------------------------------

source("R/ginsim_engine.R")
source("R/model_io.R")
source("R/analysis.R")

NTIME <- 200
NREP  <- 200
SEED  <- 2019

# nodes_corrected.txt raises IL4RA's ceiling from 1 to 2; the rules give it a
# level-2 target and other rules read `IL4RA:2`. See MODEL_NOTES.md.
nodes <- read_table_file("model/naldi2010/nodes_corrected.txt")
rules <- read_table_file("model/naldi2010/rules.txt")

NALDI_TFS <- c(Th1 = "TBX21", Th2 = "GATA3", Th17 = "RORC", iTreg = "FOXP3")

# Each condition is the standard in vitro polarising cocktail. The 2017 script
# (original/base_files/GinSimR_ND.R) built these as hand-counted position
# vectors and left IL6_e switched on in the iTreg one, which its own comment
# says should be "IL2 + TGFB + APC" -- naming the nodes avoids that class of
# slip entirely.
CONDITIONS <- list(
  "TCR + IL2"       = c("APC", "IL2_e"),
  "+ IL12 (Th1)"    = c("APC", "IL2_e", "IL12_e"),
  "+ IL4 (Th2)"     = c("APC", "IL2_e", "IL4_e"),
  "+ TGFB + IL6 (Th17)" = c("APC", "IL2_e", "TGFB_e", "IL6_e"),
  "+ TGFB (iTreg)"  = c("APC", "IL2_e", "TGFB_e")
)

EXPECTED <- c("TCR + IL2" = NA, "+ IL12 (Th1)" = "TBX21", "+ IL4 (Th2)" = "GATA3",
              "+ TGFB + IL6 (Th17)" = "RORC", "+ TGFB (iTreg)" = "FOXP3")

# Shared receptor chains and coactivators: unregulated, Input = 0, basal 1.
CONSTITUTIVE <- constitutive_nodes(nodes, rules)
message("constitutive components held on: ", paste(CONSTITUTIVE, collapse = ", "))

run_set <- function(rules, label) {
  out <- list()
  for (cond in names(CONDITIONS)) {
    message("  ", label, " / ", cond)
    state <- make_state(nodes, c(CONSTITUTIVE, CONDITIONS[[cond]]))
    out[[cond]] <- mean_activity(ginsim_run(state, nodes, rules,
                                            ntime = NTIME, nrep = NREP, seed = SEED))
  }
  out
}

# Four networks, so that a lost polarisation can be attributed to one miRNA
# rather than to "the extension" in general.
NETWORKS <- list(
  "base (no miRNA)"  = MIRNAS,
  "+ miR-34c-5p"     = "miR-155-5p",
  "+ miR-155-5p"     = "miR-34c-5p",
  "+ both"           = character(0)
)

sets <- list()
for (net in names(NETWORKS)) {
  message(net)
  sets[[net]] <- run_set(knock_out(rules, NETWORKS[[net]]), net)
}
base <- sets[["base (no miRNA)"]]
ext  <- sets[["+ both"]]

dir.create("output", showWarnings = FALSE)

build <- function(activity) {
  tbl <- data.frame(condition = names(CONDITIONS), stringsAsFactors = FALSE)
  for (tf in NALDI_TFS) {
    tbl[[tf]] <- round(vapply(activity, function(a) a[[tf]], 0), 4)
  }
  tbl$expected <- unname(EXPECTED[tbl$condition])
  baseline <- activity[["TCR + IL2"]]
  tbl$verdict <- vapply(seq_len(nrow(tbl)), function(i) {
    want <- tbl$expected[i]
    if (is.na(want)) return("baseline")
    got  <- vapply(activity[[i]][NALDI_TFS], identity, 0)
    rose <- got[[want]] > baseline[[want]] + 0.05
    top  <- names(which.max(got - baseline[NALDI_TFS])) == want
    if (rose && top) "pass" else if (rose) "partial" else "fail"
  }, "")
  tbl
}

tables <- lapply(names(sets), function(net) {
  t <- build(sets[[net]]); t$network <- net; t
})
result <- do.call(rbind, tables)
result <- result[, c("network", "condition", unname(NALDI_TFS), "expected", "verdict")]

write.csv(result, "output/naldi_polarisation.csv", row.names = FALSE)
print(result, row.names = FALSE)

cat("\npolarisations recovered:\n")
for (net in names(sets)) {
  n <- sum(result$verdict[result$network == net] == "pass")
  cat(sprintf("  %-18s %d / 4\n", net, n))
}

# Which cocktail breaks, and under which miRNA?
lost <- result[result$verdict == "fail", c("network", "condition", "expected")]
if (nrow(lost)) {
  cat("\nlost polarisations:\n"); print(lost, row.names = FALSE)
} else {
  cat("\nno polarisation lost in any network\n")
}

# --- trace the lost Th1 polarisation back to a single edge --------------------
#
# TBX21 in this network is `(TBX21 | STAT1) & !GATA3` -- it hangs off STAT1, not
# STAT4, so supplying IL12 is not by itself enough. With IFNB_e and IL27_e off,
# the only route to STAT1 is IFNGR, and IFNGR carries a `!miR-155-5p` term.
# Print the chain so the claim in VALIDATION.md is checkable, not asserted.
chain <- c("miR-155-5p", "IL12R", "STAT4", "IFNG", "IFNGR", "STAT1", "TBX21")
th1 <- data.frame(node = chain, stringsAsFactors = FALSE)
for (net in names(sets)) th1[[net]] <- round(sets[[net]][["+ IL12 (Th1)"]][chain], 3)
cat("\nIL12 condition, the route to TBX21:\n")
print(th1, row.names = FALSE)
write.csv(th1, "output/naldi_th1_chain.csv", row.names = FALSE)

to_png("figures/naldi_polarisation.png",
       plot_condition_bars(ext, unname(NALDI_TFS),
                           "Naldi 2010 + miRNAs: master regulators by polarising cocktail"),
       width = 1500)

message("done")
