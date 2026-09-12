# 03_mir34c_tf_variants.R ------------------------------------------------------
#
# miR-34c-5p's promoter was the open question of the thesis. Which transcription
# factors actually drive it? Over February and March 2017 I saved twelve rule
# tables proposing different answers -- GATA3 alone, MYC alone, TP53 + FOXO3,
# and combinations -- and ran the model under each.
#
# One thing worth knowing before reading the numbers: the twelve tables are not
# twelve one-line edits of a single network. They fall into three backbone
# generations that also differ from each other in eight other rules (IKBKG,
# IL4, STAT3, TGFBR, GATA3, BCL6, SPI1-PU1 and miR-155-5p itself), because the
# rest of the model was still being revised in parallel. Comparing a March
# hypothesis against a February one therefore confounds the promoter question
# with those revisions. The script prints the backbone group for every row, and
# comparisons are only meaningful within a group.
#
#   Rscript scripts/03_mir34c_tf_variants.R
# ------------------------------------------------------------------------------

source("R/ginsim_engine.R")
source("R/model_io.R")
source("R/analysis.R")

NTIME <- 200
NREP  <- 150
SEED  <- 2019

MIR34C <- 84   # the node whose formula each variant proposes

nodes <- read_table_file("model/nodes.txt")
state <- stimulation_state(nodes, "IL2")

variants <- list.files("model/mir34c_variants", pattern = "[.]txt$", full.names = TRUE)
labels   <- sub("[.]txt$", "", basename(variants))

tables <- lapply(variants, function(f) repair_2017_table(read_table_file(f), quiet = TRUE))
names(tables) <- labels

# The published February table (model/rules.txt), for reference.
tables[["published"]] <- repair_2017_table(read_table_file("model/rules.txt"), quiet = TRUE)

# --- group by everything except the miR-34c-5p rule ---------------------------
backbone <- vapply(tables, function(t) {
  other <- t[t[, 1] != MIR34C, ]
  paste(other[[1]], other[[3]], other[[4]], collapse = "|")
}, "")
group <- match(backbone, unique(backbone))
names(group) <- names(tables)

message(length(unique(group)), " distinct backbones across ",
        length(tables), " tables:")
for (g in sort(unique(group))) {
  message("  backbone ", g, ": ", paste(names(group)[group == g], collapse = ", "))
}

# --- run them all -------------------------------------------------------------
runs <- list()
for (label in names(tables)) {
  message("running ", label, " ...")
  rules <- tables[[label]]

  # The 6_3_17 variant treats miR-34c-5p as two-level (TP53 + FOXO3 gives strong
  # induction). Raise each node's ceiling to whatever its own rules ask for.
  variant_nodes <- nodes
  for (i in seq_len(nrow(rules))) {
    tgt <- rules[i, 1]
    variant_nodes[tgt, 3] <- max(variant_nodes[tgt, 3], rules[i, 3])
  }

  runs[[label]] <- ginsim_run(state, variant_nodes, rules,
                              ntime = NTIME, nrep = NREP, seed = SEED)
}

activity <- lapply(runs, mean_activity)

dir.create("output", showWarnings = FALSE)

result <- data.frame(
  variant     = names(tables),
  backbone    = unname(group),
  mir34c_rule = vapply(tables, function(t)
    paste(t[t[, 1] == MIR34C, 3], t[t[, 1] == MIR34C, 4],
          sep = ": ", collapse = " / "), ""),
  stringsAsFactors = FALSE
)
result$`miR-34c-5p` <- round(vapply(activity, function(a) a[["miR-34c-5p"]], 0), 4)
for (tf in unname(MASTER_TFS)) {
  result[[tf]] <- round(vapply(activity, function(a) a[[tf]], 0), 4)
}
result <- result[order(result$backbone, -result$`miR-34c-5p`), ]

write.csv(result, "output/mir34c_tf_variants.csv", row.names = FALSE)
print(result[, c("variant", "backbone", "miR-34c-5p", unname(MASTER_TFS))],
      row.names = FALSE)

# --- figure: within-backbone comparison ---------------------------------------
to_png("figures/mir34c_tf_variants.png",
       {
         ord   <- order(result$backbone, result$`miR-34c-5p`)
         shade <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
         op <- par(mar = c(4.5, 12, 3.5, 2))
         on.exit(par(op), add = TRUE)
         barplot(result$`miR-34c-5p`[ord],
                 names.arg = result$variant[ord], horiz = TRUE, las = 1,
                 col = shade[result$backbone[ord]], border = NA,
                 xlim = c(0, 1), cex.names = 0.8,
                 xlab = "miR-34c-5p occupancy over the second half of the run",
                 main = "Twelve hypotheses for the miR-34c-5p promoter")
         legend("bottomright", bty = "n", cex = 0.8, fill = shade[sort(unique(group))],
                border = NA, legend = paste("backbone", sort(unique(group))))
       },
       width = 1500, height = 950)

message("done")
