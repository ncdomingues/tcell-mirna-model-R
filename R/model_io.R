# model_io.R -------------------------------------------------------------------
#
# Loading the node/rule tables, building initial states by node NAME rather than
# by the hand-counted position vectors of the 2017 scripts, and switching miRNA
# regulation on and off.
# ------------------------------------------------------------------------------


#' Read a node or rule table
#'
#' The 2017 tables were written by hand and carry trailing tabs on a few rows,
#' which `read.delim` turns into a spurious fifth column. Keep the first four.
read_table_file <- function(path) {
  tbl <- read.delim(path, header = TRUE, sep = "\t", quote = "\"", dec = ".",
                    fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
  tbl <- tbl[, 1:4, drop = FALSE]
  tbl[[2]] <- trimws(tbl[[2]])
  tbl[[4]] <- trimws(tbl[[4]])
  tbl
}


#' Load a model: its node table plus one of its rule tables
#'
#' @param dir        directory holding `nodes.txt`
#' @param rules_file rule table filename, relative to `dir`
load_model <- function(dir = "model", rules_file = "rules.txt") {
  nodes <- read_table_file(file.path(dir, "nodes.txt"))
  rules <- read_table_file(file.path(dir, rules_file))
  names(nodes) <- c("node", "name", "max_level", "input")
  names(rules) <- c("target", "name", "value", "formula")
  list(nodes = nodes, rules = rules)
}


#' Index of a node, by name
node_id <- function(nodes, name) {
  idx <- match(name, nodes[, 2])
  if (anyNA(idx)) stop("unknown node(s): ", paste(name[is.na(idx)], collapse = ", "))
  idx
}


#' Mean trajectory of one or more nodes from a run
#'
#' Replaces `run$evol[85, ]` with `trajectory(run, "miR-155-5p")`. The magic
#' row numbers in the 2017 scripts were a standing source of mislabelled plots.
trajectory <- function(run, name) {
  run$evol[node_id(run$nodes, name), , drop = length(name) == 1L]
}


#' Build an initial state vector
#'
#' @param nodes the node table
#' @param on    node names to set to 1 (or to a named level, e.g. c(IL2 = 2))
#'
#' Every other node starts at 0.
make_state <- function(nodes, on = character(0)) {
  state <- integer(nrow(nodes))
  if (length(on) == 0) return(state)

  if (is.null(names(on))) {
    state[node_id(nodes, on)] <- 1L
  } else {
    state[node_id(nodes, names(on))] <- as.integer(on)
  }
  state
}


# The resting naive CD4+ T cell as stimulated in the thesis: anti-CD3 +
# anti-CD28 with APC-presented antigen, plus the model's unregulated
# co-factors and the constitutively present miRNA inputs. Both CD45 isoforms
# are on, matching `original/2017_model/Procedures.txt`.
BASE_STIMULUS <- c(
  "APC-Antigen", "TCR", "CD45RA", "CD45RO", "CD28", "CD3",
  "BCMcomplex",                          # BCL10 + CARD11 + MALT1
  "FOXP1", "CREBBP", "KAT2B", "CGC",     # unregulated co-factors
  "miR-181a-5p", "miR-21-5p", "miR-146b-5p", "let-7a-5p", "miR-101-3p"
)

# Polarising cytokines added on top of the base stimulus.
#
# `IL2`, `Th1` and `Th2` are exactly the three conditions scripted in 2017
# (`original/stimulations/`). `Th17` and `iTreg` are additions of mine: this
# model has no exogenous TGFB or IL6 input node, so they drive the STAT3 axis
# through the two input cytokines that do exist. Treat them as exploratory.
POLARISING_INPUTS <- list(
  IL2   = character(0),      # TCR + CD28 only; IL2 itself is rule-driven
  Th1   = "IL12",            # 2017: STATE0[37] = 1
  Th2   = c("IL12", "IL4"),  # 2017: STATE0[37] = 1; STATE0[42] = 1
  Th17  = c("IL1", "IL23"),  # added here -- STAT3 -> RORC via the IL1/IL23 inputs
  iTreg = "IL23"             # added here
)

# IL4 (node 42) is itself rule-governed, so seeding it in the Th2 condition only
# sets its starting level -- the rules take it from there. That is what the 2017
# script did and it is kept, but it is a weaker perturbation than a true input.

#' Initial state for a named stimulation condition
#'
#' @param condition one of names(POLARISING_INPUTS)
stimulation_state <- function(nodes, condition = "IL2") {
  extra <- POLARISING_INPUTS[[condition]]
  if (is.null(extra)) {
    stop("unknown condition '", condition, "'; expected one of ",
         paste(names(POLARISING_INPUTS), collapse = ", "))
  }
  make_state(nodes, c(BASE_STIMULUS, extra))
}


#' Nodes that no rule targets
#'
#' These hold whatever value the initial state gives them, forever. They come in
#' two kinds and the node table does not always tell them apart:
#'
#'  - experimental inputs -- the antigen and the exogenous cytokines, which are
#'    off unless the condition under test switches them on;
#'  - constitutive components -- shared receptor chains (CGC, IFNGR1, GP130) and
#'    coactivators (EP300, CREBBP, KAT2B), which are simply always present.
#'
#' `original/base_files/method_modelling.txt` states the convention for the
#' second kind: "Do not include inputs or unregulated nodes (with no incoming
#' edges, normally basal value = 1)". In the Naldi table those carry Input = 0,
#' so they are derivable; in my own table everything unregulated carries
#' Input = 1 and the two kinds have to be listed by hand (see BASE_STIMULUS).
unregulated_nodes <- function(nodes, rules) {
  nodes[!(nodes[, 1] %in% rules[, 1]), 2]
}


#' Constitutive components: unregulated, and not declared an experimental input
constitutive_nodes <- function(nodes, rules) {
  free <- !(nodes[, 1] %in% rules[, 1])
  nodes[free & nodes[, 4] == 0, 2]
}


#' Apply the two known transcription errors' fixes to a 2017 rule table
#'
#' Both are typing slips in tables that were maintained by hand in a text
#' editor, and both are recorded in MODEL_NOTES.md. Nothing else is touched, and
#' the function reports what it changed so a repair is never silent.
#'
#'  1. Two `IL4R` rows carry target index 18 (FOXO1) instead of 66 (IL4R). The
#'     engine targets by index, so in 2017 those rules drove FOXO1 -- and pushed
#'     it to level 2, past its declared ceiling of 1.
#'  2. Three variant tables read `TGFBRIL6R`, two node names with the operator
#'     between them missing. The neighbouring terms are all OR'd receptors, so
#'     the missing character is a `|`; only that one character is inserted.
#'
#' @param quiet suppress the message listing what was repaired
repair_2017_table <- function(rules, quiet = FALSE) {
  changes <- character(0)

  mis <- which(rules[, 1] == 18 & trimws(rules[, 2]) == "IL4R")
  if (length(mis)) {
    rules[mis, 1] <- 66
    changes <- c(changes, paste0(length(mis), " IL4R rule(s) retargeted 18 -> 66"))
  }

  typo <- grep("TGFBRIL6R", rules[, 4], fixed = TRUE)
  if (length(typo)) {
    rules[typo, 4] <- gsub("TGFBRIL6R", "TGFBR | IL6R", rules[typo, 4], fixed = TRUE)
    changes <- c(changes, paste0(length(typo), " 'TGFBRIL6R' -> 'TGFBR | IL6R'"))
  }

  if (!quiet && length(changes)) message("  repaired: ", paste(changes, collapse = "; "))
  rules
}


#' Remove a node's transcriptional rules, knocking it out of the network
#'
#' Used to build the four miRNA conditions: dropping a miRNA's own rules leaves
#' it permanently at its initial value (0), so every `!miR-...` term in the rest
#' of the table evaluates as if that miRNA were absent.
knock_out <- function(rules, name) {
  rules[!(trimws(rules[, 2]) %in% name), , drop = FALSE]
}


#' Activity at the final step
#'
#' Mean level at the last time point, normalised by the node's maximum level so
#' Boolean and multi-valued nodes are comparable.
final_activity <- function(run) {
  final <- run$evol[, ncol(run$evol)]
  stats::setNames(as.numeric(final) / run$nodes[, 3], run$nodes[, 2])
}


#' Time-averaged activity over the tail of the run
#'
#' Under TCR + IL2 stimulation this network has no fixed point -- the
#' `TCR = ... & !IL2:2` rule makes accumulating IL2 shut proximal signalling
#' down, which then restarts, so trajectories settle into a cycle rather than a
#' steady state (see MODEL_NOTES.md). A single final step is therefore a
#' snapshot at an arbitrary phase of that cycle. Averaging the ensemble mean
#' over the tail of the run measures occupancy instead, which is stable.
#'
#' The default window is the second half of the run. Under the standard 200-step
#' protocol the transient is over by about step 100 and the cycle period is
#' roughly 60 steps, so that window covers the settled behaviour across more
#' than one full cycle and is not sensitive to where the run happens to stop.
#'
#' @param tail fraction of the run to average over, from the end
mean_activity <- function(run, tail = 0.5) {
  n     <- ncol(run$evol)
  from  <- max(1L, n - as.integer(round(n * tail)) + 1L)
  level <- rowMeans(run$evol[, from:n, drop = FALSE])
  stats::setNames(as.numeric(level) / run$nodes[, 3], run$nodes[, 2])
}
