# test_engine.R ----------------------------------------------------------------
#
# Does the engine implement the notation it claims to, and are the digitized
# tables structurally what the thesis says they are?
#
#   Rscript tests/test_engine.R
#
# No testthat, no packages -- a bare R install runs this.
# ------------------------------------------------------------------------------

source("R/ginsim_engine.R")
source("R/model_io.R")
source("R/analysis.R")

.passed <- 0L
.failed <- character(0)

ok <- function(label, expr) {
  result <- tryCatch(isTRUE(expr), error = function(e) conditionMessage(e))
  if (isTRUE(result)) {
    .passed <<- .passed + 1L
    cat("  ok   ", label, "\n")
  } else {
    .failed <<- c(.failed, label)
    cat("  FAIL ", label,
        if (is.character(result)) paste0(" (", result, ")") else "", "\n")
  }
}

throws <- function(expr) {
  inherits(try(force(expr), silent = TRUE), "try-error")
}

# A three-node toy model: IN is an input, MID is two-level, OUT reads a threshold.
toy_nodes <- data.frame(node = 1:3, name = c("IN", "MID", "OUT"),
                        max_level = c(1L, 2L, 1L), input = c(1L, 0L, 0L),
                        stringsAsFactors = FALSE)
toy_rules <- data.frame(target = c(2L, 2L, 3L), name = c("MID", "MID", "OUT"),
                        value = c(2L, 1L, 1L),
                        formula = c("IN", "IN", "MID:2"),
                        stringsAsFactors = FALSE)

evaluate <- function(nodes, rules, state) {
  compiled <- compile_rules(nodes, rules)
  env <- list2env(list(S = state), parent = baseenv())
  vapply(compiled, eval, logical(1), envir = env, USE.NAMES = FALSE)
}


cat("\ncompile_rules -- name resolution\n")

ok("bare name is true at level >= 1",
   evaluate(toy_nodes, toy_rules, c(1L, 0L, 0L))[1])

ok("bare name is false at level 0",
   !evaluate(toy_nodes, toy_rules, c(0L, 0L, 0L))[1])

ok("NODE:2 is false when the node sits at level 1",
   !evaluate(toy_nodes, toy_rules, c(1L, 1L, 0L))[3])

ok("NODE:2 is true when the node sits at level 2",
   evaluate(toy_nodes, toy_rules, c(1L, 2L, 0L))[3])

# This is the 2017 engine's bug: it never substituted `:1`, so `MID:1` survived
# as `0:1` -- an R integer sequence, silently coerced to FALSE.
ok("NODE:1 means the same as the bare name (the 2017 `0:1` bug)",
   {
     r <- toy_rules; r$formula[3] <- "MID:1"
     evaluate(toy_nodes, r, c(1L, 1L, 0L))[3] &&
       !evaluate(toy_nodes, r, c(1L, 0L, 0L))[3]
   })

ok("longer names are substituted before the names they contain",
   {
     n <- data.frame(node = 1:2, name = c("IL2", "IL2R"), max_level = c(1L, 1L),
                     input = c(1L, 0L), stringsAsFactors = FALSE)
     r <- data.frame(target = 2L, name = "IL2R", value = 1L,
                     formula = "IL2R & !IL2", stringsAsFactors = FALSE)
     # If `IL2` were substituted first it would corrupt `IL2R` into `1R`/`0R`.
     evaluate(n, r, c(0L, 1L)) && !evaluate(n, r, c(1L, 1L))
   })

ok("a missing operator between two names is caught, not parsed as a call",
   throws({
     n <- data.frame(node = 1:2, name = c("TGFBR", "IL6R"), max_level = c(1L, 1L),
                     input = c(1L, 1L), stringsAsFactors = FALSE)
     r <- data.frame(target = 2L, name = "IL6R", value = 1L,
                     formula = "TGFBRIL6R", stringsAsFactors = FALSE)
     compile_rules(n, r)
   }))

ok("hyphenated names resolve (miR-34c-5p, CDKN1A-p21)",
   {
     n <- data.frame(node = 1:2, name = c("miR-34c-5p", "CDKN1A-p21"),
                     max_level = c(1L, 1L), input = c(1L, 0L),
                     stringsAsFactors = FALSE)
     r <- data.frame(target = 2L, name = "CDKN1A-p21", value = 1L,
                     formula = "!miR-34c-5p", stringsAsFactors = FALSE)
     evaluate(n, r, c(0L, 0L)) && !evaluate(n, r, c(1L, 0L))
   })

ok("TRUE is accepted as a constitutive formula",
   {
     n <- data.frame(node = 1L, name = "X", max_level = 1L, input = 0L,
                     stringsAsFactors = FALSE)
     r <- data.frame(target = 1L, name = "X", value = 1L, formula = "TRUE",
                     stringsAsFactors = FALSE)
     evaluate(n, r, 0L)
   })

ok("an undeclared name in a formula is an error, not a silent FALSE",
   throws({
     r <- toy_rules; r$formula[3] <- "GHOST"
     compile_rules(toy_nodes, r)
   }))


cat("\ncompile_rules -- operator semantics\n")

ops_nodes <- data.frame(node = 1:3, name = c("A", "B", "C"),
                        max_level = rep(1L, 3), input = rep(1L, 3),
                        stringsAsFactors = FALSE)
ops <- function(formula, state) {
  r <- data.frame(target = 1L, name = "A", value = 1L, formula = formula,
                  stringsAsFactors = FALSE)
  evaluate(ops_nodes, r, state)
}

ok("& is AND",  ops("A & B", c(1L, 1L, 0L)) && !ops("A & B", c(1L, 0L, 0L)))
ok("| is OR",   ops("A | B", c(0L, 1L, 0L)) && !ops("A | B", c(0L, 0L, 0L)))
ok("! is NOT",  ops("!A", c(0L, 0L, 0L))    && !ops("!A", c(1L, 0L, 0L)))
ok("NOT binds tighter than AND", !ops("!A & B", c(1L, 1L, 0L)))
ok("AND binds tighter than OR",   ops("A | B & C", c(1L, 0L, 0L)))
ok("parentheses override precedence",
   !ops("!(A | B)", c(1L, 0L, 0L)) && ops("!(A | B)", c(0L, 0L, 0L)))


cat("\nupdate semantics\n")

ok("a node takes the highest value among its satisfied rules",
   {
     run <- ginsim_run(c(1L, 0L, 0L), toy_nodes, toy_rules, 20, 1, seed = 1)
     run$evol["MID", ncol(run$evol)] == 2
   })

ok("a node with no satisfied rule falls to 0",
   {
     run <- ginsim_run(c(0L, 2L, 0L), toy_nodes, toy_rules, 20, 1, seed = 1)
     run$evol["MID", ncol(run$evol)] == 0
   })

ok("input nodes -- no rule targets them -- hold their value",
   {
     run <- ginsim_run(c(1L, 0L, 0L), toy_nodes, toy_rules, 20, 1, seed = 1)
     all(run$evol["IN", ] == 1)
   })

ok("updates move one level at a time, never jumping 0 -> 2",
   {
     compiled <- compile_rules(toy_nodes, toy_rules)
     targets  <- lapply(1:3, function(i) which(toy_rules[, 1] == i))
     nxt <- ginsim_next(c(1L, 0L, 0L), toy_rules, compiled, targets)
     nxt[2] == 1L
   })

ok("exactly one node changes per step",
   {
     compiled <- compile_rules(toy_nodes, toy_rules)
     targets  <- lapply(1:3, function(i) which(toy_rules[, 1] == i))
     set.seed(3)
     st <- c(1L, 0L, 0L)
     all(vapply(1:8, function(i) {
       nxt <- ginsim_next(st, toy_rules, compiled, targets)
       d <- sum(nxt != st); st <<- nxt; d <= 1L
     }, logical(1)))
   })

ok("a steady state is recorded once reached",
   {
     run <- ginsim_run(c(1L, 0L, 0L), toy_nodes, toy_rules, 30, 5, seed = 1)
     run$n_steady == 5 && all(run$ss[, "MID"] == 2)
   })

ok("the same seed reproduces the same ensemble",
   {
     a <- ginsim_run(c(1L, 0L, 0L), toy_nodes, toy_rules, 15, 10, seed = 99)
     b <- ginsim_run(c(1L, 0L, 0L), toy_nodes, toy_rules, 15, 10, seed = 99)
     identical(a$evol, b$evol)
   })


cat("\ntable validation\n")

ok("a rule targeting a node index that does not exist is rejected",
   throws({ r <- toy_rules; r$target[1] <- 99L; check_model(toy_nodes, r) }))

ok("a rule value above the node's ceiling is rejected under strict",
   throws({ r <- toy_rules; r$value[1] <- 5L; check_model(toy_nodes, r) }))

ok("the same table passes with strict = FALSE, so 2017 runs can be reproduced",
   {
     r <- toy_rules; r$value[1] <- 5L
     suppressWarnings(check_model(toy_nodes, r, strict = FALSE))
   })

ok("a state vector of the wrong length is rejected",
   throws(check_model(toy_nodes, toy_rules, state = c(1L, 0L))))

ok("a rule whose Name disagrees with the node it targets warns",
   {
     r <- toy_rules; r$name[1] <- "ELSEWHERE"
     w <- FALSE
     withCallingHandlers(check_model(toy_nodes, r),
                         warning = function(x) { w <<- TRUE; invokeRestart("muffleWarning") })
     w
   })


cat("\nmodel/ tables\n")

published <- load_model(rules_file = "rules.txt")
corrected <- load_model(rules_file = "rules_corrected.txt")

ok("the model has 91 nodes",  nrow(published$nodes) == 91)
ok("the model has 86 rules",  nrow(published$rules) == 86)
ok("both miRNAs are rule-governed, not inputs",
   all(MIRNAS %in% published$rules[, 2]))
ok("every master transcription factor is present",
   all(MASTER_TFS %in% published$nodes[, 2]))

ok("the published table is rejected under strict validation (IL4R -> node 18)",
   throws(suppressWarnings(check_model(published$nodes, published$rules))))

ok("the corrected table passes strict validation",
   suppressWarnings(check_model(corrected$nodes, corrected$rules)))

ok("the correction moves exactly two rows, and nothing else changes",
   {
     diffs <- which(published$rules$target != corrected$rules$target)
     identical(published$rules[-diffs, ], corrected$rules[-diffs, ]) &&
       length(diffs) == 2 &&
       all(published$rules$target[diffs] == 18) &&
       all(corrected$rules$target[diffs] == 66) &&
       all(published$rules$name[diffs] == "IL4R")
   })

variant_files  <- list.files("model/mir34c_variants", pattern = "[.]txt$",
                             full.names = TRUE)
variant_tables <- lapply(variant_files, function(f)
  repair_2017_table(read_table_file(f), quiet = TRUE))
names(variant_tables) <- sub("[.]txt$", "", basename(variant_files))

ok("there are twelve miR-34c-5p promoter variants", length(variant_tables) == 12)

ok("every variant compiles once repaired",
   all(vapply(variant_tables, function(r) {
     n <- published$nodes
     for (i in seq_len(nrow(r))) n[r[i, 1], 3] <- max(n[r[i, 1], 3], r[i, 3])
     length(compile_rules(n, r)) == nrow(r)
   }, logical(1))))

ok("every variant proposes a different formula for node 84",
   {
     f84 <- vapply(variant_tables, function(r)
       paste(r[r[, 1] == 84, 3], r[r[, 1] == 84, 4], collapse = "/"), "")
     length(unique(f84)) == length(f84)
   })

# Documented in MODEL_NOTES.md: the twelve are three backbone generations, not
# twelve edits of one table. If that ever stops being true, the grouping in
# scripts/03 is reading a different set of files than it was written for.
ok("the variants fall into exactly three backbone generations",
   {
     backbone <- vapply(variant_tables, function(r) {
       o <- r[r[, 1] != 84, ]
       paste(o[[1]], o[[3]], o[[4]], collapse = "|")
     }, "")
     length(unique(backbone)) == 3
   })

ok("none of the three backbones matches the published February table",
   {
     pub <- repair_2017_table(published$rules, quiet = TRUE)
     pub <- pub[pub[, 1] != 84, ]
     sig <- paste(pub[[1]], pub[[3]], pub[[4]], collapse = "|")
     backbone <- vapply(variant_tables, function(r) {
       o <- r[r[, 1] != 84, ]
       paste(o[[1]], o[[3]], o[[4]], collapse = "|")
     }, "")
     !(sig %in% backbone)
   })

ok("repair_2017_table fixes both known slips and leaves everything else alone",
   {
     before <- read_table_file("model/mir34c_variants/MYC_TP53_FOXO3.txt")
     after  <- repair_2017_table(before, quiet = TRUE)
     moved  <- which(before[[1]] != after[[1]])
     edited <- which(before[[4]] != after[[4]])
     length(moved) == 2 && all(before[moved, 2] == "IL4R") &&
       all(after[moved, 1] == 66) &&
       length(edited) == 1 && grepl("TGFBR | IL6R", after[edited, 4], fixed = TRUE) &&
       identical(before[-c(moved, edited), -1], after[-c(moved, edited), -1])
   })

ok("repair_2017_table is a no-op on an already-clean table",
   identical(corrected$rules, repair_2017_table(corrected$rules, quiet = TRUE)))

ok("the Naldi 2010 and Abou-Jaoude 2015 base tables compile",
   {
     naldi <- list(nodes = read_table_file("model/naldi2010/nodes_corrected.txt"),
                   rules = read_table_file("model/naldi2010/rules.txt"))
     aj    <- list(nodes = read_table_file("model/abou_jaoude2015/nodes.txt"),
                   rules = read_table_file("model/abou_jaoude2015/rules_corrected.txt"))
     length(compile_rules(naldi$nodes, naldi$rules)) == nrow(naldi$rules) &&
       length(compile_rules(aj$nodes, aj$rules)) == nrow(aj$rules)
   })


cat("\nhelpers\n")

ok("node_id maps names to row indices",
   all(node_id(published$nodes, c("TCR", "miR-155-5p")) == c(2L, 85L)))

ok("node_id rejects an unknown name", throws(node_id(published$nodes, "GHOST")))

ok("make_state sets only the nodes it is given",
   {
     s <- make_state(published$nodes, c("TCR", "CD28"))
     sum(s) == 2 && s[node_id(published$nodes, "TCR")] == 1
   })

ok("make_state accepts explicit levels",
   make_state(published$nodes, c(IL2 = 2))[34] == 2)

ok("knock_out removes a miRNA's own rules and nothing else",
   {
     r <- knock_out(corrected$rules, "miR-155-5p")
     nrow(r) == nrow(corrected$rules) - 1 && !("miR-155-5p" %in% r[, 2])
   })

ok("a knocked-out miRNA stays at 0 for the whole run",
   {
     r   <- knock_out(corrected$rules, MIRNAS)
     run <- ginsim_run(stimulation_state(corrected$nodes, "IL2"),
                       corrected$nodes, r, 40, 5, seed = 4)
     all(run$evol[MIRNAS, ] == 0)
   })

ok("trajectory returns a vector for one node and a matrix for several",
   {
     run <- ginsim_run(stimulation_state(corrected$nodes, "IL2"),
                       corrected$nodes, corrected$rules, 20, 3, seed = 5)
     is.null(dim(trajectory(run, "IL2"))) &&
       nrow(trajectory(run, MIRNAS)) == 2
   })

ok("mean_activity normalises multi-valued nodes onto 0-1",
   {
     run <- ginsim_run(stimulation_state(corrected$nodes, "IL2"),
                       corrected$nodes, corrected$rules, 60, 5, seed = 6)
     a <- mean_activity(run)
     all(a >= 0) && all(a <= 1)
   })

ok("stimulation_state rejects an unknown condition",
   throws(stimulation_state(corrected$nodes, "Th42")))


cat("\n", .passed, " passed, ", length(.failed), " failed\n", sep = "")
if (length(.failed)) {
  cat("failed:\n"); cat(paste0("  - ", .failed, collapse = "\n"), "\n")
  quit(status = 1)
}
