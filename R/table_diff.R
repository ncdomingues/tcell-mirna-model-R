# table_diff.R -----------------------------------------------------------------
#
# Rule-level diff between two rule tables -- specifically the 2017 working table
# and the thesis's Table 5.
#
# Two problems make a plain text diff useless here:
#
#  1. The write-up generalised some node names (`NFKB1` -> `NFKB`, `NFATC1` ->
#     `NFAT`, `AKT1` -> `AKT1-PKB`). A rename is not a rewritten rule.
#  2. The same logic can be written several ways. `CGC & IL4` on one row plus
#     `IL4` on another is the same function as `IL4 | (CGC & IL4)` on one row,
#     because the engine ORs the alternatives for a (node, value) pair.
#
# So formulas are compared by behaviour, not by text: both are evaluated over a
# large random sample of states and reported as equivalent only if they never
# disagree. That is a practical test, not a proof -- it can miss a difference
# that only shows on states the sample never visits -- so the sample size is
# reported alongside the verdict.
# ------------------------------------------------------------------------------


# Working-table name -> Table 5 name. Applied to the working table only.
NAME_ALIASES <- c(
  "NFATC1" = "NFAT",
  "NFKB1"  = "NFKB",
  "AKT1"   = "AKT1-PKB"
)


#' Rename nodes in a formula
#'
#' Token-aware, not plain substitution: `AKT1` is a prefix of `AKT1-PKB`, so a
#' `gsub` would turn an already-renamed formula into `AKT1-PKB-PKB`. Identifiers
#' are matched whole and mapped only on an exact hit.
rename_nodes <- function(f, aliases) {
  if (length(aliases) == 0) return(f)
  m <- gregexpr("[A-Za-z][A-Za-z0-9_.-]*", f)
  regmatches(f, m) <- lapply(regmatches(f, m), function(tok) {
    hit <- match(tok, names(aliases))
    ifelse(is.na(hit), tok, aliases[hit])
  })
  f
}


#' Compile a formula over a shared vocabulary
#'
#' @param vocab named integer vector: node name -> index into the state vector
compile_over <- function(formula, vocab) {
  text <- formula
  for (nm in names(vocab)[order(nchar(names(vocab)), decreasing = TRUE)]) {
    idx <- vocab[[nm]]
    for (k in 3:1) {
      text <- gsub(paste0(nm, ":", k), paste0("<", idx, "#", k, ">"), text, fixed = TRUE)
    }
    text <- gsub(nm, paste0("<", idx, "#1>"), text, fixed = TRUE)
  }
  if (grepl("[A-Za-z]", gsub("<[0-9]+#[0-9]+>|TRUE|FALSE", "", text))) {
    stop("unresolved name in '", formula, "'")
  }
  parse(text = gsub("<([0-9]+)#([0-9]+)>", "(S[\\1] >= \\2)", text))[[1]]
}


#' Do two formulas agree on every sampled state?
formulas_agree <- function(f_a, f_b, vocab, levels, n = 5000, seed = 1) {
  ea <- compile_over(f_a, vocab)
  eb <- compile_over(f_b, vocab)
  set.seed(seed)
  for (i in seq_len(n)) {
    S <- as.integer(round(runif(length(levels), 0, levels)))
    env <- list2env(list(S = S), parent = baseenv())
    if (!identical(eval(ea, env), eval(eb, env))) return(FALSE)
  }
  TRUE
}


#' Compare two rule tables node by node
#'
#' @return one row per (node, value) present in either table: the formula from
#'   each, and a verdict -- identical, equivalent, differs, only in a, only in b
diff_tables <- function(nodes_a, rules_a, nodes_b, rules_b,
                        label_a = "a", label_b = "b", n_sample = 5000) {
  collapse <- function(nodes, rules, aliases) {
    name <- trimws(nodes[rules[, 1], 2])
    if (length(aliases)) {
      hit  <- match(name, names(aliases))
      name <- ifelse(is.na(hit), name, aliases[hit])
    }
    f  <- gsub("[[:space:]]+", "", rename_nodes(rules[, 4], aliases))
    id <- paste(name, rules[, 3], sep = "@")
    # Alternatives for one (node, value) are OR'd by the engine.
    tapply(f, id, function(x) paste0("(", paste(sort(unique(x)), collapse = ")|("), ")"))
  }

  ca <- collapse(nodes_a, rules_a, NAME_ALIASES)
  cb <- collapse(nodes_b, rules_b, character(0))

  # Shared vocabulary: every node either table declares, after renaming.
  names_a <- rename_nodes(trimws(nodes_a[, 2]), NAME_ALIASES)
  lv <- c(stats::setNames(as.integer(nodes_b[, 3]), trimws(nodes_b[, 2])),
          stats::setNames(as.integer(nodes_a[, 3]), names_a))
  lv <- tapply(lv, names(lv), max)
  vocab  <- stats::setNames(seq_along(lv), names(lv))
  levels <- as.integer(lv)

  ids <- sort(union(names(ca), names(cb)))
  out <- data.frame(node = sub("@.*", "", ids),
                    value = as.integer(sub(".*@", "", ids)),
                    stringsAsFactors = FALSE)
  out[[label_a]] <- unname(ifelse(ids %in% names(ca), ca[ids], NA_character_))
  out[[label_b]] <- unname(ifelse(ids %in% names(cb), cb[ids], NA_character_))

  out$verdict <- vapply(seq_along(ids), function(i) {
    fa <- out[[label_a]][i]; fb <- out[[label_b]][i]
    if (is.na(fa)) return("only in b")
    if (is.na(fb)) return("only in a")
    if (identical(fa, fb)) return("identical")
    ok <- tryCatch(formulas_agree(fa, fb, vocab, levels, n = n_sample),
                   error = function(e) NA)
    if (is.na(ok)) "uncomparable" else if (ok) "equivalent" else "differs"
  }, "")

  out[order(factor(out$verdict, levels = c("differs", "only in a", "only in b",
                                           "equivalent", "identical", "uncomparable")),
            out$node, out$value), ]
}
