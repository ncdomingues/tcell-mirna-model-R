# ginsim_engine.R --------------------------------------------------------------
#
# Asynchronous simulator for GINsim-style multi-valued logical models.
#
# This is a documented, hardened refactor of the engine used for my 2019 PhD
# thesis (`original/base_files/GinSimR_ND.R`). The update semantics are
# unchanged -- see MODEL_NOTES.md for the three behaviour-preserving fixes.
#
# A model is two tab-separated tables:
#
#   nodes.txt   Node  Name           Max level  Input
#               7     ZAP70          2          0
#
#   rules.txt   Target  Name   Value  Formulae
#               7       ZAP70  2      LCK:2
#
# Formula syntax: `&` AND, `|` OR, `!` NOT, `NODE:k` "NODE is at level >= k".
# Nodes with no targeting rule are inputs: they hold their initial value.
# A node with several rules takes the highest Value whose formula is true.
# ------------------------------------------------------------------------------


#' Order node indices by decreasing name length
#'
#' Rule formulae are evaluated by textual substitution, so longer names must be
#' replaced first: substituting `IL2` before `IL2R` would corrupt `IL2R`.
get_order_var <- function(nodes) {
  order(nchar(nodes[, 2]), decreasing = TRUE)
}


#' Validate a node/rule table pair, and fail loudly rather than mid-simulation
#'
#' The 2017 scripts ran with `read.delim` defaults and swallowed malformed rows
#' as cryptic `eval(parse(...))` warnings. This front-loads those checks.
check_model <- function(nodes, rules, state = NULL, strict = TRUE) {
  fail <- if (strict) stop else function(...) warning(..., call. = FALSE)
  stopifnot(ncol(nodes) >= 4, ncol(rules) >= 4)

  if (anyNA(nodes[, 1:3]) || any(!nzchar(trimws(nodes[, 2])))) {
    stop("nodes table has missing Node/Name/Max level entries")
  }
  if (anyNA(rules[, c(1, 3)]) || any(!nzchar(trimws(rules[, 4])))) {
    bad <- which(is.na(rules[, 1]) | is.na(rules[, 3]) | !nzchar(trimws(rules[, 4])))
    stop("rules table has empty Target/Value/Formulae on row(s): ",
         paste(bad, collapse = ", "))
  }

  out_of_range <- rules[, 1][rules[, 1] < 1 | rules[, 1] > nrow(nodes)]
  if (length(out_of_range)) {
    stop("rules target node index/indices outside the node table: ",
         paste(unique(out_of_range), collapse = ", "))
  }

  # A rule's Name column is documentation only -- the engine targets by index.
  # Where the two disagree the table is almost certainly wrong (see the IL4R
  # case in MODEL_NOTES.md), so say so instead of silently driving the wrong node.
  declared <- trimws(rules[, 2])
  actual   <- trimws(nodes[rules[, 1], 2])
  mismatch <- which(declared != actual)
  if (length(mismatch)) {
    warning("rule Name does not match the targeted node on row(s) ",
            paste(mismatch, collapse = ", "), ": ",
            paste0(declared[mismatch], " -> node ", rules[mismatch, 1],
                   " (", actual[mismatch], ")", collapse = "; "),
            call. = FALSE)
  }

  over <- which(rules[, 3] > nodes[rules[, 1], 3])
  if (length(over)) {
    fail("rule target value exceeds the targeted node's Max level on row(s): ",
         paste(over, collapse = ", "))
  }

  if (!is.null(state)) {
    if (length(state) != nrow(nodes)) {
      stop("state vector has ", length(state), " entries but the model has ",
           nrow(nodes), " nodes")
    }
    if (any(state < 0) || any(state > nodes[, 3])) {
      stop("state vector has values outside the nodes' declared levels")
    }
  }

  invisible(TRUE)
}


#' Compile a rule table into evaluable R expressions
#'
#' The 2017 engine rebuilt every formula from text on every single update --
#' 91 nodes x 86 formulas of `gsub` per step, which is why a 200-step, 200-run
#' ensemble took the best part of an afternoon. Names resolve to node indices
#' once, here, and the simulation loop then just evaluates the expressions
#' against the state vector. The semantics are identical; only the cost changes.
#'
#' Substitution still runs longest-name-first, and writes through an
#' intermediate marker (`<index#level>`) so that the R code being
#' inserted can never itself be matched by a later, shorter node name.
#'
#' @return a list of expressions, one per rule row, each referring to `S`
compile_rules <- function(nodes, rules) {
  order_var <- get_order_var(nodes)
  text      <- trimws(rules[, 4])

  # A formula may reference a level the node table does not declare. In the
  # thesis Table 5 transcription, `STAT5 = IL2R:2 | IL4R:2` reads IL2R at level
  # 2 while every IL2R rule targets level 1, so that term can never be true --
  # a dead branch, and worth saying so out loud. Substitute up to whatever level
  # is actually referenced so the formula still compiles and evaluates false,
  # which is what the reference means.
  tokens <- unlist(regmatches(text, gregexpr("[A-Za-z][A-Za-z0-9_.-]*:[0-9]+", text)))
  referenced <- rep(0L, nrow(nodes))
  if (length(tokens)) {
    tok_name  <- sub(":[0-9]+$", "", tokens)
    tok_level <- as.integer(sub("^.*:", "", tokens))
    for (i in seq_along(tokens)) {
      idx <- match(tok_name[i], nodes[, 2])
      if (!is.na(idx)) referenced[idx] <- max(referenced[idx], tok_level[i])
    }
  }
  dead <- which(referenced > nodes[, 3])
  if (length(dead)) {
    warning("threshold reference above the node's reachable level (always ",
            "false): ",
            paste0(nodes[dead, 2], ":", referenced[dead], " but ",
                   nodes[dead, 2], " tops out at ", nodes[dead, 3],
                   collapse = "; "), call. = FALSE)
  }

  for (idx in order_var) {
    var_name  <- nodes[idx, 2]
    max_value <- max(nodes[idx, 3], referenced[idx])

    # `NAME:k` means "NAME is at level >= k". Resolve the thresholds before the
    # bare name, highest first, so `IL2R:2` is never clipped to `IL2R` + ":2".
    #
    # The 2017 engine stopped this loop at k = 2, so `LCK:1` survived
    # substitution as `0:1` -- R's integer sequence c(0, 1), not a truth value.
    # That is the source of the "number of items to replace is not a multiple of
    # replacement length" warnings in the original notes, and it silently
    # switched off the whole PI3K/AKT1 arm. See MODEL_NOTES.md.
    for (k in seq(max_value, 1)) {
      text <- gsub(paste0(var_name, ":", k),
                   paste0("<", idx, "#", k, ">"), text, fixed = TRUE)
    }
    text <- gsub(var_name, paste0("<", idx, "#1>"), text, fixed = TRUE)
  }

  # `TRUE`/`FALSE` are legitimate formulae -- the Abou-Jaoude table uses them
  # for constitutively expressed receptor chains. Anything else still carrying
  # letters is a node name the table never declared.
  residue  <- gsub("<[0-9]+#[0-9]+>|TRUE|FALSE", "", text)
  leftover <- grepl("[A-Za-z]", residue)
  if (any(leftover)) {
    stop("unresolved name(s) in formula/formulae: ",
         paste(unique(rules[leftover, 4]), collapse = "; "))
  }

  text     <- gsub("<([0-9]+)#([0-9]+)>", "(S[\\1] >= \\2)", text)
  compiled <- lapply(text, function(f) parse(text = f)[[1]])

  # Parsing is not enough. A missing operator between two node names -- the
  # `TGFBRIL6R` typo in three of the 2017 variant tables -- resolves to
  # `(S[61] >= 1)(S[60] >= 1)`, which R parses happily as a function call and
  # only rejects on evaluation. Probe every expression against an all-off and an
  # all-on state, so a broken formula surfaces here rather than mid-simulation.
  probe <- function(S) {
    for (i in seq_along(compiled)) {
      val <- tryCatch(eval(compiled[[i]]), error = function(e) e)
      if (inherits(val, "error")) {
        stop("formula does not evaluate: '", rules[i, 4], "' (",
             conditionMessage(val), ")")
      }
      if (length(val) != 1L || is.na(val) || !is.logical(val)) {
        stop("formula does not reduce to a single truth value: '",
             rules[i, 4], "'")
      }
    }
  }
  probe(integer(nrow(nodes)))
  probe(as.integer(nodes[, 3]))

  compiled
}


#' One asynchronous update step
#'
#' Computes the synchronous target state, then moves exactly one node -- picked
#' uniformly at random among those that disagree with their target -- a single
#' level towards it. This is GINsim's asynchronous, non-deterministic scheme.
#'
#' @param compiled  output of `compile_rules()`
#' @param targets   list mapping each node index to the rows targeting it
#'
#' @return the updated state vector (unchanged if the state is a steady state)
ginsim_next <- function(state, rules, compiled, targets) {
  nvar <- length(state)
  env  <- list2env(list(S = state), parent = baseenv())

  solved <- vapply(compiled, eval, logical(1), envir = env, USE.NAMES = FALSE)

  # Target value per node: the highest Value among its satisfied rules, 0 if
  # none is satisfied, and the current value for input nodes (no rules at all).
  sync_target <- integer(nvar)
  for (i in seq_len(nvar)) {
    rows <- targets[[i]]
    if (length(rows) == 0L) {
      sync_target[i] <- state[i]
    } else {
      hit <- solved[rows]
      sync_target[i] <- if (any(hit)) max(rules[rows[hit], 3]) else 0L
    }
  }

  changed <- which(state != sync_target)
  if (length(changed) > 0) {
    # `sample()` on a length-1 vector samples 1:n, so pick explicitly.
    to_change <- if (length(changed) == 1L) changed else sample(changed, 1)
    state[to_change] <- state[to_change] +
      if (sync_target[to_change] > state[to_change]) 1L else -1L
  }
  state
}


#' Run an ensemble of asynchronous simulations
#'
#' @param state   initial state vector, one entry per node
#' @param nodes   node table (Node, Name, Max level, Input)
#' @param rules   rule table (Target, Name, Value, Formulae)
#' @param ntime   time points per trajectory
#' @param nrep    independent trajectories to average
#' @param seed    optional RNG seed, for reproducible output
#' @param strict  stop on a malformed rule table. Set FALSE to reproduce the
#'                2017 runs on the verbatim `model/rules.txt`, which contains
#'                two mis-targeted rows (see MODEL_NOTES.md); they then warn
#'                and run exactly as the original engine ran them.
#'
#' @return a list with
#'   `evol`    nodes x (ntime+1) matrix, mean level per node per time point
#'   `ss`      matrix of steady states reached, one per row (see note below)
#'   `n_steady` how many of the `nrep` trajectories reached a steady state
#'   `nodes`, `ntime`, `nrep`, `seed` -- the run's provenance
#'
#' Note: `ss` may contain duplicate rows -- one per trajectory that settled.
#' Unlike the 2017 version it does not prepend the initial state, so every row
#' really is a steady state. Use `unique(run$ss)` for the distinct attractors.
ginsim_run <- function(state, nodes, rules, ntime, nrep, seed = NULL,
                       strict = TRUE) {
  check_model(nodes, rules, state, strict = strict)
  if (!is.null(seed)) set.seed(seed)

  compiled   <- compile_rules(nodes, rules)
  targets    <- lapply(seq_len(nrow(nodes)), function(i) which(rules[, 1] == i))
  nvar       <- length(state)
  state_cum  <- matrix(0, nrow = nvar, ncol = ntime + 1)
  steady     <- matrix(numeric(0), nrow = 0, ncol = nvar)

  for (rep in seq_len(nrep)) {
    traj      <- matrix(0L, nrow = nvar, ncol = ntime + 1)
    traj[, 1] <- state
    settled   <- FALSE

    for (t in seq_len(ntime)) {
      if (settled) {
        traj[, t + 1] <- traj[, t]
        next
      }
      traj[, t + 1] <- ginsim_next(traj[, t], rules, compiled, targets)
      if (all(traj[, t + 1] == traj[, t])) {
        settled <- TRUE
        steady  <- rbind(steady, traj[, t])
      }
    }
    state_cum <- state_cum + traj
  }

  state_cum <- state_cum / nrep
  rownames(state_cum) <- nodes[, 2]
  if (nrow(steady) > 0) colnames(steady) <- nodes[, 2]

  list(evol = state_cum, ss = steady, n_steady = nrow(steady),
       nodes = nodes, ntime = ntime, nrep = nrep, seed = seed)
}
