# run_all.R --------------------------------------------------------------------
#
# Everything, in order, from a clean checkout:
#
#   Rscript run_all.R
#
# Takes a few minutes on a laptop and rewrites output/ and figures/ (except
# figures/original_2017/, which is the 2017 record and is never regenerated).
# ------------------------------------------------------------------------------

scripts <- c("tests/test_engine.R",
             "scripts/01_stimulations.R",
             "scripts/02_mirna_conditions.R",
             "scripts/03_mir34c_tf_variants.R",
             "scripts/04_validate_naldi.R")

for (s in scripts) {
  cat("\n", strrep("=", 78), "\n", s, "\n", strrep("=", 78), "\n", sep = "")
  started <- Sys.time()
  source(s, echo = FALSE, local = new.env())
  cat(sprintf("-- %s took %.0fs\n", s,
              as.numeric(difftime(Sys.time(), started, units = "secs"))))
}

cat("\nall done\n")
