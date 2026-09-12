# The T-cell model in R

The original GINsim/R logical model from my PhD -- CD4+ T-cell activation and
differentiation, extended with miR-34c-5p and miR-155-5p regulation -- recovered
from the 2017 working directory, repaired, documented, and made to run again.

This is the R side of a pair. [`tcell-mirna-model`](https://github.com/ncdomingues/tcell-mirna-model)
is a from-scratch Python rebuild of the same biology from the published rule
table. **This repository is the actual 2017 code and data**: the simulation
engine, the rule tables, the twelve competing hypotheses for the miR-34c-5p
promoter, the base networks it was built on, and the figures as they came out at
the time. Thesis: *sncRNA regulatory networks in T cell activation and viral
response*, University of Lisbon, 2019.

Everything here runs on a bare R install. No packages, no GINsim, no MaBoSS.

## Running it

```bash
Rscript run_all.R
```

Around six minutes, and it rewrites `output/` and `figures/`. Individually:

```bash
Rscript tests/test_engine.R            # 52 tests, ~1s
Rscript scripts/01_stimulations.R      # five stimulation conditions
Rscript scripts/02_mirna_conditions.R  # the four-condition miRNA comparison
Rscript scripts/03_mir34c_tf_variants.R  # twelve promoter hypotheses
Rscript scripts/04_validate_naldi.R    # does the biology hold up at all?
Rscript scripts/05_compare_tables.R    # working table vs. the thesis's Table 5
```

## What's here

```
R/
  ginsim_engine.R    the asynchronous multi-valued simulator
  model_io.R         loading tables, building states by node name, knock-outs
  analysis.R         summaries and plots
model/
  nodes.txt          91 nodes, verbatim from 2017
  rules.txt          86 formulas, verbatim -- including two mis-targeted rows
  rules_corrected.txt  the same table with those two rows fixed
  mir34c_variants/   twelve competing formulas for the miR-34c-5p promoter
  thesis_table5/     Table 5 as printed in the thesis, with its annotations
  naldi2010/         Naldi et al. 2010, as transcribed and extended in 2017
  abou_jaoude2015/   Abou-Jaoude et al. 2015, likewise
scripts/             the five analyses, in order
comparison/          the Python rebuild's published numbers, for cross-checking
tests/test_engine.R  52 tests: notation, update semantics, table integrity
original/            the 2017 files, untouched -- engine, scripts, procedures
ginsim/              the .zginml network files, openable in GINsim
figures/             regenerated figures, plus figures/original_2017/
output/              CSV results
MODEL_NOTES.md       what the model is, and every repair made to it
VALIDATION.md        what has actually been checked, and what failed
COMPARISON.md        the working table against the thesis's published Table 5
```

## The model

91 nodes, 86 logical formulas. Proximal TCR signalling after Saez-Rodriguez et
al. 2007; T-helper differentiation after Naldi et al. 2010 and Abou-Jaoude et
al. 2015; the two miRNAs and their targets are mine.

Formula syntax is the thesis's own: `&` AND, `|` OR, `!` NOT, and `NODE:k` for
"NODE is at level at least k". A node with several rules takes the highest value
whose formula is true. Updates are asynchronous: each step moves exactly one
node, picked at random among those that disagree with their target, one level
towards it.

Two tab-separated tables define a model, and that format has not changed since
2017:

```
nodes.txt   Node  Name   Max level  Input
            7     ZAP70  2          0

rules.txt   Target  Name   Value  Formulae
            7       ZAP70  2      LCK:2
```

## What it found

**Under TCR + CD28 stimulation the network never settles.** Not in 200 steps,
not in any of 200 trajectories, in any condition. The cause is a rule in the
table itself: `TCR = APC-Antigen & CD3 & CD28 & !IL2:2`. IL2 accumulates,
reaches level 2, shuts TCR down, the whole proximal arm collapses, IL2 decays
and TCR restarts -- a cycle of roughly 60 update steps. That is the model's
IL2 feedback working exactly as written, but it means the 2017 notes' talk of
steady states does not apply to this network. Reported numbers are therefore
occupancy averaged over the second half of the run, not a final state.

**Both miRNAs together shift the cell from Th17 towards iTreg.** Occupancy of
each master regulator, TCR + CD28 stimulation, 200 trajectories:

| Subset | Node | no miRNA | miR-34c-5p | miR-155-5p | both |
|---|---|---|---|---|---|
| Th1 | TBX21 | 0.71 | 0.71 | 0.70 | 0.71 |
| Th2 | GATA3 | 0.00 | 0.00 | 0.00 | 0.00 |
| Th17 | RORC | **0.36** | 0.01 | **0.97** | **0.01** |
| iTreg | FOXP3 | 0.20 | 0.40 | 0.04 | **0.50** |
| Tfh | BCL6 | 0.18 | 0.00 | 0.15 | 0.00 |
| Th22 | AHR | 0.16 | 0.00 | 0.00 | 0.00 |
| Th9 | SPI1-PU1 | 0.19 | 0.42 | 0.11 | 0.48 |

miR-155-5p alone drives RORC almost to saturation. Adding miR-34c-5p collapses
it to nothing and doubles FOXP3 -- miR-34c-5p is epistatic to miR-155-5p on the
Th17 axis. Th1 is untouched by either. That Th17/iTreg shift is what the thesis
argued for, and it survives the repairs made here.

**The promoter hypothesis matters enormously, and the twelve saved tables can't
be compared directly.** They are three backbone generations, not twelve one-line
edits -- eight other rules also differ between them, because the rest of the
network was being revised at the same time. `scripts/03` prints the backbone
group for each and only compares within groups. Within a group the spread is
still large: under backbone 2, `TP53_FOXO3` puts miR-34c-5p at 0.29 while `MYC`
puts it at essentially zero, and RORC moves from 0.38 to 0.89 as a result.

**The base network reproduces textbook T-helper biology; my miR-155-5p
annotations cost it Th1.** On Naldi et al. 2010, which unlike my own table has
real exogenous cytokine inputs, all four polarising cocktails drive the right
master regulator -- 4/4. Adding miR-34c-5p keeps it at 4/4. Adding miR-155-5p
drops it to 3/4: IL12 no longer commits to Th1. The IL12 signal arrives fine --
IL12R and STAT4 are both on -- but in this network TBX21 hangs off STAT1 rather
than STAT4, and the only route to STAT1 is IFNGR, which carries a
`!miR-155-5p` term I added in 2017. The cell still makes IFN-gamma; it just
cannot hear it.

That is a failed positive control, not a curiosity: miR-155-deficient T cells
are Th2-biased and make *less* IFN-gamma, so miR-155-5p should promote Th1, not
abolish it. It is reported rather than tuned away. VALIDATION.md has the
per-node trace and what it means for the numbers above.

## The working table is not the published table

The thesis printed its rule table as Table 5, two years after the simulations
were run. `model/thesis_table5/` holds it, and `scripts/05` runs both tables
through the same engine under the same protocol, so differences are attributable
to the tables rather than to implementations. They had never been compared.

They agree on the destination — **5 of 7 subsets**, including both headline
ones: Th17 down, iTreg up. The thesis's central claim holds either way.

They disagree completely on **which miRNA gets there**. Wherever the two miRNAs
pull in different directions, the working table lands on miR-34c-5p's outcome,
4 times out of 4; Table 5 lands on miR-155-5p's, 5 times out of 5. That is why
the Python rebuild reports miR-155-5p as the dominant driver while this
repository reports miR-34c-5p as epistatic on the Th17 axis. It was never the
code — the same Table 5 run through two independent engines agrees in all seven
directions.

It traces to a single rule, and it is the one rule the thesis was trying to
determine:

```
working table   miR-34c-5p = (GATA3 | FOS | MYC | TP53 | FOXO3 | SP1) & !STAT3
thesis Table 5  miR-34c-5p = GATA3 & MYC
```

A six-way OR against a conjunction. In the working table miR-34c-5p is on 98.5%
of the time; under Table 5, 4.8%. Table 5 actually gives it *more* target edges
(16 against 7) — they just never fire. Full workings in COMPARISON.md, including
a load-bearing typo in Table 5 that kills the IL2 to STAT5 amplification loop.

## Repairs

The 2017 tables and engine carried four defects, all found by validation added
here, all documented in MODEL_NOTES.md and none fixed silently:

1. The engine never substituted `NODE:1` references, so `LCK:1` became `0:1` --
   an R integer sequence, not a truth value. This silently switched off the
   PI3K/AKT1 arm of the network and is the source of the
   `number of items to replace is not a multiple of replacement length`
   warnings recorded in the original notes.
2. Two `IL4R` rules carried target index 18, which is FOXO1. They drove the
   wrong node for two years, and pushed it past its declared ceiling.
3. The Naldi node table gives `IL4RA` a maximum level of 1 while its own rules
   assign it level 2 and other rules read `IL4RA:2`.
4. Three of the promoter variants read `TGFBRIL6R` -- two node names with the
   operator between them missing. R parses that as a function call and only
   fails on evaluation, so it never surfaced.

`model/rules.txt` is still the verbatim table, and `ginsim_run(..., strict =
FALSE)` will still run it as it ran in 2017.
