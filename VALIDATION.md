# Validation

Three questions, answered separately, because they fail independently:

1. **Does the engine implement the notation it claims to?** -- `tests/test_engine.R`
2. **Are the digitized tables structurally sound?** -- also `tests/test_engine.R`
3. **Does the network reproduce real T-cell biology?** -- `scripts/04_validate_naldi.R`

Run all three yourself:

```bash
Rscript tests/test_engine.R
Rscript scripts/04_validate_naldi.R
```

## 1 & 2. Code and tables -- 52/52 tests pass

`tests/test_engine.R` needs no packages and takes about a second.

**Name resolution.** Bare names are true at level >= 1. `NODE:k` is true at
level >= k. Longer names are substituted before names they contain, so `IL2`
cannot corrupt `IL2R`. Hyphenated names (`miR-34c-5p`, `CDKN1A-p21`) resolve.
`TRUE` is accepted as a constitutive formula. An undeclared name is an error,
not a silent false.

Two of these tests exist because the 2017 engine got them wrong:

- `NODE:1 means the same as the bare name` locks in the fix for the `0:1`
  sequence bug that switched off the PI3K arm (MODEL_NOTES.md, repair 1).
- `a missing operator between two names is caught, not parsed as a call` locks
  in the compile-time evaluation probe that catches `TGFBRIL6R` (repair 4). R
  parses that as a function call, so a syntax check alone passes it.

**Operator semantics.** AND, OR, NOT, and the precedence that follows from
evaluating in R: NOT before AND before OR, with parentheses overriding.

**Update semantics.** A node takes the highest value among its satisfied rules
and falls to 0 when none is satisfied. Nodes with no rule hold their value.
Exactly one node changes per step and only by one level, so 0 -> 2 never happens
in a single update. A steady state is recorded when reached. The same seed
reproduces the same ensemble.

**Table validation.** Out-of-range targets, values above a node's ceiling, and
wrong-length state vectors are rejected; a rule whose `Name` disagrees with the
node it targets warns. The same table passes under `strict = FALSE`, which is
what lets the verbatim 2017 table still run as it ran.

**The tables themselves.** 91 nodes and 86 rules. Both miRNAs rule-governed
rather than inputs. All seven master transcription factors present. The verbatim
table is *rejected* under strict validation and the corrected one passes, and
the correction is asserted to move exactly two rows and change nothing else.
Twelve promoter variants, each proposing a distinct formula, all compiling once
repaired, falling into exactly three backbone generations, none of which matches
the published table. Both base networks compile. `repair_2017_table()` is
asserted to fix both known slips, touch nothing else, and be a no-op on a clean
table.

## 3. Biology

### The test

Naldi et al. 2010 -- one of the two published networks this model was built on
-- was validated by its authors by showing that each polarising cytokine
cocktail drives its matching master transcription factor. That test is run here
on the 2017 transcription of that network, in four versions: base, plus each
miRNA separately, plus both. 200 trajectories, 200 steps, occupancy over the
second half.

This network is used for the test rather than my own because it has real
exogenous cytokine inputs (`IL4_e`, `IL12_e`, `TGFB_e`, `IL6_e`). My own table
does not, which is a limitation of my table, discussed below.

A condition passes when the expected regulator rises by more than 0.05 over the
TCR + IL2 baseline *and* is the largest riser among the four.

### Results

| Network | + IL12 -> TBX21 | + IL4 -> GATA3 | + TGFB + IL6 -> RORC | + TGFB -> FOXP3 | |
|---|---|---|---|---|---|
| base (no miRNA) | pass (1.00) | pass (1.00) | pass (1.00) | pass (1.00) | **4/4** |
| + miR-34c-5p | pass (1.00) | pass (1.00) | pass (1.00) | pass (1.00) | **4/4** |
| + miR-155-5p | **fail** (0.00) | pass (1.00) | pass (0.49) | pass (0.26) | **3/4** |
| + both | **fail** (0.00) | pass (0.96) | pass (0.49) | pass (0.26) | **3/4** |

**The base network is 4/4.** The 2017 transcription of Naldi et al. 2010 --
once `IL4RA`'s ceiling is corrected and the constitutive receptor chains are
held on -- reproduces every textbook polarisation. That is the strongest single
piece of evidence that the digitization is sound, and it is independent of
anything to do with miRNAs.

**miR-34c-5p costs nothing.** Still 4/4.

**miR-155-5p costs Th1, and the cause is one edge.** `scripts/04` prints the
chain (`output/naldi_th1_chain.csv`):

| node | base | + miR-34c-5p | + miR-155-5p | + both |
|---|---|---|---|---|
| miR-155-5p | 0 | 0 | 1 | 1.000 |
| IL12R | 1 | 1 | 1 | 0.835 |
| STAT4 | 1 | 1 | 1 | 0.835 |
| IFNG | 1 | 1 | 1 | 0.835 |
| **IFNGR** | 1 | 1 | **0** | **0.000** |
| **STAT1** | 1 | 1 | **0** | **0.000** |
| **TBX21** | 1 | 1 | **0** | **0.000** |

The IL12 signal arrives intact -- IL12R and STAT4 are both on. But in this
network `TBX21 = (TBX21 | STAT1) & !GATA3` hangs off STAT1, not STAT4, and with
IFNB_e and IL27_e off the only route to STAT1 is
`IFNGR = IFNGR1 & IFNGR2 & (IFNG | IFNG_e) & !miR-155-5p`. The cell still makes
IFNG; it just cannot hear it. One `!miR-155-5p` term closes the entire Th1 axis.

That is a mechanistically clean result, and it is also a **problem for the
model**, because it contradicts the literature.

### The miR-155-5p positive control, and where it fails

miR-155-5p is the immune system's best-characterised miRNA, which makes it the
natural control: whatever the model says about the poorly-characterised
miR-34c-5p is only worth reading if the model gets miR-155-5p right.

miR-155-deficient T cells are Th2-biased and produce less IFN-gamma, and
miR-155 is required for Th17 responses (Rodriguez et al. 2007; Thai et al. 2007;
O'Connell et al. 2010). So miR-155-5p should be **Th1-promoting** and
**Th17-promoting**.

| Expectation | My table (`scripts/02`) | Naldi extension (`scripts/04`) |
|---|---|---|
| Th17-promoting | **matches** -- RORC 0.36 -> 0.97 | **contradicts** -- RORC 1.00 -> 0.49 |
| Th1-promoting | uninformative -- TBX21 0.71 -> 0.70 | **contradicts** -- TBX21 1.00 -> 0.00 |

Split, and not in the model's favour. My own table gets the Th17 direction
firmly right and says nothing either way about Th1. The Naldi extension gets
both backwards, and the Th1 reversal traces to the single `!miR-155-5p` term on
IFNGR that I added in 2017 -- an edge that is defensible in isolation (miR-155
does target IFNGR-pathway components) but which, in a network where TBX21 reads
only STAT1, becomes a hard off-switch for Th1 rather than a modulation.

I am not tuning that edge away. It is what the 2017 table says, and the
disagreement is more useful recorded than smoothed over: it says that the
miR-155-5p annotations in the Naldi extension need revisiting before anything
that network says about Th1 is trusted.

### Known gaps

**GATA3 never engages in my own table.** Under every condition tested, in every
miRNA combination, GATA3 sits at 0 -- so the Th2 row of the main result table is
uninformative rather than a finding, and the Th1 half of the miR-155-5p control
cannot be run on that table at all. The cause is structural and is spelled out
in MODEL_NOTES.md: GATA3 needs STAT6, STAT6 needs IL4R, IL4R needs IL4, and IL4
needs STAT6 or GATA3, with no exogenous IL4 input to break into the loop. The
variant backbones added a STAT5 route to GATA3, which fixes it; the published
table did not.

**No steady states.** The attractor counts in
`output/stimulation_attractors.csv` are zero across all five conditions and 200
trajectories. This is a cyclic attractor, not a convergence failure, and its
cause is identified in MODEL_NOTES.md. It does mean that none of the
steady-state analysis the 2017 notes assumed is available for this network.

**The Th17 and iTreg stimulation conditions in `scripts/01` are mine.** Only
IL2, Th1 and Th2 were scripted in 2017.

## Bottom line

The engine does what it claims (52/52). The digitized Naldi network reproduces
all four textbook polarisations, which validates the transcription. miR-34c-5p's
annotations cost nothing on that test. miR-155-5p's annotations do: they abolish
Th1 commitment through a single edge, in the opposite direction to the
literature, and they reverse the Th17 direction that my own table gets right.

So the Th17 -> iTreg shift reported in the README rests on a network whose
miR-155-5p wiring fails its own positive control on a related network. That is
worth knowing before quoting the number.
