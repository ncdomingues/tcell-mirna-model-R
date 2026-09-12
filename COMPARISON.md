# Working table vs. thesis Table 5

Two rule tables for the same biology, run through the same engine under the
same protocol. Any difference in the results is therefore a difference between
the tables, not between implementations.

- **A. the 2017 working table** — `model/rules_corrected.txt`, 91 nodes, 86
  rules. The file that was actually simulated at the bench.
- **B. thesis Table 5** — `model/thesis_table5/`, 81 nodes, 89 rules. What went
  into the write-up two years later, ported from the digitization in the
  [Python rebuild](https://github.com/ncdomingues/tcell-mirna-model).

Reproduce with `Rscript scripts/05_compare_tables.R`. Protocol for both: 200
trajectories, 200 update steps, seed 2019, occupancy over the second half.
Stimulus derived identically — every node no rule targets is switched on,
except exogenous polarising cytokines. That is TCR + CD28 with autocrine IL2.

They are **not the same network**, and until now they had never been run
against each other.

## Side by side

Master transcription factors, occupancy 0–1:

| | | no miRNA | miR-34c-5p | miR-155-5p | both |
|---|---|---|---|---|---|
| **Th1** TBX21 | working | 0.700 | 0.713 | 0.702 | 0.716 |
| | Table 5 | 0.216 | 0.249 | 0.241 | 0.252 |
| **Th2** GATA3 | working | 0.000 | 0.000 | 0.000 | 0.000 |
| | Table 5 | 0.559 | 0.511 | 0.339 | 0.331 |
| **Th17** RORC | working | 0.348 | 0.010 | **0.965** | 0.015 |
| | Table 5 | 0.271 | 0.226 | 0.202 | 0.168 |
| **iTreg** FOXP3 | working | 0.217 | 0.369 | 0.043 | **0.489** |
| | Table 5 | 0.354 | 0.301 | 0.494 | 0.475 |
| **Tfh** BCL6 | working | 0.202 | 0.002 | 0.158 | 0.002 |
| | Table 5 | 0.403 | 0.345 | 0.147 | 0.169 |
| **Th22** AHR | working | 0.160 | 0.002 | 0.000 | 0.000 |
| | Table 5 | 0.478 | 0.384 | 0.275 | 0.234 |
| **Th9** SPI1-PU1 | working | 0.209 | 0.405 | 0.118 | **0.465** |
| | Table 5 | 0.349 | 0.295 | 0.110 | 0.124 |

`figures/table_comparison.png` is the same data as paired bar charts. The
working table responds like a set of switches; Table 5 responds in gradients.

## They agree on the destination

Direction of the combined-miRNA effect, relative to no miRNA:

| Subset | working | Table 5 | agree |
|---|---|---|---|
| Th1 | flat (+0.02) | flat (+0.04) | yes |
| Th2 | flat (0.00, dead) | down (−0.23) | no |
| **Th17** | **down (−0.33)** | **down (−0.10)** | **yes** |
| **iTreg** | **up (+0.27)** | **up (+0.12)** | **yes** |
| Tfh | down (−0.20) | down (−0.23) | yes |
| Th22 | down (−0.16) | down (−0.24) | yes |
| Th9 | up (+0.26) | down (−0.23) | no |

**5 of 7 agree**, and the two headline subsets are among them: both tables put
the pair of miRNAs on a Th17-down, iTreg-up trajectory. The thesis's central
claim does not depend on which table you use.

Of the two disagreements, Th2 is not really one — GATA3 is structurally dead in
the working table (no exogenous IL4 to break into the STAT6/IL4R/IL4 loop, see
MODEL_NOTES.md), so that row is uninformative rather than contradictory. Th9 is
a genuine contradiction: SPI1-PU1 goes up in one table and down in the other.

## They disagree completely on which miRNA is in charge

This is the real finding. With both miRNAs present, which solo outcome does the
pair land on?

| Subset | working table | thesis Table 5 |
|---|---|---|
| Th1 | solo effects alike | solo effects alike |
| Th2 | solo effects alike | miR-155-5p |
| Th17 | **miR-34c-5p** | solo effects alike |
| iTreg | **miR-34c-5p** | miR-155-5p |
| Tfh | **miR-34c-5p** | miR-155-5p |
| Th22 | solo effects alike | miR-155-5p |
| Th9 | **miR-34c-5p** | miR-155-5p |

Wherever the two miRNAs pull differently, the working table lands on
miR-34c-5p's outcome — every time, 4 for 4. Table 5 lands on miR-155-5p's —
every time, 5 for 5. No exceptions in either direction.

That is why the Python rebuild's headline ("miR-155-5p is the dominant driver,
miR-34c-5p's individual effect is close to baseline") and this repository's
("miR-34c-5p is epistatic to miR-155-5p on the Th17 axis") disagree. It is not
a methodological artefact of Python versus R, or of the two simulation engines.
It is the tables.

## It comes down to one rule

The obvious explanation would be that Table 5 gave miR-34c-5p fewer targets.
The opposite is true — it has **more than twice as many**:

| | miR-34c-5p terms | miR-155-5p terms |
|---|---|---|
| 2017 working table | 7 | 11 |
| thesis Table 5 | 16 | 12 |

So the thing to measure is whether the miRNA is ever expressed at all:

| | working table | Table 5 |
|---|---|---|
| miR-34c-5p occupancy | **0.985** | **0.048** |
| miR-155-5p occupancy | 1.000 | 0.247 |
| STAT3 | 0.015 | 0.254 |
| RORC | 0.015 | 0.168 |

In Table 5, miR-34c-5p is essentially never on — so its sixteen target edges
never fire, and miR-155-5p is left to do the gating alone. The cause is the
promoter rule, the one thing the thesis was actually trying to determine:

```
working table   miR-34c-5p = (GATA3 | FOS | MYC | TP53 | FOXO3 | SP1) & !STAT3
thesis Table 5  miR-34c-5p:2 = GATA3 & MYC & (TP53 | FOXO3 | SP1)
                miR-34c-5p:1 = GATA3 & MYC
```

A six-way OR against a conjunction. TP53 sits at 1.0 in the working table, so
the OR is satisfied continuously and miR-34c-5p is on whenever STAT3 is off.
Table 5 needs GATA3 **and** MYC together, and those sit at 0.34 and 0.25.

The working table then closes a toggle that Table 5 does not:
`STAT3 = (…) & !miR-34c-5p` and `miR-34c-5p = (…) & !STAT3` are mutual
repressors. miR-34c-5p wins, STAT3 is pinned at 0.015, and RORC — which needs
STAT3 — collapses with it. That single loop is the entire Th17 result.

The twelve files in `model/mir34c_variants/` are hypotheses about exactly this
rule. Table 5 is closest to the conjunctive March 2017 line (`10_3_17`,
`1_3_17`, `6_3_17`); the working table is the broad February OR. The write-up
adopted the later, more restrictive hypothesis, and that choice — not any
change in the miRNA's annotated targets — is what flipped which miRNA runs the
network.

## A typo in Table 5

The two tables also differ in ways that are not deliberate. `scripts/05` diffs
them rule by rule, comparing formulas by **behaviour** rather than text — both
are evaluated over 5,000 random states and called equivalent only if they never
disagree, which stops a rename or a regrouping being reported as a rewrite.

```
identical    40      differs        38
equivalent    3      only in A       2      only in B      5
```

One of the two rules present only in the working table matters:

```
working table   IL2R:2 = IL2 & CGC & STAT5:2
thesis Table 5  IL2R:1 = IL2 & CGC & STAT5:2
```

Table 5 assigns that rule to level 1 instead of 2 — where it is redundant, since
`IL2R:1 = IL2 & CGC` already covers it. But `STAT5:2 = IL2R:2 | IL4R:2` reads
IL2R at level 2, which nothing in Table 5 can now produce. The IL2 → STAT5
high-level amplification loop is dead, and the engine says so:

```
Warning: threshold reference above the node's reachable level (always false):
IL2R:2 but IL2R tops out at 1
```

This is in the thesis itself, not in the transcription — `Results_PhDthesis.docx`
row 202 reads `IL2R | 1 | IL2 & CGC & STAT5:2`, while its own annotation on
that row says "A high level of IL2R promotes a high level of STAT5, required to
activate cell proliferation". The intent was level 2. The working table has it
right. It is a typing slip in the published table, and it is load-bearing.

## The Python rebuild's own numbers

For reference, `comparison/python_rebuild_phenotypes.csv` holds what the Python
rebuild published — same Table 5, but a different engine and protocol (400
runs, 10 sweeps, activation frequency at the final state, percentages):

| | no miRNA | miR-34c-5p | miR-155-5p | both |
|---|---|---|---|---|
| Th1 | 27.0 | 30.0 | 36.2 | 38.2 |
| Th2 | 44.0 | 45.2 | 15.5 | 13.8 |
| Th17 | 34.8 | 29.8 | 24.8 | 20.8 |
| iTreg | 23.8 | 22.8 | 31.5 | 30.2 |
| Tfh | 50.5 | 43.5 | 22.2 | 19.2 |
| Th22 | 40.2 | 30.5 | 14.0 | 9.8 |
| Th9 | 37.8 | 33.2 | 4.5 | 3.0 |

Do not read these against the R numbers as if the scales matched — they do not.
Read the directions. `scripts/05` does that check itself and writes it to
`output/table5_vs_python_rebuild.csv`: all seven agree in sign with this
repository's Table 5 run, and miR-155-5p is the dominant driver in both, 7
subsets out of 7. Two independent
engines, two independent protocols, one table, same answer: the divergence
between the repositories was never the code.
