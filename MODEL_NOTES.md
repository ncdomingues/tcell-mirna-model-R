# Modelling notes

What the model is, what was changed to make it run, and what was left alone.

The rule is simple: **nothing in the 2017 material was edited in place.**
`model/nodes.txt` and `model/rules.txt` are the tables as they were, and
everything under `original/` is byte-for-byte what was in the thesis directory.
Corrections live in separate files, or in `repair_2017_table()`, which reports
every change it makes.

## The network

91 nodes, 86 logical formulas, built in three layers:

- proximal TCR signalling from Saez-Rodriguez et al. 2007 (LCK, ZAP70, PI3K,
  the NFAT and NFKB arms);
- T-helper differentiation from Naldi et al. 2010 and Abou-Jaoude et al. 2015
  (the STAT/cytokine receptor layer and the master transcription factors);
- miR-34c-5p and miR-155-5p, their transcriptional inputs and their targets,
  which is the part that is mine.

Five other miRNAs (miR-181a-5p, miR-21-5p, miR-146b-5p, let-7a-5p, miR-101-3p)
are present as constitutive inputs and are never regulated by the network.

### Semantics

`&` AND, `|` OR, `!` NOT, `NODE:k` "NODE is at level at least k". Precedence is
R's: NOT before AND before OR. This matters -- `A | B & C` is `A | (B & C)` --
and the tables were written in R's notation and evaluated by `eval(parse(...))`,
so R's precedence is not an interpretation, it is what the formulas meant.

A node with several rules takes the highest `Value` whose formula is true, and 0
if none is. A node that no rule targets keeps its initial value forever.

Updates are asynchronous and non-deterministic: compute every node's target,
then move exactly one node -- chosen uniformly at random among those that
disagree with their target -- a single level towards it. One node per step, never
a jump from 0 to 2.

### Inputs versus constitutive components

Both are "nodes with no rule", and the `Input` column does not separate them
consistently across the three tables. `original/base_files/method_modelling.txt`
gives the convention: unregulated nodes with no incoming edges have basal value
1. In the Naldi table those carry `Input = 0` and `constitutive_nodes()` derives
them. In my own table everything unregulated carries `Input = 1`, so the split
is made explicitly in `BASE_STIMULUS` -- APC-Antigen, TCR, CD3, CD28 and the
CD45 isoforms as the stimulus; FOXP1, CREBBP, KAT2B, CGC and the five
constitutive miRNAs as always-present components.

Getting this wrong is silent and total: with the co-factors off, the Naldi
network produces zero activity in every condition, which is how the first run of
`scripts/04` came out.

## The four repairs

### 1. `NODE:1` was never substituted (engine)

The 2017 engine resolved threshold references by looping down to `k = 2`, then
substituted the bare name. So `LCK:1` was left with its `:1` intact and became
`0:1` or `1:1`. In R those are integer sequences: `c(0, 1)` and `1`. Assigning a
length-2 value into `solvedrules[i]` produced

```
number of items to replace is not a multiple of replacement length
```

-- the warning recorded verbatim in `original/base_files/method_modelling.txt`,
and R kept the first element, `0`. Every `NODE:1` reference therefore read as
false whenever the node sat at level 1.

Three rules were affected: `ZAP70 = LCK:1`, `PI3K = LCK:2 | LCK:1`, and
`STAT5 = !IL2R:2 & !IL4R:2 & (IL4R:1 | IL2R:1)`. Since PI3K is the gateway to
PIP3, AKT1, GSK3B and CDKN1A, the entire PI3K arm of the network was dead in the
2017 runs. With `ZAP70 = LCK:1` also broken, ZAP70 could only reach level 2 via
`LCK:2`, which needs CD45RO.

`compile_rules()` now resolves `NODE:k` for every k down to 1. The test
`NODE:1 means the same as the bare name` locks it in.

This is the one repair that changes the model's behaviour substantially, and it
is a fix to the simulator, not to the biology: `NODE:1` and the bare name have
always meant the same thing in GINsim notation.

### 2. Two `IL4R` rules target node 18 (`model/rules.txt`)

```
18	IL4R	2	CGC & IL4 & STAT5
18	IL4R	1	CGC & IL4
66	IL4R	1	IL4
```

Node 18 is FOXO1; IL4R is node 66. The engine targets by index and ignores the
`Name` column, so these two rules drove FOXO1 -- on top of FOXO1's own rule --
and the level-2 row pushed it to 2 against a declared ceiling of 1.

`model/rules_corrected.txt` retargets those two rows to 66 and changes nothing
else. The third row, `66 IL4R 1 IL4`, is redundant with the corrected
`66 IL4R 1 CGC & IL4` given that CGC is constitutive, but it is in the table so
it stays.

Under TCR + CD28 stimulation alone this changes nothing measurable, because IL4
never turns on and neither placement of the rule ever fires. It shows up in the
Th2 condition, where seeding IL4 does engage the receptor.

`check_model()` rejects the verbatim table under `strict = TRUE` and warns under
`strict = FALSE`, which is how `ginsim_run(..., strict = FALSE)` still
reproduces the 2017 behaviour.

### 3. `IL4RA`'s ceiling in the Naldi node table

`model/naldi2010/nodes.txt` declares `37 IL4RA 1`, but the rule table assigns
`37 IL4RA 2 STAT5:2` and two other rules read `IL4RA:2`. Under the 2017 engine a
node declared Boolean skipped the threshold branch entirely, so `IL4RA:2` was
substituted as `1:2` or `0:2` -- another integer sequence, another silent false.

`model/naldi2010/nodes_corrected.txt` sets the ceiling to 2.

### 4. `TGFBRIL6R` (three promoter variants)

`GATA3_MYC_TP53_FOXO3`, `GATA3_TP53_FOXO3` and `MYC_TP53_FOXO3` all contain

```
58	STAT3	1	(TGFBRIL6R & IL23R) & (IL10R | IL1R) & !miR-34c-5p
```

Two node names with the operator between them missing. After substitution this
becomes `(S[61] >= 1)(S[60] >= 1)`, which R parses happily as a function call
and rejects only on evaluation -- so a purely syntactic check passes it.

Every surrounding term is an OR'd receptor, so the missing character is `|`.
`repair_2017_table()` inserts that one character and nothing else.
`compile_rules()` now evaluates every expression against an all-off and an
all-on state at compile time, so this class of error cannot slip through again.

## The cycle

Under TCR + CD28 stimulation none of 200 trajectories reaches a steady state in
200 steps, in any of the five conditions. This is not a failure to converge; it
is a cyclic attractor, and it is written into the table:

```
2	TCR	1	APC-Antigen & CD3 & CD28 & !IL2:2
```

IL2 is autocrine here -- NFAT and NFKB drive it, STAT5 takes it to level 2 --
and at level 2 it shuts TCR off. LCK, ZAP70, the NFAT arm and IL2 itself then
collapse, TCR comes back, and the whole thing repeats with a period of roughly
60 update steps. `figures/trajectories_IL2.png` shows TBX21 and IL2 oscillating
in antiphase.

Biologically this is IL2-mediated downmodulation of TCR signalling, and it is a
deliberate feature of the rule. Practically it means:

- there is nothing useful in `run$ss` for this model, and the attractor counts
  in `output/stimulation_attractors.csv` are all zero;
- a single final state is a snapshot at an arbitrary phase, so every reported
  number is `mean_activity()`: the ensemble mean averaged over the second half
  of the run, which covers the settled behaviour across more than one cycle.

The miRNAs, TP53, MDM2, SIRT1 and BCMcomplex sit outside the cycle -- they latch
and stay latched. The oscillating set is the proximal signalling arm plus IL2,
STAT5, TBX21, FOXP3 and SPI1-PU1.

## Known limitations, not repaired

**GATA3 and the Th2 arm never engage** in my own table under any condition
tested. GATA3 needs STAT6, STAT6 needs IL4R, IL4R needs IL4, and IL4 needs
STAT6 or GATA3 -- a closed loop with no exogenous entry point, since this table
has no `IL4_e` input node. Seeding IL4 in the Th2 condition, as the 2017 script
did, only sets a starting value that the rules immediately overwrite. The
variant backbones fix this by adding a STAT5 route to GATA3
(`(GATA3 | STAT6 | STAT5) & NFKB1 & !miR-34c-5p`), which is why GATA3 is
non-zero throughout `scripts/03` and zero for the published table.

**Th17 and iTreg conditions are mine, not the thesis's.** Only IL2, Th1 and Th2
were scripted in 2017. This table has no exogenous TGFB or IL6 input, so the two
added conditions drive the STAT3 axis through IL1 and IL23, which do exist as
inputs. They are labelled as additions in `POLARISING_INPUTS`.

**The twelve promoter variants are three networks, not one.** Besides the
miR-34c-5p rule they differ in eight others: IKBKG, IL4, STAT3, TGFBR, GATA3,
BCL6, SPI1-PU1 and miR-155-5p. None of the three backbones matches the published
February table. A cross-group comparison measures the backbone revision as much
as the promoter hypothesis, so `scripts/03` reports the group and the test suite
asserts there are exactly three.

**`Value` ordering is not assumed.** The engine takes the maximum over satisfied
rules rather than reading rows top-down, so a table that lists its level-1 rule
before its level-2 rule behaves identically. Several of the 2017 tables do.
