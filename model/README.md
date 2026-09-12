# The tables

A model is two tab-separated files. The format has not changed since 2017 and
is described in `original/base_files/method_modelling.txt`.

```
nodes.txt   Node  Name   Max level  Input
            7     ZAP70  2          0

rules.txt   Target  Name   Value  Formulae
            7       ZAP70  2      LCK:2
```

`Target` is a row index into the node table -- **the engine targets by index and
ignores `Name`**, which is how two `IL4R` rules drove FOXO1 for two years. A node
with several rules takes the highest `Value` whose formula is true; a node with
no rule holds its initial value. Formula syntax: `&` AND, `|` OR, `!` NOT,
`NODE:k` "NODE is at level at least k".

## My model, 2017

| file | what |
|---|---|
| `nodes.txt` | 91 nodes, verbatim |
| `rules.txt` | 86 formulas, verbatim -- fails strict validation, on purpose |
| `rules_corrected.txt` | the same, with two `IL4R` rows retargeted 18 -> 66 |

`rules.txt` is kept exactly as it was so the 2017 runs can be reproduced with
`ginsim_run(..., strict = FALSE)`. Everything in `scripts/` uses the corrected
table.

## `mir34c_variants/`

Twelve rule tables from February and March 2017, each proposing a different
formula for the miR-34c-5p promoter (node 84). Named for the transcription
factors they invoke, or for the date they were saved.

They are **not** twelve one-line edits of one network -- they are three backbone
generations that also differ in eight other rules, and none matches
`rules.txt`. `scripts/03_mir34c_tf_variants.R` groups them and only compares
within a group. Three of them carry a `TGFBRIL6R` typo that
`repair_2017_table()` fixes. See MODEL_NOTES.md.

## Base networks

`naldi2010/` and `abou_jaoude2015/` are the published networks this model was
built on, as transcribed in 2017 and already extended there with the two miRNAs.

| file | what |
|---|---|
| `naldi2010/nodes.txt` | 82 nodes, verbatim -- `IL4RA` ceiling is wrong |
| `naldi2010/nodes_corrected.txt` | `IL4RA` raised from 1 to 2, as its rules require |
| `naldi2010/rules.txt` | 63 formulas, verbatim |
| `abou_jaoude2015/nodes.txt` | 118 nodes, verbatim |
| `abou_jaoude2015/rules.txt` | 103 formulas, verbatim -- one stale `Name` label |
| `abou_jaoude2015/rules_corrected.txt` | that label changed `CREB` -> `CREB1` |

The Naldi network carries real exogenous cytokine inputs, which is why the
biological validation in `scripts/04` runs on it rather than on my own table.

Sources: Naldi et al. (2010) *PLoS Comput Biol* 6:e1000912; Abou-Jaoude et al.
(2015) *Front Bioeng Biotechnol* 3:86; proximal signalling after Saez-Rodriguez
et al. (2007) *PLoS Comput Biol* 3:e163.
