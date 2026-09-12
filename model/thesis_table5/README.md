# Thesis Table 5

The rule table as printed in the thesis (`Results_PhDthesis.docx`, Table 5):
**76 rule-governed nodes, 89 rules, 5 inputs**.

This is *not* the table that was simulated in 2017. `../rules.txt` is. Table 5
was assembled for the write-up two years later, and the two differ in ways that
change the model's behaviour — see [COMPARISON.md](../../COMPARISON.md).

## Provenance

Ported from `src/rules_data.py` in the
[Python rebuild](https://github.com/ncdomingues/tcell-mirna-model), which
transcribed Table 5 from the thesis document. Converting that list of
`(node, value, formula, description)` tuples into this repository's
tab-separated format is mechanical: node indices are assigned inputs-first then
in order of first appearance, and each node's Max level is the highest value any
of its rules targets.

| file | what |
|---|---|
| `nodes.txt` | 81 nodes — 5 inputs, 76 rule-governed |
| `rules.txt` | 89 rules |
| `annotations.tsv` | each rule with the literature rationale from Table 5 |

`annotations.tsv` is the reason this table is worth keeping alongside the
working one: it carries the published justification for every edge, which the
2017 file never had.

## One thing to know before using it

`STAT5:2 = IL2R:2 | IL4R:2` reads IL2R at level 2, but every IL2R rule in this
table targets level 1, so that branch can never fire. Loading it warns:

```
threshold reference above the node's reachable level (always false):
IL2R:2 but IL2R tops out at 1
```

The slip is in the thesis, not in the transcription — Table 5 gives
`IL2R | 1 | IL2 & CGC & STAT5:2` while its own annotation on the same row says a
high level of IL2R is what promotes a high level of STAT5. The 2017 working
table has the same rule at level 2. It is left as published; the warning is the
documentation.
