# Results from the Python rebuild

Vendored, not generated. These are the numbers published by
[ncdomingues/tcell-mirna-model](https://github.com/ncdomingues/tcell-mirna-model),
copied here so `scripts/05_compare_tables.R` can cross-check against them
without needing that repository checked out.

| file | what |
|---|---|
| `python_rebuild_phenotypes.csv` | master regulator activation frequency (%), four miRNA conditions |
| `python_rebuild_metadata.json` | its run parameters and node lists |

That run used the same rule table as `model/thesis_table5/` but a different
engine and protocol: 400 trajectories, 10 sweeps, activation frequency at the
final state rather than occupancy over a window. **The magnitudes are not
comparable to this repository's** — only the directions are, and that is all
`scripts/05` compares them on.
