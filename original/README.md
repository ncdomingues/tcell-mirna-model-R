# The 2017 files, untouched

Nothing in this directory has been edited. It is the working material from the
thesis directory as it stood, kept for provenance and because the commentary in
it is part of the record.

- `base_files/GinSimR_ND.R` -- the original simulation engine, plus the
  Naldi and Abou-Jaoude driver code. `R/ginsim_engine.R` is the documented
  refactor of the `getordervar` / `ginsimnxt` / `ginsimrun` functions here.
- `base_files/method_modelling.txt` -- the method note: table format, the basal
  value convention for unregulated nodes, and the
  `number of items to replace is not a multiple of replacement length` warning
  that turned out to be a real bug (MODEL_NOTES.md, repair 1).
- `2017_model/` -- the simulation and plotting scripts for my own network,
  and `Procedures.txt`, which records the initial state and the 100- versus
  200-time-point runs.
- `stimulations/` -- the IL2, Th1 and Th2 conditions as originally scripted,
  with their hand-counted position vectors.
- `mir34c_tf_variations/` -- the plotting scripts from the February and March
  2017 promoter work, plus the `formulae.txt` and `nodes_values.txt` that
  directory carried, which differ from the published table.

The `.RData` and `.Rhistory` files from those directories are deliberately not
included: they are large, binary, and hold no information the tables and scripts
do not.
