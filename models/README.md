# Run-length calibrations

`revio.run-model.tsv` is the existing LA estimator. The two corroboration models
are byte-identical to the frozen HG002 models used in the native comparisons:

| File | SHA256 |
| --- | --- |
| `lao-corroboration.run-model.tsv` | `2ca8148df1d222f3b011c85858aa538f9c76dcf824303deeec677cd5c4271cb7` |
| `ul-corroboration.run-model.tsv` | `cba574b3c1c779cdc24a2ece937485062d8c1920edee7533c9c44782e46a8af7` |

Supply the latter with `--experiment-corroboration-model`, **not** `--tier-model`:
they resolve uncertain first-tier calls without replacing ordinary lower-tier
modal estimation. Tier labels must match the command's `--tier` labels.
Calibration is dataset/chemistry/basecaller-specific, not a universal ONT model.

Runtime requires only the TSV files; no metadata directory is needed. Detailed
training provenance remains outside Git in the workspace's
`hg002-test/calibration20/models/` and
`hg002-test/ul-priority-calibration/models/`. No genomic datasets accompany these
small parameter files. See `experiments/README.md` for evaluation limitations.
