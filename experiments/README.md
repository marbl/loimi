# Isolated consensus experiments

## Terminal homopolymer filter

`--experiment-terminal-hp-filter` excludes query-terminal run observations
from interior template positions before evidence accumulation and tier-depth
qualification. It applies equally to every tier and placement class, including
borrowed evidence. It does not remove reads from alignment export or change the
polishing vote policy. `UNFLANKED_RUN` identifies excluded observations in the
per-read diagnostic. Exclusions contribute to `filtered_hpc_votes`.

The small policy is isolated in `src/experiments/hp_anchors.hpp`; removing its
CLI option and guarded call from `main.cpp` removes the experiment. Default
behavior remains unchanged. This first-stage experiment only checks that query
runs exist on each required side, NOT that a multi-base alignment anchor is
good. Preserve the legacy outward-facing template-end behavior and singleton
passthrough; do not interpret results as validation of a new telomere policy.
Clipped partial-run protection and local HPC alignment filtering still apply.

`hp-anchors` tests interior/boundary contracts. `native-anchors` checks both
query-terminal directions, before-tier-selection depth changes, fallback to
the next tier, and the diagnostic status in real synthetic cnspacks. The
biological test compares normal polished outputs across the ten frozen windows;
the same-tier borrowing flag is enabled in BOTH arms. Run-length-only sampled
evidence replay can additionally assess absence of LAO, but is not an end-to-end
LA/UL assembly validation and does not apply borrowing.

## Native per-read evidence tracing

For a small set of `TIG RUN` coordinates, add `--evidence-sites sites.tsv
--evidence-read-votes`. The extra `.run-read-votes.tsv` streams one row per
aligned placement/reference run, including tier, declared/resolved placement
count, orientation, HPC span, alignment operation, observed run length, and
eligibility status. It does not change consensus votes or retain a chromosome
of per-read diagnostics. Histogram and per-read exports include lower tiers
even when a higher tier has already determined the call.

`ELIGIBLE` means the observation passed evidence checks, not that its tier won.
`DIAGNOSTIC_EXCLUDED` means the explicit diagnostic read-ID exclusion omitted
it from histograms, not from consensus. Other statuses distinguish nonmatching
operations, N template bases, clipped partial runs, and local alignment filters.
Insertions consume no template run and have no separate row. Whole-read failures
are in `.read-decisions.tsv`; reads not spanning a site have no site row.
Query raw coordinates are zero-based half-open in the oriented full input read;
HPC coordinates are zero-based half-open in the supplied template. For deletions,
length is zero and the query interval is empty. Native regression coverage
checks histogram counts, coordinates, and unchanged consensus with tracing on.

Branch: `codex/local-consensus-experiments`. Existing uncommitted work was
preserved when creating this branch; branch creation is not a snapshot of those
edits. No dataset or experimental output belongs in Git.

## Boundaries and removal

- `src/experiments/borrowing.cpp` owns the bounded borrowing calculation. The
  default model reader and length estimator remain in `src/models.cpp`.
  Borrowing is disabled unless `--experiment-consistent-ambiguous` is supplied.
  Removing it completely also requires removing its three declarations in
  `RunModel`, CLI option, guarded evidence/audit hooks in `main.cpp`, CMake source
  entry, and borrowing tests. Disabling the flag requires none of those edits.
- The local graph is a **diagnostic executable**, not a production dependency.
  It is built only with `HPC_BUILD_GRAPH_PROBE=ON`. Its upstream graph source
  remains in the external Verkko checkout with its original license. Set the
  option OFF to exclude the graph; no caller changes are needed.
- Python experiments here never edit assemblies or train production models.
  Their outputs go to an explicitly supplied directory. They can be removed
  independently, along with the `experiment-classification` CTest entry.
- Additive retention remains a separate opt-in policy. Do not combine it with
  borrowing or interpret results from one policy as validation of the other.

## Sequence-context classification

`classify_remaining.py --native-root PATH --output FILE.json` reads the paired
ten-window artifacts. It excludes overlapping reference blocks in either arm
and classifies remaining borrowing-arm differences by HP association and exact
short tandem repeats. Classification is descriptive, not a mechanistic verdict.
`test_classification.py` provides data-free boundary and repeat tests.

## Matched-read local graph

Run `local_graph.py` with `--native-root`, `--input-root`, `--probe`, `--model`,
and `--output`. The input directory needs `part002_fast.fasta`, its matching
indexed `part002_fast.cram`, and `part002.fasta` (mature consensus). The probe
argument points to `build/hp-graph-probe`; its sibling `apply-run-model` is used
for the unchanged LA model. The native root needs the classification above.
Use `--require-mature-agreement` for a second, explicitly selected diagnostic
subset where Q100 and mature Verkko agree. Keep this separate from the broader
set, which includes differences shared by both callers.

The script selects up to two anchored remaining HP errors and two nearby correct
controls per window. For each, it extracts approximately 1 kb read chunks from
unique primary CRAM placements, using the uploaded assembly as the alignment
reference. Selection never filters reads by agreement with truth. It compares
the LA-calibrated or fallback modal length call and Verkko's local graph using
identical HPC-eligible reads. The second arm removes LAO before selection, using
LA then UL. All selected read names, local references, graphs and histograms
are retained for audit. The graph diagnostic's sparse fallback is the first
nonempty tier, not the last tier.
Chunks exceeding the alignment error cutoff are recorded as `REJECT` and
excluded from both estimators. Forward/reverse histograms and read strand
assignments are retained; they are diagnostics, not a new strand-based filter.

These truth-selected, exact-flank-resolvable cases are diagnostic, not an
unbiased quality estimate. Spanning chunk requirements alter local coverage;
MAPQ60 proxies native uniqueness; no ambiguous-read fallback or full Verkko
read-promotion replay is attempted. The frozen three-tier LA model is retained;
UL uses mode, not an unvalidated new calibration. Correct controls are essential
because a graph can fix one allele while breaking another.

Do not integrate a local-graph correction into production solely because this
diagnostic improves. First define a truth-independent trigger and acceptance
rule, test unchanged controls and absent-LAO data, then validate end-to-end on
additional held-out windows. Keep any future backend behind its own opt-in flag
and separate translation unit; do not add graph-specific branches to the
ordinary run-length estimator.
# Aligned-query flank experiment

`--experiment-hp-flank-bases N` requires N raw bases on each side of a run
within the oriented, aligned, clipped query and immediately adjacent matching
HPC runs. It does not require perfect matching throughout those N bases.
Only the outward-facing condition is waived at a supplied template endpoint;
end extension and singleton handling retain their separate policies.
The filter changes vote eligibility before tier qualification, not alignment
export. Existing local HPC error filtering still applies. The default is off.

The ten-window flank/retention experiment uses this option at 1,000 and 3,000
bases, additive higher-tier unique retention in every arm, and no borrowing.
Both LA/LAO/UL and LA/UL are run through native consensus and polishing.
It is an experiment, not a validated default. Run-distance validation lives in
`tests/hp_anchors.cpp`; `tests/anchor_native.py` exercises native tier fallback.

