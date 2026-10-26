# Isolated consensus experiments

## Boundary-aware flank filtering (2026-10-26)

The 500-bp filter now counts usable query overhangs removed by template-boundary
alignment when measuring distance to a read end. Previously it measured only
the template-clipped query and erased terminal run evidence, which could make
all extension anchors fail. The outward exact-template-run waiver is unchanged.
Overhangs contribute distance only: adjacent HPC matches, local alignment-quality
filtering, partial-run protection, tier ordering, and extension support thresholds
are unchanged. Packaged/trimming clips are excluded; internal alignment endpoints
get no bonus. Region-mode artificial right prefix boundaries get no bonus either.
The change applies even with end extension disabled, so the terminal draft itself
can be called correctly. It does not exempt a broad terminal region from filtering.

`native-end-flanks` tests forward/reverse boundary reads, both/left/right/disabled
extension, 500-bp filtering of high-priority short interior reads, complete run
recovery, and CRAM validation. Unit tests exercise exact distance thresholds and
missing adjacent anchors. The paired biological tail test is kept outside Git
in `hg002-test/end-flank-regression/`; results there distinguish synthetic
both-end guarantees from recovery of the real right telomere.

Completed biological regression with 500-bp flanks and HP realignment enabled,
corroboration disabled to match the reported run: right extension 0 -> 1,934 bp,
79 anchored molecules. All 310 unsupported right-terminal runs regained evidence;
315 changed run records are confined to the final 316 HPC runs. Interior run
evidence is identical. Total unsupported fixture runs fell 705 -> 395; this does
not claim resolution of the artificial left fixture boundary. `seqtk telo` changed
from no terminal call to a 2,080-bp right telomere. CRAM quickcheck and calmd passed
with 34,151 records. All 16 Release and sanitizer suites pass, and the old binary
fails the new native regression. Full-chromosome QV and actual left-telomere
recovery were not retested. Full results: ignored `end-flank-regression/validation.json`.

## Full-chromosome trial with linked-event grouping

The combined HP policies and linked-event grouping are present on this branch.
Use a fresh Release build; the historical `build/` executable was frozen for
older comparisons and should not be assumed current. For example, from this repo:

```sh
cmake -S . -B build-full-chromosome -DCMAKE_BUILD_TYPE=Release -DVERKKO_ROOT=/path/to/verkko
cmake --build build-full-chromosome -j 8
ctest --test-dir build-full-chromosome --output-on-failure

./build-full-chromosome/hpc-consensus --threads 32 \
  --import /path/to/part002.cnspack \
  --template-hpc /path/to/part002_template.fasta \
  --tier 'LA=hifi_m*' --tier-model LA=models/revio.run-model.tsv \
  --tier 'LAO=hifi_[0-9a-ln-z]*' --tier 'UL=ont_*' \
  --experiment-corroboration-model LAO=models/lao-corroboration.run-model.tsv \
  --experiment-corroboration-model UL=models/ul-corroboration.run-model.tsv \
  --experiment-hp-realign --experiment-hp-flank-bases 500 \
  --polish-non-hpc --output /path/to/part002_trial \
  --cram /path/to/part002_trial.cram
```

Add `--tig-map /path/to/packages.tigName_to_ID.map` if template names need mapping
to packaged tig IDs. Do **not** pass `--experiment-master-events`: linked-event
grouping is the branch default, and that flag disables it for baseline controls.
Do not add retention or borrowing to this tested combination. The example keeps
normal end extension; add `--no-extend-ends` to match the internal-window tests.
CRAM is optional. Neither model metadata nor bulk experimental outputs need to
be committed. Include the new source/header/test/model files in the commit;
avoid blanket staging of local name maps or benchmark data.

The three HP policies improved aggregate differences on ten historical windows
with master event extraction. Linked grouping has separate targeted validation;
the exact four-change combination has native synthetic/CRAM regression coverage,
but its full-chromosome accuracy and resource costs remain to be measured.

## Combined native follow-up (2026-10-23)

`hg002-test/run_combined_native.py` runs three arms on the same ten historical
validation windows: corroboration plus HP realignment, 500-bp flanks alone, and
all three. `build-combined-native` preserves the previous experiment's frozen
binary and outputs. Each window first reproduces the frozen master baseline;
the prior baseline sequences/calls are reused with recorded hashes. Master is
not modified. LA/LAO/UL, LA calibration, ordinary polishing, master event
extraction, no retention, and no end extension remain fixed.

Composition is explicit: initial corroboration precedes expansion. The HP pass
then rebuilds evidence and reapplies the same corroboration rule only when the
first unique tier qualifies. This avoids automatically reverting corroborated
lengths to the first-tier-only estimate, while allowing genuinely changed local
evidence to change the decision. First-stage votes are never counted twice.
With 500-bp flanks enabled, the rebuilt HP observations must also have >=500
aligned raw query bases on both sides; the HP pass cannot bypass the end filter.
The local realignment context remains 64 bases per side. No thresholds are
trained or selected from this validation run. Borrowing/retention/mean/terminal
filter combinations remain rejected; the three requested policies can compose.

Unit tests cover corroboration after rebuilding, rejection of non-unique lower
votes for corroboration, and raw-query flanks. Native fixtures cover both
combinations and 500-bp flanks, including CRAM. Results go only to ignored
`hg002-test/combined-native/`. The final shared mask includes the previous four
independent policies, permitting comparison of combined and individual benefits
on identical truth positions. Changed-locus audits mark unresolved flanks rather
than treating them as successes. Two concurrent four-thread jobs bound active
genome tasks; runtime/RSS remain descriptive rather than controlled benchmarks.

## Independent native comparison (2026-10-23)

The workspace runner `hg002-test/run_four_native.py` compares four independent
arms to the **current hpc-master executable**, frozen by SHA256. The ten historical
validation windows use LA/LAO/UL, the existing LA Revio model, ordinary non-HPC
polishing, no end extension, and no additive retention. All four active arms use
`build-four-native/hpc-consensus` after a per-window unmodified-output parity
check against master. The initial trial found that master accepts borrowing/flank
flags but has their implementation inactive; those candidate outputs are explicitly
discarded. Only completed, hash-verified master baselines are reused. The branch's
linked-event normalization is another difference: `--experiment-master-events`
selects the frozen master extractor in `src/experiments/master_events.*` in all
experimental arms, including the parity control. These are not pristine
holdout windows. Nothing changes production defaults or the master checkout.

1. `--experiment-consistent-ambiguous`: existing bounded same-tier borrowing.
2. `--experiment-hp-flank-bases 1000`: existing aligned-query span requirement.
3. `--experiment-corroboration-model LABEL=FILE`: frozen replay rule, limited to
   qualified unique first-tier calls with model top-two log-likelihood gap <=3.
   Consult the first lower unique tier with >=5 observations and model coverage
   for both candidates. Switch only if per-read log-likelihood advantage >=0.5
   and direct alternate support >=0.6. Stop after that tier, including on a veto.
   These models **never** replace ordinary lower-tier modal calls. Do not pool
   corrected and raw ONT. The experiment requires support five and a first-tier
   model, not any particular read-name label. `.corroboration.tsv` audits decisions;
   `UNIQUE_CORROBORATED` marks changed run calls. Depth/peak remain original-tier
   histogram diagnostics, not corroborating alternate counts.
4. `--experiment-hp-realign`: one-pass **re-estimation after local realignment**,
   not a new candidate-sequence likelihood model. Select qualified runs where
   the existing model differs from mode or the existing uncertainty predicate
   fires (model LR <100; unmodeled modal fraction <90%). Locally align each
   spanning read against the expanded draft with 64 raw bases on either side,
   requiring 12 exact terminal matches and matching immediate non-run flanks.
   Rebuild complete pure-run length histograms, using unique tiers first, then
   same-tier pooled evidence, with the existing models/mode and depth threshold.
   No sparse post-polishing edits. Never add these observations to initial votes.
   Limits: run <=1000 bases, query chunk <=4096 bases, local error cap 0.15;
   reference/read ends lacking full context are skipped. Orientation is already
   normalized by the caller. Exact anchors are **not** uniqueness guarantees.
   Non-HPC edits win conflicts. Accepted run replacements are included in the
   existing edit VCF (historically named `.non-hpc-edits.vcf`) for coordinate
   lifting, and separately audited in `.hp-realignment.tsv`; the summary's
   `non_hpc_edits` excludes them. Final CRAM uses the usual lifted realignment
   path. Fast expansion-only mode is unchanged; this flag requires polishing.

New policies live in `src/experiments/refinement.*`, with unit and native fixture
tests (including CRAM when samtools is available). The independent study used
one policy per arm; supported compositions are described above. Histograms are sparse by uncertain target; the
local pass retains one read's decoded trace/query-coordinate map at a time.
There is still O(number of candidate sites × represented tier histograms) evidence
storage; it is not a hard chromosome-memory cap. Local alignment is currently
serial so its performance cost can be measured before optimization.

The runner uses at most two concurrent four-thread consensus jobs. Each window
runs master, a parity control, then the four independent candidates. Each is scored on a
paired depth-one truth mask. `score_four_native.py` additionally intersects the
masks across **all five arms**, reporting per-window differences, elapsed job
seconds, and peak RSS. Timings from concurrent jobs are not controlled speed
benchmarks. Pair-specific results cannot be compared as though their denominators
were identical. Fix/regression audit must be interpreted separately from aggregate
base-count changes; a single maternal truth does not establish switch accuracy.

Outputs are in ignored `hg002-test/four-native-v2/`; previous experiments are
preserved, and the invalid initial trial is marked in `four-native/DISCARDED.txt`.
The new policies are hypotheses under test, not validated improvements.

## Corroboration, confidence, and local-anchor pilot (2026-10-22)

`corroboration.py` replays the frozen 20-window filtered evidence, with additive
retention, without borrowing or new flanking. Ten training windows select a
rule before ten historical validation windows are opened. It consults only the
first lower tier with >=5 unique observations and model coverage for LA's top
two candidates. Candidate switches require a small LA log-likelihood gap, a
lower-tier per-observation likelihood advantage, and direct alternate support.
The grid is LA gaps {1,3,10}, lower advantages {0.1,0.5}, support fractions
{0.6,0.8}; training minimizes inverse-sampling-weighted wrong runs, breaking ties
by fewer changes. Baseline is a candidate. No LAO+UL pooling is permitted.
Each tier's likelihood uses its own frozen model. Non-LA baseline calls remain
modal; the corroborating model does not silently become a new baseline caller.

The same pass measures likelihood-based confidence and training-fitted error
rates by margin, base, depth, called length, and model/mode disagreement.
Fine bins with <100 training sites back off to margin-only bins; error rates
use Beta(1,1) smoothing. Report Brier scores and actual/expected errors, not
uncalibrated likelihood ratios as posterior confidence guarantees. Model fitting
and confidence fitting share historical training windows; their training error
is optimistic. Validation windows were already inspected in prior work.

`local_anchor_analysis.py` collects native per-read diagnostics from training
w00/w02 and validation w01/w11 using a separate analysis binary. It samples
one tenth of existing truth-labelable loci via a fixed integer hash and includes
1-kb-changed labeled sites in a separately reported diagnostic subset. That
subset must not be interpreted as unbiased accuracy. The replay tests 2,5,10,20
consecutive matching HPC runs on each side against baseline and 1-kb aligned
span. Native baseline call parity is mandatory before scoring. It reports tier
switches, loss of LA unique depth5, strand-mode disagreement, fixes and breaks.
Training selects a policy before validation analysis. These are local alignment
anchors, not proof of unique genomic placement. Polishing is not replayed.

`--evidence-read-votes` now additionally reports consecutive left/right matching
HPC runs (capped at32) and aligned raw-query bases on each side. Collection
does not change votes or consensus and bounds work per requested site. The
separate `build-quality-analysis` leaves the frozen `build/hpc-consensus` intact.
Data-free replay checks: `python3 experiments/test_quality_analysis.py`.

Artifacts live outside Git under `hg002-test/corroboration-confidence` and
`hg002-test/local-anchor-pilot`. Short truth runs were originally sampled 1/100;
weighted counts are estimates, not whole-sequence QV or phase-switch accuracy.
No production defaults change, and promising replays require native validation.

Broader local-anchor replay: `local_anchor_analysis.py --validation-all
--all-labels` evaluates all existing labels in the ten historical validation
windows, without the pilot's extra subsampling. Baseline, 1-kb span, and fixed
10/20-run local anchors are the primary comparisons; 2/5-run results remain
secondary diagnostics. No threshold is selected from validation. Reports split
the two pilot validation windows from the eight additional windows, and retain
per-window fixes, breaks, tier changes, and mandatory native baseline parity.
Outputs use a fresh `hg002-test/local-anchor-broad` directory; previous pilot
results and all production defaults are preserved. This remains run-length
replay, not an end-to-end polishing or phase-switch test.

First completed corroboration replay: training chose LA gap<=3, lower-tier
per-observation advantage>=0.5, alternate fraction>=0.6 for LA/LAO/UL; LA/UL
selected baseline. On 734,161 validation labels, LA/LAO/UL changed 73 calls:
51 fixes and 22 breaks, reducing wrong runs 307->278 (weighted estimate
1396->1367). LA/UL retained 1003 wrong runs with additive retention. No window
worsened in aggregate for the selected three-tier rule, but individual
regressions preclude a regression-free claim. This is not native polished QV.

Confidence audit covers 691,430 model-eligible unique-LA validation loci, with
104 errors. Likelihood-normalized confidence predicts 57.49 errors; the simple
training-fitted bins predict 133.84. Brier score worsens slightly from
0.000118127 to 0.000119654. Thus these fitted bins are not a validated confidence
improvement. Model/mode disagreement concentrates 40 errors in 1,151 loci;
agreement contains 64 errors in 690,279 loci. This is a useful diagnostic
stratum, not grounds to replace the model with mode. Local-anchor pilot results
remain pending in the separate artifact directory.

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

