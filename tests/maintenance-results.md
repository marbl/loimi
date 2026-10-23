# Maintenance validation — 2026-10-14

## Alignment performance experiment — 2026-10-21

Separate Release executable: `build-perf/hpc-consensus`. The original
`build/hpc-consensus` was not rebuilt; its SHA256 remains
`1d542f3881d318674eb0224cd227d255f3169e9a8e9f972d97e825d8e27ea4f6`.
Changes: bound initial raw Edlib search by the existing valid raw trace's edit
count; parallelize final-reference CRAM realignment in bounded batches with
ordered serial emission. No voting or scoring policy changes.

Four-thread polishing+CRAM comparisons on the 224 Mb-site fixture:

| Order | Before wall | After wall | Before OS peak RSS | After OS peak RSS |
| --- | ---: | ---: | ---: | ---: |
| Before, after | 21.91 s | 14.86 s | 296.2 MiB | 297.5 MiB |
| After, before | 25.54 s | 15.55 s | 310.1 MiB | 283.7 MiB |

Wall time decreased 32–39%; user CPU decreased from 65.61/66.87 seconds to
43.26/43.29 seconds. Other quality runs and verification work were active, so
these are preliminary measurements under contention, not isolated performance
bounds. This fixture has zero accepted polishing edits and primarily exercises
the tighter bound, not final-reference realignment. All text artifacts and
decoded SAM records match byte-for-byte; CRAM integrity checks pass.
The first pair predates the cleanup that also releases failed/cropped traces;
the reverse-order pair uses the final executable (hash recorded in its report).

Release and ASan/UBSan suites each pass all 12 tests. New coverage includes
1,000 seeded random bound/trace comparisons, a 20-kb noisy traceback, unchanged
error-rate rejection, and native serial/parallel CRAM parity with both one-read
and multi-read batches. These checks do not establish whole-chromosome QV or RAM.
Artifacts: `../hg002-test/perf-alignment/` relative to the repository root.

## Earlier maintenance baseline

Baseline: the executable immediately before this maintenance pass, **not** the
mature Verkko caller or the older high-memory prototype. Same native libraries,
Release optimization, three ordered tiers, LA Revio model, four threads.
Runs were sequential; the 224 Mb-site comparison was repeated in reverse order.
These are local macOS measurements, not full-chromosome benchmarks.

| Fixture / mode | Before wall time | After wall time | Before peak RSS | After peak RSS |
| --- | ---: | ---: | ---: | ---: |
| Terminal 10.3 Mb, fast consensus | 167.69 s | 170.15 s | 2958 MiB | 2900 MiB |
| 224 Mb-site, 0.72 Mb, polish + CRAM | 12.64 s | 12.50 s | 282 MiB | 289 MiB |
| Same fixture, reversed execution order | 12.47 s | 12.38 s | 291 MiB | 300 MiB |
| 199 Mb-site, 1.09 Mb, polish + CRAM | 16.33 s | 16.59 s | 301 MiB | 307 MiB |

Peak RSS is the OS `time -l` measurement, distinct from the caller's intermediate
high-water logs and from macOS “peak memory footprint.” The first two fixtures
were measured before the final build-system cleanup; the final 199 Mb-site
comparison used the final rebuilt binary. The worker exception-propagation
change is included in the terminal and 199 Mb-site runs.

Interpretation: the larger fast fixture used about **2% less peak RSS**, but ran
about **1.5% slower**. Small polishing/export fixtures showed approximately
unchanged runtime and slightly higher OS peak RSS. There is **no demonstrated
end-to-end speedup or universal peak-memory reduction** from this pass. Small
differences need more repetitions before attributing them to a particular code
change. Alignment work remains dominant; its scoring/ties were not changed.

The code does eliminate large transient reference copies during edit application
and draft FASTA writing, transfers the draft buffer when unchanged, avoids
per-shift insertion-allele allocation, and removes redundant eligibility checks.
These are local allocation savings; they do not necessarily affect a peak set
earlier by HPC evidence or later by samtools. For a 244 Mb reference, removing a
single full-sequence temporary avoids approximately 244 MB at that particular
stage, not a demonstrated 244 MB reduction in process peak RSS.

## Correctness

- All four CTest suites pass in Release and in a separate Debug build with
  AddressSanitizer and UndefinedBehaviorSanitizer. Linked prebuilt native code
  is not instrumented.
- All three biological fixtures have byte-identical final FASTA and diagnostic
  outputs relative to baseline. Both polishing fixtures also have identical
  draft FASTA, edit VCF, and decoded SAM records (including order and tags).
- CRAM integrity checks pass. Existing native integration tests additionally
  validate final-reference NM/MD using samtools calmd.
- Identical sequence means identical input to QV evaluation; meryl was not
  rerun for this behavior-preserving cleanup. Full-chromosome QV/performance
  has not been remeasured.
- Resource exhaustion is intentionally different: worker allocation failures
  now fail the run after joining workers, rather than silently dropping reads.

Local reproducibility artifacts are in `../hg002-test/maintenance-comparison/`
relative to the repository root: `tail-fast`, `site224-final`, `site224-repeat`,
and `site199`. Each contains logs and `comparison.json` with commands and
matching artifact lists. The final comparison also records executable hashes.
Datasets and generated outputs are not part of the repository. Reproduce with
`tests/compare_versions.py`; see `MAINTAINING.md` for usage and caveats.
