# Maintenance validation — 2026-10-14

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
