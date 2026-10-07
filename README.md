# hpc-consensus

Standalone C++ prototype for expanding a fixed homopolymer-compressed (HPC)
template using reads and placements from Verkko `.cnspack` files. It reuses the
native Canu/Verkko package reader; it does not use POA or remap reads genome-wide.

Status: functional prototype, October 6, 2026. Synthetic native-package tests
pass; biological accuracy and genome-scale resource use have not yet been
validated. The existing Verkko pipeline is not modified.

## Build

Requires a C++20 compiler, CMake 3.18+, pkg-config, Python 3 (test configuration),
and a built Verkko checkout containing `build/lib/libverkko.a`. Link dependencies
are OpenSSL, libcurl, liblzma, zlib, bzip2, threads, and the OpenMP runtime used by
the Verkko build. Use a compatible compiler, architecture, and store ABI.
`samtools` on PATH is required only for CRAM export and its integration tests.

From this directory:

```sh
cmake -S . -B build \
  -DVERKKO_ROOT=/Users/korens/Documents/Software/verkko \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
build/hpc-consensus --help
```

The executable is `build/hpc-consensus`. It is standalone at execution time
apart from its dynamically linked system libraries and, for CRAM, `samtools`.
Building still requires Verkko headers and its static library.

## Inputs

- `--import`: one native `.cnspack`, potentially containing multiple tigs.
  Reads, stored names, orientations, clipping, and HPC placement intervals come
  from this file. A separate read FASTA or sequence store is not needed at runtime.
- `--template-hpc`: plain FASTA of HPC path sequences. The current package writer
  does not embed these sequences, so this companion file is normally required.
  Embedded tig bases are accepted as an HPC template if present. If both sources
  are supplied, they must agree.
- Input template names are `piece000001`-style names or numeric packaged tig IDs.
  All output references use `piece` plus the tig ID padded to at least six digits,
  consistently in FASTA, CRAM, VCF, and diagnostics. No tig-name map is needed;
  `--tig-map` has been removed. Arbitrary graph segment names are not automatically
  associated with numeric tigs; scfmap-based input association is not implemented.
- Templates must be A/C/G/T/N, homopolymer compressed, and correspond to the
  packaged tig. Their offset is inferred from reads as described below.
- `--tier` specifies name globs, highest priority first. Patterns match the
  complete name stored in the package; the first match wins. Some package writers
  prepend library names, so inspect the actual packaged naming convention.

Native packages are ABI-dependent, trusted inputs, not a portable or hardened
exchange format. Use packages compatible with the linked Verkko build.

## Run

Fast uncompression, with illustrative name patterns (replace them for your data):

```sh
build/hpc-consensus \
  --import part.cnspack \
  --template-hpc paths.hpc.fasta \
  --tier 'hifi=hifi_*' --tier 'ont=ont_*' \
  --threads 4 --output results/fast
```

Create the output parent directory first. Quote globs to prevent shell expansion.
Unmatched eligible reads cause an error by default; use `--unmatched drop` to
discard them, or add an explicit last tier such as `--tier 'other=*'`.
Bare patterns are allowed and receive automatic labels.

Correction and alignment output are independent opt-ins. Append either or both:

```sh
--polish-non-hpc
--cram results/final.cram
```

Use a fresh output prefix for each run. Existing output files are refused.
The default does not perform raw-base realignment, correction, or CRAM creation.
It can still align compressed read sequences within their known layout intervals
to establish run correspondence when they do not exactly match the template.

## Calling rules

Template offsets are inferred independently for each tig. The tool tries up to
32 uniquely placed, eligible reads, ordered by tier and then proximity to the
layout start, preferring reads away from possible overhangs. Full compressed
reads are aligned to the template; tied optimal placements are rejected.
With six or more eligible reads, six placements are required, split into discovery and validation sets.
Each set must support a consistent offset cluster (at least three reads and
75% of estimates within a 20-HPC-base span); the two estimates must agree
within 20 bases. The combined cluster's lower median is used. This is a
heuristic consistency check, not a calibrated confidence probability.
Unresolved tigs fail explicitly. Reads crossing the inferred template boundary
are aligned semi-globally to their overlapping HPC template interval. Overhangs
are clipped at aligned run boundaries and represented as soft clips in CRAM;
they do not extend or invent template sequence. Left, right, and two-sided
overhangs are supported. Ambiguous best endpoints or excessive alignment error
produce `HPC_ALIGNMENT_FAILED`; reads wholly outside the template produce
`OUTSIDE_TEMPLATE`. Existing raw clipping is retained. If it cuts through a run
(the adjacent clipped base is identical), that partial run cannot vote for a
run length. The old global `--template-offset` option has been removed.
`PREFIX.template-offsets.tsv` records offsets, attempts, placements, validation
support, tolerance, and method (`VALIDATED`, `SPARSE`, or `SINGLE_READ`).
With fewer than six eligible reads in a multi-read tig, every eligible read is
aligned and the available unique placements are clustered without held-out
validation. At least one placement is required; inconsistent estimates still
fail the 75% cluster agreement check. `SPARSE` explicitly indicates weaker
evidence, not validated confidence.

A tig containing exactly one packaged read bypasses inference and consensus:
its oriented, clipped sequence is output directly, without requiring a template.
It must still be uniquely placed and match an accepted tier. Run support
thresholds do not apply to passthrough. Optional CRAM preserves the original
oriented read and soft clips, with a perfect match to the emitted sequence.
Polishing leaves singleton sequences unchanged. No offset is inferred (reported
as `.`), and run-call coordinates describe the passthrough sequence.

Only placements with packaged placement count exactly one and no distinct
alternative in the input package are eligible for consensus and offset inference.
Ambiguous reads are excluded from consensus, but their distinct placements are
aligned and exported when CRAM is requested. The legacy `skipConsensus` flag is ignored:
explicit user tiers replace its template-building priority policy. Any number
of ordered tiers is supported. All eligible reads can be aligned and exported
to CRAM, even if their tier never wins a consensus call. Repeated IDs are allowed
within and across tigs. Identical entries (ID, coordinates, strand, and clipping
within a tig) are deduplicated. A package prepass counts distinct occurrences.
Uniqueness is trusted from the layout, not independently inferred by mapping.

For each template run, only matching aligned read bases contribute a run length.
The first tier with at least `--min-run-support` observations wins; lower tiers
cannot outvote it. The called length is the empirical histogram mode. Tied modes
are resolved by minimum total absolute distance to observations, then smaller
length. This is not a fitted or smoothed distribution model. Quality priority is
entirely user-specified by tiers, not per-base quality scores.

Without eligible evidence a run gets length one and is reported as
`NO_EVIDENCE`; `--missing-run error` instead stops the run. N runs do not collect
run-length evidence. Fast-mode output preserves the input HPC base sequence.

Optional correction calls substitutions and short non-HPC events from raw traces.
It selects the first tier with sufficient spanning depth and requires alternate
support and fraction thresholds. If both strands span the event in that tier,
both must support it. Terminal events, pure run-length changes, overlapping or
adjacent accepted edits, and alternate alleles containing N are excluded.
There is one correction pass, not iterative polishing.

## Options and defaults

| Option | Default | Meaning |
| --- | --- | --- |
| `--threads N` | 1 | Parallel per-read alignment work within a tig; also passed to samtools |
| `--batch-bases N` | 8000000 | Approximate raw-base budget for each HPC alignment batch; a single read may exceed it |
| `--min-run-support N` | 1 | Minimum total observations for a tier to call a run |
| `--missing-run one\|error` | one | Policy when no tier qualifies |
| `--unmatched error\|drop` | error | Policy for eligible reads matching no tier |
| `--max-error F` | 0.15 | Unit edit-distance cap, rounded up, relative to the longer sequence |
| `--polish-non-hpc` | off | Enable conservative correction |
| `--min-polish-support N` | 5 | Minimum spanning depth and alternate count in the chosen tier |
| `--min-polish-fraction F` | 0.9 | Alternate fraction; must exceed 0.5 |
| `--max-polish-indel N` | 500 | Maximum reference or alternate event length, including complex events |
| `--cram FILE` | absent | Write sorted, indexed, embedded-reference CRAM |

## Outputs

For `--output PREFIX`:

- `PREFIX.fasta`: final expanded consensus, corrected only if requested.
- `PREFIX.run-calls.tsv`: zero-based run index and half-open **draft** base
  coordinates, base, called length, winning tier, tier depth, and modal count.
- `PREFIX.read-decisions.tsv`: tig, packaged read ID/name, assigned tier, and
  status (`ACCEPTED`, `AMBIGUOUS_PLACEMENT`, `UNMATCHED_TIER`,
  `HPC_ALIGNMENT_FAILED`, or `FINAL_ALIGNMENT_FAILED`). `ACCEPTED` does not mean
  the read's tier won at every run.
- `PREFIX.summary.json`: tig/read counts, accepted HPC evidence reads, dropped
  reads, exact HPC matches, unsupported runs, edits, raw alignment calls, and
  CRAM records. `consensus_reads` counts successful HPC evidence reads even if
  their later raw alignment fails; it need not equal the final CRAM count.
- With correction: `PREFIX.draft.fasta` and `PREFIX.non-hpc-edits.vcf`.
  VCF positions are one-based against the saved **draft**, not the final FASTA.
  Run-call coordinates also remain in draft space after correction.
- With CRAM: the requested CRAM, its `.crai`, and `PREFIX.fasta.fai`.

CRAM contains at most one primary alignment per accepted packaged read, with
stored names, orientation, soft clipping, tier RG, and final-reference NM/MD
tags. The first successfully emitted occurrence in package/layout order is
primary; later occurrences are secondary (flag 256), regardless of output sort
order. MAPQ is 60 for unique placements and rounded `-10 log10(1 - 1/n)` for
`n > 1` equally good placements (2:3, 3:2, 4:1). `n` is the greater of the
stored placement count and distinct package occurrences, before alignment
filtering. Unknown stored multiplicity with one observed placement gets MAPQ 0.
These are multiplicity-based conventions, not calibrated alignment confidence.
Primary selection is scoped to one invocation; independently produced packages
must reconcile primary flags when merged. QUAL is unavailable, and AS is negative unit edit
distance, not an affine-gap score. Failed alignments are omitted and reported.
References are embedded so decoding does not depend on the temporary FASTA.

Outputs are staged in `PREFIX.work.*`; failures retain this directory for
diagnosis. Publication of multiple files is not transactional. Keep the CRAM
destination on the same filesystem as the prefix for rename-based publication.

## Tests and E. coli handoff

See [tests/README.md](tests/README.md) for automated test coverage and the small
bundle needed for biological testing. The first E. coli test completed:
automatic inference recovered offset 5608 and reproduced the manually aligned
fast consensus exactly (99.9951% FastANI versus 99.9985% for the existing assembly;
the prototype test covers only the main tig). This is not production validation.

## Prototype limitations

- Processes one tig at a time, retaining decoded read sequences, compact traces,
  and per-run evidence. HPC alignment scratch space is batched, not retained for
  all reads. Total memory is still not a hard-capped/windowed implementation.
- Uses unit-cost edlib global alignment within known intervals, not affine-gap
  alignment or remapping. Raw repair currently realigns the whole read layout
  interval if its HPC trace has errors, and again if correction touches it.
- Raw run-expanded alignments are valid but are not guaranteed optimal.
- Conservatism of event filtering is not a calibrated error probability;
  complex repeats, clipped boundaries, systematic errors, and mixed alleles
  require biological evaluation.
- No genome-scale runtime, memory, or consensus-QV claims yet. No automatic
  package discovery, pipeline integration, or multi-package CRAM merge.

The broader [design proposal](../CONSENSUS_V2_DESIGN.md) includes future work;
this README describes the implemented prototype.

## Memory behavior

The package prepass reads/discards one native read at a time. Main loading keeps
decoded sequences but not native packed blobs; same-read/same-orientation
placements share sequence storage. Opposite orientations may require two copies.
HPC alignment operates in batches of approximately `--batch-bases` raw bases
(also capped at 4096 placements); one read may exceed that budget. Per-read run
arrays and expanded traces are discarded immediately after evidence accumulation.
Only compact CIGAR traces survive to later stages. Raw trace construction
re-encodes one read per worker, but does not repeat its HPC alignment.

Run histograms use a 12-byte inline bin plus sparse overflow for alternative
lengths. At default support, only the winning tier is retained; with
`--min-run-support > 1`, bins are allocated for every specified tier so fallback
remains exact. For 170 million runs, inline storage is about 2.04 GB per retained
tier, plus sparse overflow. The raw-coordinate prefix vector adds about 1.36 GB.
Decoded reads and polishing candidates remain data-dependent memory costs.

Polishing references existing sequences/compact traces instead of copying them.
Candidates failing the existing size/boundary/homopolymer filters are discarded
before storage; indexed interval endpoints provide exact spanning depths.
Stage messages report process high-water RSS in GiB, not current RSS or a memory
limit. `samtools sort` is a separate process with its own memory consumption.

The chromosome-sized user dataset has not been benchmarked locally. E. coli
polish+CRAM before/after measurements and sequence-equivalence checks are recorded
in `../ecoli-test/memory-results.txt`; small-test savings do not guarantee the
same ratio at chromosome scale.
