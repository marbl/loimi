# hpc-consensus

## License and contributing

This project's original source and documentation are licensed under
[MIT No Attribution (MIT-0)](LICENSE), using the
[SPDX MIT-0 license text](https://spdx.org/licenses/MIT-0.html).
Linked dependencies, including Verkko/Canu and their bundled libraries, retain
their own licenses; this license does not relicense them or replace their
distribution requirements.

See [CREDITS.md](CREDITS.md) for algorithmic influences, upstream developers,
libraries actually used, evaluation tools, and suggested scientific citations.

Repository engineering expectations are in [AGENTS.md](AGENTS.md). C++ formatting
and the module map, invariants and regression workflow are described in
[MAINTAINING.md](MAINTAINING.md). Formatting
is defined by `.clang-format`; use `clang-format -i src/*.cpp src/*.hpp tests/core.cpp`
after code edits, then rebuild and run `make test` as described below. Comments should
explain invariants, ownership, coordinate conventions, and non-obvious decisions,
not simply narrate statements.

Standalone C++ prototype for expanding a fixed homopolymer-compressed (HPC)
template using reads and placements from Verkko `.cnspack` files. It reuses the
native Canu/Verkko package reader; it does not use POA or remap reads genome-wide.

Status: prototype. Recent user-reported human chromosome runs reached about
Merqury Q66 with HP realignment and 500-bp flank filtering, versus Q62.7 without
those options. These are results from the experimental checkout, not a new
full-chromosome benchmark of this port. Accuracy remains dataset-dependent.
The existing Verkko pipeline is not modified.

## Default calling workflow

Expansion-only runs use the 500-bp HP flank filter, then write the expanded
template without raw realignment. Add `--polish-non-hpc` to enable linked-event
polishing, HP realignment, and coverage-supported extension at both ends.
CRAM export remains optional via `--cram FILE`.

The promoted features no longer need experimental flags:

| Feature | Default | Disable |
| --- | --- | --- |
| HP query flank filter | 500 raw bases on each side | `--no-hp-flank-filter` |
| HP local realignment | on during polishing | `--no-hp-realign` |
| Linked-event polishing | on during polishing | `--no-linked-events` |
| End extension | on during polishing | `--no-extend-ends` |

Use `--hp-flank-bases N` to change the flank distance. `--hp-realign` and
`--linked-events` explicitly enable their defaults; the last occurrence of
these enable/disable controls wins. Only the standard option names are
accepted. No corroboration model is used; optional ordinary tier models still apply.

## Build with Make

Requires a C++20 compiler with OpenMP, GNU Make, pkg-config, OpenSSL, libcurl,
liblzma, zlib and bzip2 development libraries. Build Verkko first: this tool
links its native `libverkko.a` and must use a compatible compiler/store ABI.
Python 3 is required for tests; samtools enables CRAM tests and export.

When installed as `verkko/src/hpc-consensus`, the Makefile finds headers in
`../canu/src` and the library in `../../build/lib/libverkko.a` automatically:

```sh
make -j 4
make test
```

The standalone executable is `build/bin/hpc-consensus`; objects are in `build/obj`.
For a separate checkout, supply the Verkko repository root:

```sh
make -j 4 VERKKO_ROOT=/path/to/verkko
make test VERKKO_ROOT=/path/to/verkko
```

`CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, `LDLIBS`, `PKG_CONFIG`, and
`PYTHON` can be overridden. Compiler selection is left to Make and the caller;
there is no compiler discovery. Set `CXX` in the environment or on the command
line when needed, using a C++20 compiler with OpenMP compatible with Verkko.
The standard flag defaults to `-std=c++2a`, the spelling accepted by GCC 9
and newer GCC/Clang. Override `CXX_STANDARD` if needed; the Makefile does not
select or replace your compiler.
`OPENMP_FLAGS` and `OPENMP_LIBS` default to `-fopenmp` for compile and link.
Set `VERKKO_BUILD` for a custom parent build directory, `VERKKO_LIB` for an
explicit archive. Use a fresh `BUILD_DIR` or
`make clean` when changing compiler or flags.

### Parent Verkko integration

The parent build should depend on its completed `libverkko.a`, then invoke
the submodule. An example recipe, to adapt in Verkko's makefiles:

```make
# Dependency: $(TARGET_DIR)/lib/libverkko.a
$(MAKE) -C hpc-consensus CXX="$(CXX)" \
  VERKKO_BUILD="$(abspath $(TARGET_DIR))"
```

The submodule never invokes the parent build, fetches dependencies, or edits
Verkko. This avoids recursive build cycles. The executable remains at the
submodule's `build/bin/hpc-consensus`; the parent controls packaging or copying
it elsewhere. Models remain optional external files.

## Inputs

- `--import`: one native `.cnspack`, potentially containing multiple tigs.
  Reads, stored names, orientations, clipping, and HPC placement intervals come
  from this file. A separate read FASTA or sequence store is not needed at runtime.
- `--template-hpc`: plain FASTA of HPC path sequences. The current package writer
  does not embed these sequences, so this companion file is normally required.
  Embedded tig bases are accepted as an HPC template if present. If both sources
  are supplied, they must agree.
- Template FASTA names must encode the packaged tig ID: numeric IDs,
  `tig00000109`, and legacy `piece000109` names are recognized directly.
  The numeric suffix must match the internal ID; arbitrary names and differently
  numbered pieces are not translated. No tig-name map is required or accepted.
  Automatic scfmap/GFA-path association is not implemented.
  All output references use `tig` plus the internal ID zero-padded to eight digits
  (`tig00000109`), consistently in FASTA, CRAM, VCF, and diagnostics. Larger IDs
  retain all digits rather than being truncated.
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
build/bin/hpc-consensus \
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

## Regional iteration

To process only the first 10 Mb of internal tig 109, add:

```sh
--region tig00000109:0-10000000 --region-flank 100000
```

Coordinates are zero-based, half-open **uncompressed draft** bases, not HPC
runs, template-name suffixes, or mature Verkko consensus coordinates. Polishing
translates these boundaries through accepted edits; final length can therefore
differ from the requested length. An end beyond the contig is clamped; a start
beyond it is an error. Only one region is accepted per invocation.

The implementation counts placements across the entire package, scans again to
retain at most 32 ranked offset anchors, then retains only reads overlapping a
conservative HPC prefix. Since every run expands to at least one base, END plus
the flank bounds the required number of HPC runs. Full overlapping placements
are retained as context. Raw alignment/polishing is restricted to reads overlapping
the requested draft interval plus the flank. The flank is used conservatively in
HPC runs during read selection and in raw bases during polishing selection.

This is not just output truncation, but it is also not random access: all native
read records are scanned three times, and complete template strings and placement
metadata are still loaded. For a nonzero start, run calling must process the
prefix through END to establish raw coordinates. Long boundary-spanning reads can
extend the working prefix. No chromosome-scale memory bound has been measured.

FASTA, CRAM and edit VCF use the canonical tig name with crop-relative coordinates;
CRAM contains overlapping alignments soft-clipped to the final crop, with original
placement multiplicities retained. Edits touching a crop boundary are not applied.
`.region.tsv` records the draft crop, corresponding polished-prefix boundaries,
and processed HPC prefix length. Multi-read `.run-calls.tsv` keeps **prefix** run
and draft coordinates, including context; singleton run calls are crop-relative.
Summary counts describe processed evidence/context, not exclusively emitted bases.
Separate regions of the same tig have the same name: keep their outputs separate.

## Calling rules

### Offset alignment in HPC space

Layout coordinates and template coordinates may have different origins.
Each aligned read estimates their translation in compressed runs; this is not
a raw-base offset or a global parameter shared across tigs.
For a whole-read anchor, `offset = layout_begin - aligned_template_begin`.
The caller then maps a layout coordinate with `template_position =
layout_position - offset`. Internal anchors account for their position within
the read before estimating the same translation.

#### Whole-read discovery and validation

Template offsets are inferred independently for each tig. The tool tries up to
32 uniquely placed, eligible reads, ordered by tier and then proximity to the
layout start, preferring reads away from possible overhangs. Full compressed
reads are aligned to the template; tied optimal placements are rejected.
With six or more eligible reads, six placements are required, split into discovery and validation sets.
Each set must support a consistent offset cluster (at least three reads and
75% of estimates within a 20-HPC-base span); the two estimates must agree
within 20 bases. The combined cluster's lower median is used. This is a
heuristic consistency check, not a calibrated confidence probability.

#### Internal anchors and layout-overlap rescue

If whole-read inference fails, the fallback tests five internal chunks per
candidate read (up to 5,000 HPC bases each). Only chunks occurring exactly once
in the template provide anchors. At least 75% of informative reads must agree
within a 20-base offset span, with at least three supporting reads (one for
sparse layouts). Each read counts once. Five agreeing reads out of six can
therefore qualify; requiring six agreeing reads would defeat the outlier rule.
A distant cluster with at least half as much support is rejected. At least
three supporting reads (one for sparse layouts) must also pass predicted-overlap
alignment, allowing overhangs at template ends. Noisy reads without exact
internal anchors may still fail.

If these unique-read paths fail, a final `LAYOUT_OVERLAP` rescue samples each
user-ordered tier, using a deterministic sample of up to 64 distinct read IDs
per tier. Known non-unique placements are allowed for coordinate inference;
this does not promote their consensus votes or MAPQ. The shorter of the read
and template is aligned inside the longer. Failed partial overlaps can instead
use up to 64 internal chunks (up to 2,000 HPC bases), requiring separated,
collinear matches that occur once in the template. Tied alignments at distant
loci and anchors shorter than 128 HPC bases are rejected. Local optimal
endpoints within 256 HPC bases describe one uncertain locus. Noisy matches
pinned to either sequence boundary are checked with internal anchors instead
of assuming they represent a complete overlap.
If the initial exact internal chain fails, chunk widths are halved down to
128 HPC bases. Each chain still requires separated, unique template matches
and subsequent fixed-overlap validation; shorter chunks only propose origins.

The rescue evaluates shared candidates against all sampled reads, allowing
lower tiers to corroborate a higher-tier origin. Local proposals are grouped
within 256 HPC bases and represented by the per-ID lower median in the highest
available tier. More than 128 candidate loci fails explicitly. Each candidate
is scored using **global alignment of the layout-predicted query interval**;
it cannot relocate to a different repeat copy inside an enclosing read.

#### Resolving competing offset candidates

Each molecule contributes in its highest tier with usable overlap scores.
Repeated IDs and recognized `hifi_UUID` / `ont_UUID` corrected/raw pairs count
once; other names use their packaged IDs, without guessing biological aliases.
The comparison cost is edit fraction plus `0.01 * (1 - overlap / longest_overlap)`
to avoid rewarding tiny perfect overlaps. Candidates within 0.001 of a
molecule's best cost share its support. Each tier is evaluated independently in
user order, with unique support breaking an equal count within that tier.
The first decisive tier wins: its leading candidate must meet the independent
molecule minimum (or unanimous informative support within that tier), have at
least 60% of that tier's informative molecules, and
have no equally supported distant candidate in that tier. A split, tied, or
insufficiently corroborated tier falls through to the next tier. Lower tiers
can corroborate or resolve an inconclusive higher tier, but cannot override a
decisive higher tier. A candidate requires
one molecule for a one-molecule sample, two for samples of two or three, and
three for larger samples (`min(3, n / 2 + 1)`), counting independent support
across tiers. The selected tier is reported in the rescue log.
Unanimity permits a single informative molecule when other sampled reads have
no usable scores, but never bypasses the competing-distant-offset check.
This explicitly relaxes the independent-support minimum, not the alignment gate.
Nearby
representatives within 512 HPC bases are
treated as local coordinate uncertainty, not competing repeat copies.
These are heuristic consistency checks, not calibrated confidence probabilities;
consensus tier voting and MAPQ are unchanged. Region mode rereads the
package with the same bounded ID sample only if its initial anchors fail.
If every tier remains unresolved, a final repeat-specific tie-break verifies
an exact period spanning at least two copies across the **entire** HPC template.
Only equally supported candidates separated by whole repeat units (allowing
at most 20 HPC bases of phase drift, further capped at one quarter-period) can
use this fallback. The smallest valid offset wins that tie. A better-supported
later offset, including one resolved by a later tier, always takes precedence.
Partially repetitive templates and incompatible repeat phases still fail.
Period discovery tests at most 128 prefix matches and uses no chromosome-sized
scratch array. This is a deterministic coordinate convention for equivalent
repeat copies, not evidence that a particular biological copy was identified.

If these checks still fail, the caller rescans the **same candidate origins**
with bounded endpoint freedom. Small layout drift can otherwise penalize the
correct repeat copy simply because predicted clipping endpoints are inexact.
This failure-only pass requires at least 1,000 HPC bases of predicted overlap,
trims up to 256 bases (at most one eighth of the overlap) from each reference
end, and aligns that core inside the original predicted query overlap. Both
aligned endpoints must remain within the same local uncertainty window.
It applies the same molecule counting, tier order, majority, and ambiguity
checks; it does not lower voting thresholds or search the rest of a long read.
Successful rescues are recorded as `LOCAL_LAYOUT` in `template-offsets.tsv`
and `local_drift=1` in the log. They add bounded per-read alignment work only
for otherwise failing tigs, with no persistent per-read sequence copies.

#### Boundary reads and small tigs

Unresolved tigs fail explicitly. Reads crossing the inferred template boundary
are aligned semi-globally to their overlapping HPC template interval. Overhangs
are clipped at aligned run boundaries and represented as soft clips in CRAM;
they do not extend the template during this initial alignment. Eligible
soft-clipped overhangs can subsequently support end extension. Left, right, and two-sided
overhangs are supported. For repeat-rescued origins, query clipping is aligned
near the layout-predicted position with at most 256 HPC bases of adjustment,
strictly less than half the verified period. Optimal endpoints within that
window are resolved by proximity to the prediction, then smaller query start;
an alignment cannot jump to another repeat copy. Local-drift rescues also
constrain query clipping to within 256 HPC bases of the layout prediction. Other offset methods keep
their existing alignment rules. Ambiguous best endpoints or excessive alignment error
produce `HPC_ALIGNMENT_FAILED`; reads wholly outside the template produce
`OUTSIDE_TEMPLATE`. Existing raw clipping is retained. If it cuts through a run
(the adjacent clipped base is identical), that partial run cannot vote for a
run length. The old global `--template-offset` option has been removed.
`PREFIX.template-offsets.tsv` records offsets, attempts, placements, validation
support, tolerance, and method (`VALIDATED`, `SPARSE`, `INTERNAL_OVERLAP`,
`LAYOUT_OVERLAP`, `REPEAT_LAYOUT`, or `SINGLE_READ`). For `LAYOUT_OVERLAP`, placements and
validation support count supporting molecules, and tolerance reports the
256-HPC-base proposal span. For `INTERNAL_OVERLAP`, they count supporting reads
and successful coordinate-constrained overlap checks, not held-out reads.
With fewer than six eligible reads in a multi-read tig, every eligible read is
aligned and the available unique placements are clustered without held-out
validation. At least one placement is required; inconsistent estimates still
fail the 75% cluster agreement check. `SPARSE` explicitly indicates weaker
evidence, not validated confidence.

A tig containing exactly one packaged read bypasses inference and consensus:
its oriented, clipped sequence is output directly, without requiring a template.
It must match an accepted tier and have a known placement count. Run support
thresholds do not apply to passthrough. Optional CRAM preserves the original
oriented read and soft clips, with a perfect match to the emitted sequence.
Polishing leaves singleton sequences unchanged. No offset is inferred (reported
as `.`), and run-call coordinates describe the passthrough sequence.

### HPC run-length calling

#### Tier selection and evidence

Placements with packaged placement count exactly one and no distinct alternative
in the input package are unique. Run calling searches all unique tiers first,
requiring at least five observations by default. If none qualifies, it searches
the tiers again, pooling unique and known ambiguous observations **within each
tier**, never across technologies. Thus qualifying unique UL beats pooled LA.
If no pooled tier qualifies, the first tier with any evidence supplies a sparse
fallback. `--unique-only` excludes ambiguous votes, retaining the same qualifying-
then-sparse selection among unique tiers. Optional non-HPC polishing retains its
separate unique-then-ambiguous groups and polishing thresholds; its policy is not
changed by the run-calling threshold. `--unique-only` also filters polishing for
controlled comparisons; `--ambiguous-fallback` explicitly selects the default.
Unknown placement counts (zero) remain export-only. Entirely non-unique
multi-read tigs can use the layout-overlap rescue when their offsets are
unambiguous. `--unique-only` still excludes these reads from consensus votes,
but not from offset inference.
All distinct eligible placements are exported when CRAM is requested.
The legacy `skipConsensus` flag is ignored:
explicit user tiers replace its template-building priority policy. Any number
of ordered tiers is supported. All eligible reads can be aligned and exported
to CRAM, even if their tier never wins a consensus call. Repeated IDs are allowed
within and across tigs. Identical entries (ID, coordinates, strand, and clipping
within a tig) are deduplicated. A package prepass counts distinct occurrences.
Uniqueness is trusted from the layout, not independently inferred by mapping.

#### Run-length estimator

For each template run, only matching aligned read bases contribute a run length.
The first qualifying tier in the search above wins; lower tiers cannot outvote
it. `--min-run-support` is a preference threshold, not a hard evidence discard:
sparse fallback prevents covered runs becoming unsupported. Setting it to one
restores preference for even one unique observation. The called length is the
empirical histogram mode. Tied modes
are resolved by minimum total absolute distance to observations, then smaller
length. Without a tier model this is an empirical vote, not a fitted predictor. Quality priority is
entirely user-specified by tiers, not per-base quality scores.

`.run-calls.tsv` reports `placement_class` as UNIQUE, AMBIGUOUS, POOLED (both
placement classes), or NONE. Summary JSON counts `ambiguous_run_calls`,
including POOLED calls. Depth and peak refer to the selected histogram,
including both classes when pooled; peak is the modal observation count.

#### Regional fixtures

For native regional fixtures, the build also includes `extract-regions`:

```sh
build/bin/extract-regions input.cnspack 109 3424 regions.tsv
```

Each regions.tsv row is `output_prefix HPC_begin HPC_end` in original-template
zero-based half-open coordinates. The offset must be independently verified.
The utility extends each window to retain full overlapping reads, copies native
read records verbatim, rebases layout positions, and preserves original ambiguity
(the native placement-count field saturates at 31). Output `.span.tsv` gives the
template slice to supply; `.placements.tsv` retains original positions for
coverage diagnostics. It does not automatically extract template/truth FASTAs.

Without eligible evidence a run gets length one and is reported as
`NO_EVIDENCE`; `--missing-run error` instead stops the run. N runs do not collect
run-length evidence. Fast-mode output preserves the input HPC base sequence.

### Raw-base polishing and linked events

Enable this stage with `--polish-non-hpc`. Reads are aligned to the expanded
draft. By default, changes on the same read are grouped when left-normalizing
one reaches an earlier event. The combined allele retains intervening matching
bases; nearby unrelated variants are not automatically merged. This gives
equivalent split and contiguous CIGAR representations a shared vote without
assembling an allele from different reads. `--no-linked-events` restores
independent normalization for comparison.

Optional correction calls substitutions and short non-HPC events from raw traces.
It selects the first tier with sufficient spanning depth and requires alternate
support and fraction thresholds. If both strands span the event in that tier,
both must support it. Terminal events, pure run-length changes, overlapping or
adjacent accepted edits, and alternate alleles containing N are excluded.
There is one correction pass, not iterative polishing.

### Homopolymer flank filtering and local re-estimation

The flank filter requires 500 raw query bases and an immediately adjacent
matching HPC run on each side of an interior observation. It does not require
500 perfectly matching bases. Usable overhang clipped only by a true template
boundary counts toward the distance, so the filter preserves end-extension
support; package trimming does not count. Template endpoints waive only the
outward-facing requirement. Singleton passthrough remains unchanged.

HP realignment runs by default during polishing. It rebuilds evidence for uncertain or
model-adjusted runs using 64-base local alignment context and 12 exact terminal
matches, then reapplies the configured tier models. Unique tiers must qualify
independently before within-tier ambiguous fallback is considered. Rebuilt
observations also obey the requested query flank distance. Targets lacking
context, runs longer than 1,000 bases, and query windows longer than 4,096 bases
are skipped. Ordinary polishing edits take priority over overlapping HP edits.

Accepted HP edits are recorded in `.hp-realignment.tsv` in zero-based, half-open
raw draft coordinates, including any left extension. The existing
`.non-hpc-edits.vcf` includes all accepted edits for coordinate lifting, despite
its historical name; the summary's `non_hpc_edits` count excludes HP edits.
Final CRAM alignments, when requested, use the corrected final sequence.
Polishing enables end extension unless `--no-extend-ends` is supplied.

The implementation is isolated in `src/hp_refinement.*` and `src/hp_flanks.hpp`.
It retains sparse candidate histograms and temporary per-read coordinate maps,
not chromosome-length arrays for every read. Local realignment is currently
serial and adds runtime proportional to read overlap with candidate sites.
Tests cover model uncertainty, unique/ambiguous fallback, flank boundaries,
actual post-draft edits, CRAM validity, and extension on both ends.

### Coverage-supported end recovery

Polishing now first recovers read overhangs clipped solely because the HPC
template ended. This works at both tig ends; `--extend-ends` enables the same
step without internal polishing, and `--no-extend-ends` preserves the previous
boundary behavior. Fast uncompression without either option is unchanged.
Packaged raw clipping is respected, and singleton passthrough is unchanged.

Only reads anchored to the existing terminal sequence can establish an extension
(at least 20 aligned operations, with the terminal 100-operation error fraction
at most `--max-error`). A molecule counts once per end, even with duplicate
placements; unknown placement counts cannot vote. A high-priority group with
at least `--min-end-support` molecules supplies a seed path. Up to three bounded
rounds of anchored HPC alignment correct that path's substitutions, insertions,
and deletions. Each path position uses the first tier with sufficient depth,
trying unique evidence before same-tier pooling; a strict majority chooses the
allele, and an unresolved tie retains the seed allele. A qualifying higher-tier
conflict cannot be overridden by lower-tier evidence. Partial terminal runs do
not vote. Repeated paths terminate refinement rather than oscillating.

Tied prefix endpoints are resolved to the shortest optimal endpoint locally:
the read is already anchored to this contig end, so these are not alternative
genomic placements. Global template-offset uniqueness checks are unchanged.
Reads are realigned to the corrected HPC path to collect run lengths, using
unique-first then same-tier pooled selection and the supplied calibration
models. Unlike interior calling, every
added run must reach the end-support threshold: there is no sparse fallback.
The last run of each tail is excluded because it may be truncated. Recovery
stops at the first unsupported corrected-path run or the raw-base cap.

This is conservative seeded extension, not de novo assembly of arbitrary tails:
it may still stop early where path alignment or read support is insufficient,
and does not establish a telomere's exact biological repeat count. Read
overhangs, including chimeric tails, still require inspection in difficult cases.
The cap bounds terminal-tail scratch, not total process RSS. Recovered tails
are aligned to the added sequence and remaining/unresolved bases stay soft-clipped.
Left extension shifts all retained alignment coordinates; CRAM NM/MD and contig
lengths use the extended, polished sequence. Neither extension changes the
existing interior sequence itself; ordinary polishing can still edit it.

`.end-extensions.tsv` records the original draft length, bases added on each
side, and the number of anchored candidate molecules; summary JSON includes
`left_extension_bases` and `right_extension_bases`.
`.run-calls.tsv` remains in the **original, unextended template/draft coordinates**:
add `left_added` to its raw coordinates before interpreting them in the extended
draft. `.draft.fasta` and polishing VCF coordinates refer to the extended draft.
Extension is separate from VCF edits. `--region` never extends an internal crop
boundary; only a crop reaching a true template end can include an extension.
For manually extracted fixtures, use `--end-extension-side right` or `left` to
avoid treating the fixture's artificial opposite edge as a biological end.

## Options and defaults

| Option | Default | Meaning |
| --- | --- | --- |
| `--threads N` | 1 | Parallel per-read alignment work within a tig; also passed to samtools |
| `--batch-bases N` | 8000000 | Approximate raw-base budget for each HPC alignment batch; a single read may exceed it |
| `--min-run-support N` | 5 | Preferred tier depth; retain sparse evidence if none qualifies |
| `--extend-ends` | with polishing | Recover coverage-supported template overhangs at both ends |
| `--no-extend-ends` | off | Disable end recovery, including during polishing |
| `--min-end-support N` | 5 | Distinct molecules required at every added run |
| `--max-end-extension N` | 10000 | Raw-base cap per end, and per-read tail bound |
| `--end-extension-side both\|left\|right` | both | Restrict recovery, useful for extracted fixtures |
| `--missing-run one\|error` | one | Policy when no tier qualifies |
| `--unmatched error\|drop` | error | Policy for eligible reads matching no tier |
| `--max-error F` | 0.15 | Unit edit-distance cap, rounded up, relative to the longer sequence |
| `--polish-non-hpc` | off | Enable raw-base correction, HP refinement and end recovery |
| `--hp-flank-bases N` | 500 | Raw query context required for HP observations |
| `--no-hp-flank-filter` | off | Disable that flank requirement |
| `--no-hp-realign` | off | Disable default HP refinement during polishing |
| `--no-linked-events` | off | Disable default linked-event grouping |
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

The repository-local [design proposal](CONSENSUS_V2_DESIGN.md) includes current
implementation status, historical rationale, and future work;
this README describes the implemented prototype.

## Optional per-tier run-length calibration

Attach a separate model to any named tier; tiers without a model keep the mode
caller. Models do not alter tier order, unique-placement preference, eligibility,
or singleton pass-through. HP realignment reapplies the selected tier model to
its rebuilt evidence; ordinary non-HPC polishing does not use this model.

This is a **discrete conditional-probability model**, not linear regression,
nonlinear regression, or a neural network. For each base and possible true run
length, a table stores `P(observed_length | true_length, base, tier)`. Given a
histogram with counts `n[o]`, the caller selects the true length `L` maximizing
`sum_o n[o] * log P(o | L, base, tier)`: maximum-likelihood inference with no
length prior. Reads are treated as independent observations for scoring.

The supplied calibration script estimates these categorical probabilities from
truth-labeled observations, normalizing each locus to equal weight. It smooths
toward the A/T or C/G base-family distribution and adds a small uniform outlier
component, so unseen observations do not have zero probability. The model file
contains the final probability tables; no fitting occurs during consensus.

```bash
--tier 'LA=hifi_m*' --tier 'LAO=hifi_[0-9a-ln-z]*' --tier 'UL=ont_*' \
--tier-model LA=LA.run-model.tsv \
--tier-model LAO=LAO.run-model.tsv \
--tier-model UL=UL.run-model.tsv
```

These are **dataset/chemistry-specific experimental models**, not universal
technology defaults. Validate a model after changing the read-generation or
correction pipeline. The repository includes `models/revio.run-model.tsv`; it is used only when
explicitly selected with `--tier-model`, not automatically. Mode remains the fallback
below the model's minimum depth, for an unmodeled base, or when the original
modal length is outside that base's calibrated range. The summary reports `model_changed_runs`, and
`.models.tsv` records tier-to-file associations. Run-call `peak` remains the
modal observation count, **not** a calibrated confidence score.

Model files are whitespace-separated text, with optional `#` comment lines:

```text
HPC_RUN_MODEL_V1 OBSERVED_MAX MINIMUM_DEPTH
BASE TRUE_LENGTH TRAINING_LOCI P_OBS_0 P_OBS_1 ... P_OBS_MAX
```

For each A/C/G/T base present, rows must enumerate contiguous true lengths from
one, in ascending order. Each row has `OBSERVED_MAX+1` finite positive
probabilities summing to one. The final observed bin includes all lengths at
least `OBSERVED_MAX`. Prediction maximizes the sum of histogram counts times
log conditional probabilities; ties favor the shorter length. There is no
implicit length prior. Correlated reads require separate confidence calibration.
Keep training provenance and the model checksum alongside each model file.

For calibration, `--evidence-sites FILE` accepts `TIG RUN` rows (canonical output
tig names, zero-based **input-template HPC** indexes), exporting
`.run-evidence.tsv`: `tig,run,base,tier,placement_class,length,count` as TSV.
This observes accepted, non-truncated HPC votes **before tier prioritization**,
including lower tiers. Unique and ambiguous votes remain separate; unknown
placement counts remain ineligible. `--unique-only` also excludes ambiguous
votes from this export. `--evidence-exclude-reads FILE` takes one numeric read ID
per whitespace-delimited token and excludes them from the diagnostic only,
useful for preventing training/validation read overlap. Neither switch changes
the consensus. Singleton pass-through produces no histogram evidence. Use
native extracted fixtures rather than `--region` with this export.

Evidence export adds one byte per template run plus sparse histograms only at
requested sites. Sample sites for calibration rather than exporting every run
of a whole chromosome. Model storage is shared across runs; homogeneous
histograms use precomputed single-allele predictions.

`scripts/calibrate_runs.py` is an experimental LA/LAO/UL trainer/evaluator:

```bash
python3 scripts/calibrate_runs.py --data /path/to/calibration-data train
python3 scripts/calibrate_runs.py --data /path/to/calibration-data validate
```

The dataset's `manifest.json` contains `windows` entries with `name` and `split`
(`train` or `validation`). Each named subdirectory contains `labels.tsv`
(`run`, `base`, `truth_length`), the native `.run-evidence.tsv` export named
`consensus.run-evidence.tsv`, `consensus.summary.json` as the completion marker,
and `label-summary.json` describing labeling coverage/exclusions. Labels must
come from trusted, unambiguous truth alignments in the same HPC coordinate
system. Prevent read sharing across the split before exporting evidence.
Exported sites without an accepted truth label are ignored, allowing stricter
truth masks without rerunning the read alignments.

The trainer uses unique placements, one unit of histogram mass per locus,
base-family smoothing, and a small outlier component. It requires 30 labeled
loci per base/true length and stores only contiguous supported ranges from one.
The generated models require at least five observations for prediction. These
are explicit experimental choices recorded in `models/metadata.json`, not
universal recommended settings. The current sampling contract is all labeled
truth runs of length >=10 plus 1/100 shorter runs; weighted summaries are valid
only for that design and do not account for excluded/unmappable truth regions.

Training reads only `train` windows and refuses to overwrite frozen models.
Validation checks model checksums and split names, writes `validation.json`,
and never refits. It reports individual tiers, ambiguity classes, depth/length
strata, and LA-first/LAO-first priority simulations. These are sampled run
errors, not full-sequence QV. Keep an additional independent dataset for final
assessment when comparing or tuning models against these validation results.

For predictor parity checks, `build/bin/apply-run-model MODEL` reads stdin rows
`BASE LENGTH:COUNT [LENGTH:COUNT ...]` and prints the mode and calibrated call.
This diagnostic uses the production C++ predictor without rerunning alignments.

## Memory behavior

The package prepass reads/discards one native read at a time. Main loading keeps
decoded sequences but not native packed blobs; same-read/same-orientation
placements share sequence storage. Opposite orientations may require two copies.
HPC alignment operates in batches of approximately `--batch-bases` raw bases
(also capped at 4096 placements); one read may exceed that budget. Per-read run
arrays and expanded traces are discarded immediately after evidence accumulation.
Only compact CIGAR traces survive to later stages. Raw trace construction
re-encodes one read per worker, but does not repeat its HPC alignment.

Run calling processes one tier/placement class at a time without reordering CRAM
records. Dense accumulation uses an 8-byte bin, a 16-byte selected call and two
one-byte flags per run, independent of tier count: about 4.42 GB for 170 million
runs, plus sparse overflow. Below-threshold unique histograms are retained for
later same-tier pooling as 16-byte records per distinct observed length. This
storage is data-dependent and can be substantial: with threshold five, at most
four records per run per tier are needed, plus vector spare capacity. Raising
the threshold increases this bound. This is not a process memory cap.
Accumulators and pending evidence are released before allocating the
raw-coordinate prefix vector (about 1.36 GB for 170 million runs).
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

## Local HPC vote filtering

Run-length voting now rejects locally poor parts of otherwise acceptable HPC
alignments. By default, each matching operation must have at most 15% X/I/D
operations in a centered window of radius 500 (up to 1,001 HPC alignment
operations). Terminal windows are truncated, and insertions count in both the
error numerator and operation denominator. This is separate from `--max-error`,
which controls whole-read alignment acceptance.

Use `--hpc-filter-radius N` and `--hpc-filter-error F` to adjust the gate, or
`--no-hpc-filter` to reproduce the prior voting behavior. These settings apply
to all tiers; they do not change priority, multiplicity, model parameters, or
mode tie breaking. The filter is not the experimental ambiguous-read pooling
policy. Singleton pass-through remains unchanged.

Rejected positions contribute neither consensus votes nor diagnostic evidence
histograms. Summary fields `filtered_hpc_votes` and `locally_filtered_reads`
count rejected otherwise-eligible matching operations and affected read
placements. If every vote is rejected, the existing lower-tier fallback and
missing-evidence policy apply. Reads remain available for optional raw alignment,
polishing and CRAM export; the filter does not clip or remove CRAM records.

The rolling calculation is linear in trace length and uses one temporary byte
per operation for one read during voting, not chromosome-wide per-read masks.
Exact traces skip allocation. Clean portions of damaged reads remain eligible.
This addresses locally non-collinear tails that pass the whole-read cutoff;
the defaults are not a guarantee of optimal filtering for every technology.

In the HG002 35-Mb diagnostic fixture, filtering reduced the flanked hotspot
from 4,182 incorrect runs (6,424 differing bases) to two runs (three bases),
with no changed run lengths elsewhere in the ~11-Mb fixture. The remaining
errors include an existing long-run disagreement and a terminal-read tie
exposed by removing bad evidence. This is a regional validation, not a new
whole-chromosome QV measurement or cross-technology benchmark.
