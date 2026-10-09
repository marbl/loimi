# Design status and implementation guide

Updated: 2026-10-09

This file includes the original 2026-10-02 proposal below for rationale and
research context. It is historical: where it conflicts with this status section,
the [README](README.md), CLI help, and tested implementation describe current
behavior. Proposed options and guarantees below are not automatically supported.

## Current implementation versus original proposal

- Optional `--tier-model LABEL=FILE` associates calibrated conditional run-length
  likelihoods with each technology tier. Mode remains the default/fallback;
  models do not change tier eligibility. Sampled all-tier HPC evidence export
  supports separate training and held-out validation. Model formats and limits
  are documented in README; models are not universal chemistry defaults.
- Consensus now prefers all unique tiers before falling back through non-unique
  tiers, rather than leaving ambiguous-only regions at run length one.
  `--unique-only` retains the old policy for comparisons. Offset inference still
  requires unique anchors. The run-length mode remains the default; a rounded
  mean option is explicitly experimental and not a validated replacement.
- Regional iteration supports one `--region tigID:BEGIN-END` in uncompressed
  draft coordinates, with context and cropped FASTA/CRAM/VCF output. It retains
  prefix-overlapping reads after separate offset discovery, while preserving
  full-package multiplicity. It still scans the package sequentially and calls
  the prefix through END; see README for coordinate and memory limitations.
- C++20, native Verkko package readers, and unit-cost edlib alignment are
  implemented. Affine-gap ksw2/WFA backends and direct reuse of packaged CIGARs
  remain proposals, not implemented paths.
- `--tig-map` maps `internal_ID FASTA_name` rows to template sequences; without
  it, numeric IDs, `tig00000001`, and legacy `piece000001` names are accepted.
  Outputs always use `tig` plus the internal ID zero-padded to eight digits.
  Automatic scfmap/GFA-path association and modified package writers are pending.
- Offsets are inferred independently per tig. A single eligible unique read
  passes through in layout orientation with clipping; sparse multi-read tigs
  align all eligible reads without held-out validation. Larger sets use
  discovery/validation clusters. Confidence is heuristic.
- Boundary-crossing reads are retained with semi-global HPC overlap alignment
  and raw soft clips. Partially clipped runs cannot donate truncated lengths.
- The legacy `skipConsensus` flag is ignored. Explicit ordered wildcard tiers
  determine priority; all unique tiers precede the non-unique fallback group.
- Distinct repeated placements are allowed in CRAM. The first successfully
  emitted occurrence is primary; later ones are secondary with equal
  multiplicity-adjusted MAPQ. This supersedes the original single-record policy.
  Identical entries are deduplicated. Primary bookkeeping is per invocation.
- Run lengths use the empirical mode, then absolute-deviation and shorter-length
  tie breaks. Per-base quality weighting and sidecar lengths are not implemented.
- Correction and CRAM are independent opt-ins. Correction defaults to at least
  five supporting reads, 90% allele fraction, and a 500-base maximum event.
  It excludes pure run-length changes and conflicting/adjacent edits.
  Additional stability/margin checks described below remain proposals.
- Raw repair currently realigns the complete affected read interval, not only
  small local blocks. samtools performs CRAM sort/index/validation; the module
  does not directly implement an htslib writer.
- FASTA, TSV, JSON, and optional draft FASTA/VCF diagnostics are plain files.
  Compressed/indexed VCF and the proposed extra BED reports are not implemented.
- Alignment scratch is batched, sequences are shared, traces are compact,
  histograms use inline bins plus sparse overflow, and polishing borrows
  alignments with indexed coverage queries. This is not the proposed fully
  windowed, disk-backed, hard-bounded-memory design: reads, histograms, and
  candidates still scale with tig/input size.
- Synthetic regressions and E. coli comparisons exist. They are not a release
  qualification panel. There is no demonstrated chromosome-scale memory ceiling,
  global-QV equivalence, or basis yet to replace the mature Verkko caller by default.

## Source layout

- `src/main.cpp`: CLI, native package loading, placement bookkeeping, per-tig
  inference, batched evidence accumulation, staging, and optional CRAM export.
- `src/core.hpp` / `src/core.cpp`: RLE/coordinate contracts, pairwise alignment,
  modal calling, conservative edits, coordinate lifting, and NM/MD validation.
- `tests/core.cpp`: algorithmic and coordinate regressions.
- `tests/integration.py`: native-package and optional CRAM end-to-end tests.
- `tests/prepare_ecoli.py`: revision-specific historical fixture preparation,
  not a general graph-to-template converter.

Follow [AGENTS.md](AGENTS.md) for maintenance expectations. The repository's
original code and documentation are [MIT-0](LICENSE); dependencies retain their
own licenses.

---

## Historical proposal (2026-10-02)

The following text preserves the original proposal, including superseded design
choices, unimplemented features, and anticipated performance claims. Use the
status above and README for the current contract.

# Standalone HPC consensus/uncompression tool: design proposal

Status: revised design, 2026-10-02

Working name: `hpc-consensus` (the name is deliberately not Verkko-specific)

## Executive decision

Build a standalone **HPC decompressor**, not a general consensus/POA engine.

The graph has already decided the sequence of homopolymer-compressed bases. The tool's primary job is therefore to infer one positive run length for every base in that fixed HPC sequence and expand it. Read realignment is optional.

The central invariant of the uncompression stage is:

```text
homopolymer_compress(uncompressed_draft) == input_hpc_template
```

The tool may optionally make a small, audited set of high-confidence non-homopolymer corrections and may optionally write final read alignments. Both are separate opt-in stages: the fast default stops after run-length expansion and writes only the uncompressed consensus and run-call diagnostics.

The core v1 path will:

1. read the existing Verkko `.cnspack` container, including HPC layout coordinates, uncompressed read sequences, and original read names;
2. accept at most one placement per read and discard ambiguous reads before pileup construction;
3. assign reads to user-ordered quality tiers using explicit command-line name wildcards;
4. align/project reads in HPC space and collect observed lengths for individual template runs;
5. call each run from the **peak (mode) of its observed run-length distribution**, not a basewise majority, mean, or median;
6. expand the fixed HPC template;
7. write the uncompressed FASTA and stop.

Two independent opt-ins extend that result:

- `--polish-non-hpc` uses base-space alignments to apply one conservative, non-iterative correction pass and emits an edit VCF;
- `--cram FILE` emits one final primary alignment per accepted read in a coordinate-sorted, indexed, reference-embedded CRAM.

This is substantially simpler and faster than contig-wide template stitching followed by iterative POA. It also makes the biological decisions explicit: the graph owns base order; uniquely placed reads own run lengths.

## Why POA is not the main algorithm

Recent POA work has made sequence-to-graph alignment markedly faster:

- [POASTA](https://academic.oup.com/bioinformatics/article/41/1/btae757/7942505) (2025) uses A* search, exact-match extension, and graph-topology pruning for exact affine-gap POA.
- [minipoa](https://genome.cshlp.org/content/36/9/1849) (2026) uses seed-chain-align heuristics, adaptive/static banding, and SIMD; it reports up to a fivefold speedup over abPOA and much lower memory use.
- [Theseus](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag561/8742314) (July 2026) uses sparse wavefronts for exact affine-gap alignment to arbitrary directed graphs and reports 1.6-17.6x speedups over POASTA in its evaluated MSA problems.

Those algorithms solve a harder problem than this module has. They are useful when the consensus base path is unknown and multiple sequence alternatives must be represented. Here the assembly graph has already supplied that path. Rebuilding a mutable POA would spend time rediscovering base order, introduce progressive-order effects, and create a second authority that can disagree with the assembly graph.

The relevant precedents point toward direct run-space processing instead:

- [Shasta](https://www.nature.com/articles/s41587-020-0503-6) assembles in run-length space and, in its raw assembly consensus, chooses the repeat count observed most frequently among contributing reads. This directly supports a modal run-length caller.
- [MBG](https://github.com/maickrau/MBG) already emits read paths that are exact in HPC space and may differ only in homopolymer lengths. Its `HPCConsensus` code projects read positions directly to unitig positions and its `ConsensusMaker` accumulates expansion counts.
- [LJA](https://github.com/AntonBankevich/LJA) uncompresses an HPC assembly by compressing reads, using layout-restricted banded alignment, and collecting per-position run lengths. LJA uses a median for ordinary homopolymers and local SPOA for selected complex repeats; the proposed tool adopts the efficient fixed-backbone structure but uses the required distribution peak and does not make local POA part of v1.
- Current [Verkko](https://github.com/marbl/verkko) already provides an HPC graph/path representation and final CRAM packaging. The reusable ideas are graph-to-contig projection, partitioning, and htslib output—not Canu's read-tiling template construction.

## Scope

### In scope for v1

- A fixed A/C/G/T/N HPC template built from the assembly graph.
- Existing Verkko `.cnspack` input containing `tgTig` layouts and packaged raw reads.
- User-ordered tier wildcards matched against original packaged read names.
- One accepted layout placement per read, as recorded in the package.
- Direct run projection when the package contains an HPC CIGAR or the read/template slices match exactly.
- Fast banded pairwise HPC alignment when the layout gives the exact interval but not a per-run CIGAR.
- Per-run modal length calls with deterministic tie breaking.
- Optional conservative pileup-based correction of substitutions and short non-homopolymer indels.
- Optional final CRAM/CRAI.
- FASTA, run-level confidence, and read-decision reports.
- A tier-aware CLI whose only read-name semantics are the explicit wildcard rules supplied by the user.

### Explicitly out of scope for v1

- General variant calling, haplotype discovery, or graph-branch reconstruction.
- Iterative polishing or corrections without an auditable evidence record.
- General POA/MSA.
- Genome-wide remapping of reads after consensus.
- Retaining alternative placements in the final alignment file.
- Learned or platform-calibrated polishing.

The fixed HPC skeleton remains authoritative during uncompression. The later correction pass may change it only at an explicit high-confidence edit recorded in the output VCF. Complex or competing changes remain graph/path correction problems.

## Input contract

### 1. Existing `.cnspack` container

`--import 7-consensus/packages/part001.cnspack`

The tool consumes the same sequential binary package currently passed to `utgcns`. Each package entry contains:

- a serialized `tgTig` with tig ID, HPC layout length, child read IDs, orientation, compressed `min/max` placement coordinates, clipping, `skipConsensus`, and `_numPlacement`;
- one redundant/template read record required by the existing package format; and
- one `sqRead` record per child, containing the uncompressed sequence and its stored original name.

The reader ignores the redundant first read after validating it, de-duplicates packaged reads by read ID, and retains the original `sqRead_name()` for tier matching and output. The `tgTig` placement is authoritative: the tool never searches for another contig or locus.

Reads with `_numPlacement != 1`, `skipConsensus`, or another upstream ambiguous/discard marker do not contribute. They are omitted from optional CRAM output and reported using stable reason codes. The upstream layout remains responsible for selecting or discarding placements; this tool does not repeat that decision.

`.cnspack` is a native struct-based format rather than a portable, self-describing interchange. The standalone executable must therefore build against the same in-tree Verkko/Canu store compatibility code used to create the package and reject unsupported record layouts cleanly. No conversion to PAF/GAF, FASTA/FASTQ, or a separate read manifest is part of the normal path.

### 2. Graph-derived HPC template

The current `tgTig` serialization has fields for consensus bases, but Verkko's current `layoutToPackage` path leaves `_basesLen == 0`: it packages layout length and children, not the graph-derived HPC base string. Run-length inference still requires those bases.

The preferred integration is to populate the existing `tgTig::_bases` field with the graph-derived HPC template when building the package. This changes neither the `.cnspack` container nor its binary record layout and gives the new tool a single-file input.

For backward compatibility with already-created packages whose tig bases are empty, accept:

```text
--template-hpc assembly.hpc.fasta
```

The template records are keyed to package tig IDs (or resolved through the existing `tigName_to_ID.map`). If both embedded and external templates exist, require exact agreement. Validate that adjacent symbols differ, apart from the documented `N` policy, and that template length equals the packaged HPC layout length.

An optional run-length sidecar can supply fallback lengths for uncovered runs. Without it, the canonical inverse is length 1 and the call is marked `NO_EVIDENCE`.

### 3. Ordered tiers from explicit name wildcards

Read type is not hard-coded or inferred automatically. The user supplies repeated tier rules in descending quality order:

```text
--tier 'hifi=hifi_*'
--tier 'duplex=duplex_*'
--tier 'ont=ont_*'
--tier 'fallback=*'
```

Patterns use documented shell-style glob semantics and must normally be quoted so the shell does not expand them. They are matched against the full original packaged read name. The first matching rule wins; its order is the tier priority. `LABEL=GLOB` supplies a diagnostic label; a bare `GLOB` is also valid and receives the automatic label `tier0`, `tier1`, and so on.

Unmatched reads are an error by default so a spelling mistake cannot silently discard coverage. `--unmatched drop` explicitly discards them. Report the count of unmatched reads and reads matching multiple rules; multiple matches are permitted but resolved by the first-match rule.

The internal read key remains the package read ID plus original name. No built-in rule assigns meaning to prefixes such as `hifi_` or `ont_`; those prefixes matter only when the user explicitly supplies matching patterns.

## Evidence priority: tier order first, uniqueness required

All observations entering the caller come from uniquely placed reads. Ambiguous reads have already been removed.

For template run `i`, let `H[t][k]` be the number of eligible observations of length `k` from tier `t`. The caller scans tiers in user-specified order and selects the first tier with sufficient eligible observations at that run:

```text
winning_tier(i) = first t for which sum_k H[t][k] >= min_run_support
```

The default `min_run_support` is 1, matching the stated policy literally: use the highest-quality tier whenever it has a unique read covering the run; otherwise fall back to the next tier, and so on. The threshold is configurable for users who prefer a no-call over a single-read call.

Only the winning tier votes on the run length. Lower-priority tiers cannot outvote it, regardless of their depth. They may be retained in diagnostics as concordance evidence. If the winning tier has no eligible observation for the next run, selection starts again from the highest-priority tier; fallback is per run, not per contig or window.

This is not a composite numerical score. It is a lexicographic policy with an easily testable outcome.

## Algorithm

### Phase A: encode once in run space

For each oriented read, construct:

```text
run_bases[]       // one A/C/G/T/N symbol per run
run_lengths[]     // positive integer, no 8-bit/255 cap
raw_offsets[]     // prefix sum into the original read
run_qualities[]   // optional aggregate from FASTQ qualities
```

Packaged qualities, when present, are used for input filtering and tie breaking, not as a way for a lower-priority tier to overcome a higher-priority tier.

### Phase B: obtain the read-to-template run trace

There are two paths, chosen automatically.

#### Fast path: direct packaged-layout projection

Use the packaged child orientation, clipping, and HPC `min/max` coordinates to select the exact read and template slices. If an HPC child CIGAR/delta is present, use it directly. Otherwise, when the compressed read slice is identical to the template slice, correspondence is the trivial one-to-one offset map. Neither case needs dynamic programming.

The projection output is a compact run-space trace:

```text
read_run  template_run  op(=,X,I,D)
```

#### Fallback: recover correspondence inside the layout interval

If the package supplies the compressed read placement but the two HPC slices are not identical and no usable child trace is present:

1. extract exactly the oriented read slice and target HPC interval named by the layout;
2. align those two fixed slices with an affine-gap, end-anchored, banded pairwise aligner;
3. begin with the diagonal implied by their compressed coordinates and widen the band geometrically only if the optimum touches a band edge;
4. stop at a configurable maximum and reject the read if correspondence cannot be recovered reliably.

`ksw2` is the pragmatic first backend: LJA demonstrates this pattern, it is small, fast, SIMD-enabled, and suitable for a C++ standalone executable. A WFA2 backend can be benchmarked later. This operation does not discover or score alternative placements; it only recovers indel/mismatch correspondence within the placement already fixed by the packaged layout. It remains pairwise sequence alignment to a fixed string, not sequence-to-POA alignment.

An observation contributes to a template run only when the HPC bases match in a one-to-one aligned column. HPC substitutions and indels are logged but do not mutate the skeleton and do not donate a run length to a neighboring run. An optional minimum exact-match flank can make this eligibility filter more conservative without adding special-case control flow.

### Phase C: call run lengths by distribution peak

For each template run, build a sparse histogram for each tier:

```text
H[t][k] += 1  for every eligible uniquely placed read run of observed length k
```

Choose the first tier satisfying `min_run_support`, then call:

```text
called_length = argmax_k H[winning_tier][k]
```

This is the unweighted empirical mode: the peak of the run-length distribution. It is deliberately not:

- a per-base majority vote, which biases long calls because one read contributes many repeated bases;
- a mean, which can create a length no read supports and is sensitive to outliers;
- a median, which ignores whether the distribution has a sharp peak; or
- a pooled vote across technologies, which allows depth to defeat the user's quality order.

Deterministic ties within the winning tier are resolved by:

1. greater summed, capped run quality when FASTQ qualities exist;
2. smaller total absolute deviation from all eligible lengths in that tier;
3. the optional template-sidecar length; then
4. the shorter length.

Quality can therefore break an exact count tie, but cannot move a unique histogram peak.

For every run, record:

- winning tier and its priority rank;
- observation depth;
- called length and peak count;
- runner-up length/count;
- peak fraction `peak_count / depth`;
- strand counts;
- fallback/no-evidence status; and
- counts of rejected observations by reason.

### Phase D: expand without window reconciliation

Emit each HPC symbol exactly `called_length` times. A prefix sum of called lengths maps every HPC run boundary to a final base-space coordinate.

Processing is partitioned into large run-space chunks. Chunks use alignment halos but have disjoint owned runs, so no sequence overlap or POA-window reconciliation is needed. Memory is bounded by histograms for the owned chunk plus active read traces.

When polishing or CRAM is requested, the tool caches each accepted read's compact HPC trace in a temporary target-sorted stream. The fast/default mode consumes traces into run histograms and discards them immediately, avoiding the extra temporary I/O.

## Optional post-uncompression stages

The fast path ends after Phase D. No base-space alignment is performed and no CRAM is created unless at least one post-uncompression feature is requested.

The four supported combinations are:

| Non-HPC polish | CRAM | Work after expansion |
|---|---|---|
| off | off | None; write FASTA and diagnostics |
| on | off | Build provisional alignments, call/apply edits, write FASTA and VCF |
| off | on | Render alignments once against the uncompressed draft and write CRAM |
| on | on | Build provisional alignments, call/apply edits, then re-render affected alignments and write CRAM |

### Phase E: optionally build base-space alignments and collect residual evidence

The cached HPC trace and the final run-length prefix sums already determine the read's locus and nearly all of its final alignment:

- a matching read/template run emits `min(read_length, consensus_length)` exact matches plus the residual insertion or deletion;
- an HPC insertion/deletion emits the full raw run as `I`/`D`;
- blocks containing HPC mismatches or adjacent mixed operations are locally realigned in raw base space with the same pairwise backend;
- soft clipping represents read sequence outside the accepted layout interval.

This produces a provisional alignment to the expanded draft while avoiding a second seed-and-search mapping pass. A layout-interval-restricted pairwise rescue is allowed when the lifted trace fails CIGAR validation, but the tool never searches the rest of the contig, other contigs, or alternative paths.

While creating this provisional alignment, normalize every residual mismatch and indel into an allele observation. Classify and exclude changes that only alter the length of an existing homopolymer run: those have already been decided by the modal run-length caller. The remaining observations are substitutions and non-homopolymer indels that would alter the HPC base path.

If only CRAM is requested, the provisional alignment is already final and the tool skips directly to CRAM writing. Residual allele aggregation is enabled only by `--polish-non-hpc`.

### Phase F: optionally make only obvious non-HPC corrections

Use the same tier-priority rule as run-length calling. At each candidate locus, consider uniquely placed observations from the first tier with sufficient usable depth; lower tiers cannot outvote it.

A correction is accepted only when all of the following configurable checks pass:

- one normalized alternate allele has the required observation count and fraction;
- its support exceeds both reference support and the next alternate by the required margin;
- support is present on both strands when both strands cover the locus;
- the allele remains unchanged after local left/right realignment and normalization;
- the supporting columns are outside clipped, band-edge, or otherwise unstable alignment regions; and
- the candidate does not overlap another accepted or competing edit.

The default correction classes are substitutions and short non-homopolymer indels. Maximum indel length, minimum support, minimum allele fraction, and winning margin are user-visible/profile parameters rather than scattered constants. A homopolymer-only indel is never reconsidered here.

Apply all non-overlapping accepted edits simultaneously in one pass and write them to `assembly.non-hpc-edits.vcf.gz`, including winning tier, depth, allele counts, strand counts, confidence, and whether the edit changes the HPC string. There is no polishing iteration: conflicts and marginal sites remain unchanged and are reported as low confidence.

The result is:

```text
uncompressed_draft = expand(input_hpc_template, modal_run_calls)
final_sequence      = apply(accepted_non_hpc_edits, uncompressed_draft)
```

### Phase G: optionally re-render final alignments and write CRAM

When both optional stages are enabled, the provisional CIGAR cannot be written as the final CRAM record after the reference changes: CIGAR coordinates, `NM`, `MD`, `AS`, and the CRAM reference checksum must describe the corrected sequence.

Use the accepted-edit prefix map to lift unchanged parts of each cached trace to final coordinates. Re-align only small intervals overlapping accepted edits, then recompute the complete CIGAR and tags. This is a cheap second rendering of already placed reads, not a new mapping or consensus pass.

The output contract is stricter than a conventional mapper:

- at most one primary alignment record per accepted read;
- no secondary records (`0x100`) and no supplementary records (`0x800`) by default;
- reads that cannot produce one defensible primary alignment are omitted and reported;
- ambiguous/discarded reads are absent, not emitted as unmapped records;
- all emitted CIGARs consume the correct query/reference lengths and receive recomputed `NM`, `MD`, and `AS` tags;
- the original read name is preserved and every record has an `RG` tag identifying its tier.

The single-record default is intentional for downstream tools that interpret pileup depth. A future explicit `--allow-supplementary-within-placement` option may represent validated chimeric/split reads, but it is not part of the default contract.

## Standalone CLI

With the HPC template embedded in a newly built package, fast consensus is:

```text
hpc-consensus call \
  --import part001.cnspack \
  --tier 'hifi=hifi_*' \
  --tier 'ont=ont_*' \
  --output assembly \
  --threads 32
```

For an existing package without embedded tig bases, add the companion template:

```text
--template-hpc assembly.hpc.fasta
```

The small public configuration surface should be:

```text
--import FILE                 existing .cnspack input
--template-hpc FILE           required only if package tig bases are empty
--tier 'GLOB'                 repeat in descending quality order; LABEL=GLOB is optional
--unmatched error|drop        default error
--min-run-support N           default 1
--missing-run one|sidecar|error
--max-depth-per-tier N        optional deterministic cap
--alignment-band N           initial fallback band
--max-alignment-band N       reject beyond this point
--min-aligned-fraction F
--polish-non-hpc             off by default
--min-polish-support N
--min-polish-fraction F
--min-polish-margin F
--max-polish-indel N
--cram FILE                   off when omitted
--threads N
--tmp-dir DIR
--keep-traces                debugging only
```

Defaults live in one versioned configuration object and are written to the metrics file. Names affect priority only through the explicit ordered wildcard rules supplied on this invocation; there is no built-in prefix or technology behavior in the caller.

## Outputs

Required for the fast/default path:

```text
assembly.fasta
assembly.run-calls.tsv.gz
assembly.summary.json
```

With `--polish-non-hpc`:

```text
assembly.non-hpc-edits.vcf.gz
assembly.non-hpc-edits.vcf.gz.tbi
```

With `--cram assembly.cram`:

```text
assembly.cram
assembly.cram.crai
```

When requested, `assembly.cram` is coordinate sorted, has one `@RG` per tier, records the command/config digest in `@PG`, and embeds the final reference (`embed_ref=1`) so it remains portable. The tool also writes the final FASTA's `.fai` and sequence dictionaries/checksums needed to validate CRAM reference identity.

Recommended diagnostics:

```text
assembly.low-confidence.bed.gz
assembly.read-decisions.tsv.gz
assembly.hpc-disagreements.bed.gz
```

`read-decisions.tsv.gz` contains one row per input read with matched tier/rule, placement count/group, accepted/dropped state, alignment status, consensus observation counts, and a stable reason code. This provides read accounting without polluting CRAM with ambiguous or unmapped records.

## Implementation structure

C++17/20 is the lowest-friction fit with Verkko, LJA/MBG concepts, SIMD pairwise alignment, and htslib:

```text
hpc-consensus/
  cli/             call, inspect
  cnspack/         compatible tgTig/sqRead streaming reader
  template/        embedded or companion HPC template lookup
  tiers/           ordered original-name glob matching
  layout/          packaged placement and unique-only validation
  rle/             encoding, qualities, raw/run coordinate maps
  projection/      graph fast path and banded fallback alignment
  caller/          ordered-tier selection and sparse modal histograms
  polish/          normalized residual alleles and conservative edit calls
  expand/          FASTA emission and run/base prefix maps
  alignments/      trace lifting, local repair, CIGAR/tag validation
  cram/            htslib writing, external sort/merge, CRAI
  metrics/         run calls, low confidence, decisions, summary
  tests/           synthetic and replay corpus
```

Dependencies should remain small:

- the in-tree Verkko/Canu store and utility compatibility layer for `.cnspack` reading;
- htslib for companion FASTA and optional SAM/CRAM writing/indexing;
- ksw2 for pairwise alignment;
- zlib/libdeflate as inherited through htslib.

Do not copy Canu's consensus implementation into the new tool. The `.cnspack` compatibility layer should be limited to the minimum package/store records already used by Verkko and reviewed file-by-file for licensing; the consensus, caller, and alignment code remain independent. Canu is useful algorithmic context, but its repository includes GPL-covered components while Verkko's project code is public domain and LJA is BSD-3-Clause.

## How the existing implementations inform this design

| Existing code | Keep | Do not carry forward |
|---|---|---|
| Verkko layout/package | Native `.cnspack`, HPC coordinates, original names, uncompressed reads, partitioning | Hard-coded prefix meanings; random assignment of equal placements |
| MBG `HPCConsensus` | Exact HPC path projection; direct per-position expansion evidence | Graph-internal data structures as the only input format |
| LJA `homopolish` | Compress reads once; layout-restricted banded alignment; per-run evidence; local alignment only where needed | Median run caller; fixed vote-array limits; general local SPOA in v1 |
| Canu `utgcns` | Layout coordinates; robust alignment/CIGAR checks; final alignment metadata | Read-tiling template construction; retries and exception ladder; name-based quality; iterative DAG consensus |
| Shasta | Empirical modal repeat-count call | Whole assembler/marker-graph machinery |

## Parallelism and expected complexity

Let `R` be the total number of accepted read runs, `T` the number of template runs, and `b` the fallback alignment band.

- Directly projected reads cost `O(R)`.
- Fallback alignments cost approximately `O(R * b)` rather than pairwise full-matrix or POA cost.
- Modal calling and expansion cost `O(R + T)` with small sparse histograms.
- Residual evidence collection and trace lifting are `O(R)` plus small raw-space repair blocks.

Contigs and owned run chunks are independent. Threads should use work stealing over chunks. When CRAM is requested, each chunk writes deterministic target-ordered intermediate records and a final k-way merge produces coordinate-sorted CRAM without holding all alignments in memory.

The main speed win is structural, not a micro-optimization: identical packaged HPC slices and any packaged traces avoid alignment, all other reads are confined to their known layout interval, no POA graph is built, no consensus iteration is run, and no genome-wide final mapping is performed. The final post-edit work only re-renders cached traces and locally realigns reads touching accepted edits.

## Failure behavior

The tool should fail or flag, never silently improvise:

- non-HPC template: fail with the first offending coordinate;
- more than one accepted placement group for a read: fail because the layout contract was violated;
- reads marked ambiguous by the layout: ignore and record `AMBIGUOUS_PLACEMENT` when a rejection table is supplied;
- fallback alignment exceeds band/identity limits: drop and record `HPC_ALIGNMENT_FAILED`;
- no run evidence: use the configured `one`, `sidecar`, or `error` policy and mark the interval;
- empty or invalid winning-tier histogram: fail the chunk as an internal invariant violation;
- invalid lifted CIGAR when polishing or CRAM is requested: attempt one layout-interval-restricted rescue, then omit the read with `FINAL_ALIGNMENT_FAILED`;
- CRAM/reference checksum mismatch during optional CRAM validation: fail the run.

No random choices are permitted. Depth capping uses a stable hash of `(tier_id, read_id, read_name, tig_id, placement_start)` and spatial bins, so results are invariant to file and thread order.

## Verification plan

### Unit and property tests

- RLE encode/decode is lossless, including long runs greater than 255.
- Reverse-complement twice preserves bases, lengths, qualities, and offsets.
- `HPC(expand(template, calls)) == template` for all generated cases.
- When polishing is enabled, every difference between the uncompressed draft and final sequence has exactly one normalized VCF record with supporting evidence.
- `[5, 5, 5, 6, 6, 20]` calls 5.
- A higher-priority tier with any sufficient unique support cannot be outvoted by lower tiers.
- If tier 0 has no unique covering read, tier 1 is used; if both are absent, tier 2 is used.
- Ambiguous reads contribute zero observations and, when CRAM is requested, zero CRAM records.
- Input-order and thread-count permutations produce byte-identical FASTA/run calls and, when requested, semantically identical CRAM.
- When CRAM is requested, every emitted read has exactly one primary and no secondary/supplementary alignment.
- When CRAM is requested, lifted CIGAR query/reference consumption, `NM`, and `MD` agree with an independent validator.
- Disabling non-HPC polishing makes the final sequence identical to the invariant-preserving uncompressed draft.
- With neither opt-in, no base-space aligner or CRAM writer is invoked.

### Synthetic biological cases

- unimodal, tied, skewed, and bimodal homopolymer distributions;
- name-selected tiers with deliberately contradictory run lengths;
- long homopolymers, contig ends, reverse-strand reads, clipping, and `N` runs;
- HPC substitutions/indels adjacent to long runs;
- duplicated loci with one, two, or no unique placements;
- reads crossing chunk boundaries;
- low/no-coverage runs and sidecar fallback.

### Assembly benchmarks

Compare the new tool with current Verkko/Canu consensus, and include LJA/MBG behavior where formats permit:

- Merqury/yak/meryl QV and k-mer completeness;
- substitutions and indels against trusted haplotype references;
- homopolymer length error stratified by base, true length, coverage, and winning tier;
- uncompressed-draft invariant violations (must be zero) and audited final HPC changes;
- consensus changes caused by input permutation (must be zero);
- discarded ambiguous-read counts and repeat-region pileup depth;
- CRAM validity, one-placement compliance, mapped accepted-read fraction, and CIGAR edit distance;
- wall time, CPU time, peak RSS, temporary disk, and direct-projection fraction.

The key ablations are mode versus median/mean, direct projection versus fallback alignment, strict ordered tiers versus pooled evidence, and trace-lifted CRAM versus an independent contig-constrained remap.

## Migration into Verkko

### Milestone 1: standalone replay prototype

Implement:

1. streaming import of existing `.cnspack` records and original read names;
2. embedded/companion HPC template validation;
3. ordered tier wildcards and one-placement validation;
4. RLE encoding and banded HPC alignment;
5. strict ordered-tier modal calling;
6. FASTA, run calls, decisions, and summary output.

Use existing `.cnspack` partitions plus a companion graph-derived HPC FASTA. Do not alter the Verkko pipeline during this milestone.

This milestone is the complete fast consensus mode; it has no correction or CRAM requirement.

### Milestone 2: single-file package path

Populate the already-serialized `tgTig::_bases` field with each graph-derived HPC template during package construction, so `--template-hpc` is no longer needed for new runs. Preserve backward compatibility with old packages. Use any packaged child CIGAR/delta when available and benchmark the fraction of reads requiring dynamic programming.

### Milestone 3: optional correction and CRAM

Add base-space residual evidence, the conservative edit caller, trace-derived final alignments, and CRAM/CRAI behind the independent `--polish-non-hpc` and `--cram` switches.

### Milestone 4: opt-in Verkko integration

Add a pipeline switch such as:

```text
--consensus-engine hpc-consensus
```

The integration builds the HPC templates from final graph paths, stores them in the existing tig-bases slot while creating each `.cnspack`, invokes the standalone executable with ordered tier wildcards, and always consumes its FASTA; it consumes CRAM only when that output was requested. Keep the current Canu path as a comparison fallback during validation.

### Milestone 5: default after evidence

Make the new path default only after it meets these release gates:

- no regression in global QV, completeness, or phasing on the replay panel;
- statistically supported improvement or non-regression in homopolymer length accuracy;
- zero HPC compression-invariant failures;
- deterministic output;
- one-placement CRAM contract verified whenever CRAM output is enabled; and
- a material runtime/RSS reduction relative to `utgcns`.

## Recommended v1

The smallest useful implementation is not an abPOA wrapper. It is:

```text
embedded or companion fixed HPC template
  + existing .cnspack layout and raw reads
  + explicit ordered name-wildcard tiers
  + direct projection / narrow banded pairwise fallback
  + per-tier run-length histograms
  + first-available-tier mode calls
  + optional audited non-HPC correction
  + optional trace-derived single-placement CRAM
```

That design directly matches the data Verkko already has, implements the requested homopolymer and read-priority semantics without built-in name conventions, and removes the expensive template-building problem instead of optimizing it.
