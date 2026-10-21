# Testing hpc-consensus

## Automated tests

After building, run from the module directory:

```sh
ctest --test-dir build --output-on-failure
```

`core` tests run encoding (including long homopolymers), reverse complement,
name glob ordering, run-mode/tie selection, tier fallback/support, alignment
operation conventions, NM/MD generation, substitution and insertion correction,
homopolymer exclusion, and coordinate lifting.

Boundary regressions cover left/right/two-sided HPC overhang clipping, retained
boundary reads in native packages, soft-clipped CRAM records, and final-reference
NM/MD validation with `samtools calmd`.

Memory-refactor regressions compare one-read-sized and default batches, including
multi-tier fallback with minimum run support greater than one. Core tests cover
packed traces and compare indexed polishing depths with a brute-force scan over
100 deterministic randomized read sets. Biological equivalence is checked against
a preserved pre-refactor executable on E. coli with polishing and CRAM enabled.

`native-package` uses the selected Verkko checkout's `sqStoreCreate` and
`layoutToPackage` executables to generate real packages in a temporary directory.
It checks tier priority, ambiguous-read removal, reverse orientation, deterministic
thread counts, unmatched-read and existing-output errors, all four combinations
of correction/CRAM, and an actual incorrect-template-base correction. With
samtools available, it decodes CRAM, checks single primary placements and indexes,
and checks corrected-reference tags with `samtools calmd`.

CRAM checks are skipped if samtools is absent; a passing test in that environment
does not validate export. These synthetic fixtures are regression checks, not
biological accuracy or performance benchmarks.

## E. coli test bundle

Please provide a directory (or paths) containing:

1. The `.cnspack` file(s) to process, from the compatible Verkko build.
2. The corresponding HPC path/template FASTA in the exact packaged layout
   coordinate system. Use numeric tig IDs, or include the tig-name-to-ID map.
3. The ordered read-name wildcard tiers, including how unmatched reads should
   be handled. Example names from each technology help verify the patterns.
4. Optionally, the existing Verkko consensus for a direct baseline comparison.
5. Optionally, a trusted same-strain reference for accuracy assessment, including
   plasmids if relevant. A different-strain reference cannot distinguish all
   real variation from consensus errors.

Also identify the Verkko revision/build that made the packages and any known
layout trimming or coordinate shifts. A raw read FASTA is not needed because
the reads are in the package. Keep input files unchanged; runs use fresh prefixes.

Once supplied, the planned checks are:

- Confirm template IDs, lengths, HPC form, and placement compatibility.
- Run fast consensus; inspect unsupported runs, read rejection reasons, and
  tier usage, and verify that recompressing the output recovers the template.
- Run correction only, CRAM only, and both, recording elapsed time and peak
  memory with a platform-appropriate external timing tool.
- Compare sequences across modes and audit every reported correction against
  the saved draft. CRAM-only must not change the consensus.
- Decode/index-check CRAM, verify unique primary read names, CIGAR consumption,
  and NM/MD against the final FASTA.
- Compare with the baseline/reference when available, reporting substitutions,
  indels, and homopolymer errors separately and accounting for circular origin,
  orientation, and genuine strain differences.

This checklist describes how to repeat a biological evaluation; synthetic test
success alone is not a biological quality claim. Keep benchmark inputs and
results outside the tracked test suite.

## Maintenance and equivalence tests

`maintenance` runs 1,000 deterministic differential cases for insertion
normalization and edit application, covering packed/expanded traces, long repeat
shifts, substitutions, deletions and insertions. It also checks malformed edit
intervals, numeric parsing and FASTA wrapping. Checks stay enabled in Release.
Worker tests verify exactly-once execution and joined propagation of allocation
failures. Sanitizer builds are described in `../MAINTAINING.md`.

`linked-events` checks equivalent split/contiguous CIGAR representations,
combined allele size limits, rejection of synthetic haplotypes assembled from
different reads, and 2,000 deterministic randomized trace-to-sequence roundtrips.
It verifies that jointly normalizing linked events preserves the actual read
sequence, including packed traces, substitutions, insertions and deletions.

`compare_versions.py` accepts two executables and common native input arguments.
It compares all generated text artifacts and decoded CRAM records, retaining
timing logs and a JSON report in a new directory. See `../MAINTAINING.md` for
usage, module boundaries, coordinate contracts and benchmarking cautions.

Both biological end directions have synthetic regressions: molecule-count
support, reverse orientation, automatic/explicit/disabled extension, side flags,
internal regional boundaries and CRAM NM/MD validation.

## Optional matched-read graph diagnostic

Configure with `-DHPC_BUILD_GRAPH_PROBE=ON` to build `hp-graph-probe` and enable
its synthetic `graph-probe` CTest. This requires installed Boost headers; pass
`-DHPC_BOOST_INCLUDE_DIR=/path/to/include` if needed. The target compiles the
actual external Verkko `AlnGraphBoost.C`, without modifying or relicensing it.
The production executable does not link this graph implementation.

Usage: `hp-graph-probe reference.txt reads.tsv run_begin run_end`. The reference
is a single uppercase sequence; each TSV row is `NAME TIER ORIENTED_SEQUENCE`
with diagnostic tiers LA/UL/LAO. Inputs must already be unique, fully spanning
read chunks with distinct names. Run coordinates are zero-based, half-open.
Use small local windows, not chromosome-length sequences. The output records
per-read HPC/raw evidence, selected read IDs, and local graph sequences.

The probe compares exclusive fallback with additive fallback that keeps every
higher-priority unique read. Qualification remains five reads in the fallback
tier, not five cumulatively. It also emits all-unique and LA+LAO mixture controls.
Every estimator within a policy uses the same HPC-eligible read IDs. Raw alleles
with interrupted or mismatching immediate flanks are reported, not silently
classified as clean runs. Graph inputs retain those reads.

This is not a complete utgcns replay: it uses fixed local read chunks and a fixed
raw template, one graph round, no whole-read promotion, and no ambiguous reads.
Its results isolate local evidence/algorithm differences but cannot establish
which reads contributed to a separately generated mature consensus. Never
interpret truth-selected diagnostic loci as genome-wide QV validation.
# Bounded borrowing regression coverage

The `core` tests also cover internal offset anchors, repeated-template rejection,
one-read-one-vote clustering, competing clusters, and signed sparse offsets.
`native-offset` generates one native package with two tigs whose reads overhang
both template ends. The second has five agreeing reads, one offset outlier,
and one uninformative read. It checks whole-package and region recovery with
`--no-extend-ends`, so success on just the first tig is insufficient.
Neither test needs external genome data.

`borrow` checks unmodeled confidence protection and bounded ambiguous support.
`native-borrow` builds a synthetic native cnspack and checks that the opt-in
policy changes an uncertain run while the default retains its original call.
`apply-borrow` is a histogram driver for independent Python/C++ parity checks;
it is not a biological quality benchmark. Biological paired runs must compare
polished outputs on common aligned reference spans, with mixed-evidence masks
selected without truth alleles. A single-haplotype truth does not measure phase
switches.
