# Maintaining hpc-consensus

Start with `CONSENSUS_V2_DESIGN.md` for biological policy and `README.md` for
the CLI contract. This guide explains where to make changes without accidentally
changing that policy.

## Module map

| File | Responsibility |
| --- | --- |
| `src/options.hpp`, `src/options.cpp` | Configuration defaults, CLI parsing and validation |
| `src/io.hpp`, `src/io.cpp` | FASTA, template-name map, literal-argument subprocess execution |
| `src/core.hpp` | Algorithm interfaces, coordinate conventions and ownership contracts |
| `src/core.cpp` | RLE, alignment adapters, local HPC filtering and empirical run estimator |
| `src/models.cpp` | Calibration file validation and likelihood estimator |
| `src/offsets.hpp`, `src/offsets.cpp` | HPC origin anchors, repeat-aware independent-read clustering and overlap validation |
| `src/polish.cpp` | Raw-base candidate normalization, voting, edits and alignment tags |
| `src/ends.cpp` | Coverage-supported extension at biological template ends |
| `src/parallel.hpp` | Bounded worker scheduling and joined exception propagation |
| `src/main.cpp` | Native package ownership and pipeline stage orchestration |

The pipeline still owns the tier-by-tier dense/sparse evidence accumulator in
`main.cpp`. It is intentionally not replaced with a generic graph framework or
one pileup object per chromosome base. Further extraction should preserve its
storage lifetimes and original placement order, not just its final sequence.

## Invariants to preserve

- HPC run indexes, raw draft coordinates, and final sequence coordinates are
  different spaces. Internal ranges are zero-based, half-open; SAM/VCF positions
  are one-based. Run diagnostics describe the original unextended draft.
- Run calling tries qualifying unique tiers in user order, then qualifying
  pooled evidence within each tier, then sparse fallback. No cross-tier pooling.
  End extension deliberately has no sparse fallback. Polishing has its own
  documented allele-fraction/strand gates; do not silently unify the policies.
- Read names determine only the configured tier. Placement counts determine
  uniqueness. Repeated placements share immutable bases; the first successfully
  exported occurrence is primary. Sorting the processing index is safe; sorting
  the original placement vector changes this output contract.
- `Alignment::bases()` includes soft clips and is oriented to the reference.
  `ops` may be packed or expanded according to `packed_ops`. Read RLE and expanded
  traces are temporary, not chromosome-lifetime objects.
- Models, alignment scoring, tie-breaking, filtering, and voting thresholds
  affect QV. Treat changes there as algorithm changes, not cleanup.
- cnspack is native ABI-dependent. Keep using the matching Verkko reader and
  build, and do not change third-party source as part of local refactoring.

## Development loop

1. Preserve a known-good executable before a behavior-preserving change.
2. Add a small deterministic regression for the affected contract. Prefer
   comparing optimized code with a simple reference algorithm over timing tests.
3. Format touched C++ with `.clang-format`, build, and run
   `ctest --test-dir build --output-on-failure`. C++ checks throw even in Release.
4. Use `tests/compare_versions.py` on a representative native fixture, with
   polishing and CRAM when relevant. The script requires a fresh work directory,
   preserves logs, compares text artifacts and decoded SAM records byte-for-byte,
   and checks CRAM integrity. Do not compare compressed CRAM bytes or path-bearing
   headers as a substitute for alignment equivalence.
5. For a biological algorithm change, additionally measure QV and compare to
   independent truth. Matching a baseline executable only proves equivalence,
   not that either sequence is biologically correct.

For memory/undefined-behavior checks, configure a separate build with
`-DHPC_ENABLE_SANITIZERS=ON -DCMAKE_BUILD_TYPE=Debug`, retaining the same
`VERKKO_ROOT` and compiler configuration. Run CTest there as well. This
instruments this project's code, not the already-built Verkko archive. A
production-only build can set `-DBUILD_TESTING=OFF` to omit test binaries and the
Python requirement.

Example (common input/tier flags follow the final `--`):

```sh
python3 tests/compare_versions.py \
  --before build/hpc-consensus-before-maintenance \
  --after build/hpc-consensus --work /path/to/new-comparison --cram -- \
  --import fixture.cnspack --template-hpc fixture.fasta \
  --tier 'LA=hifi_m*' --tier-model LA=models/revio.run-model.tsv \
  --tier 'LAO=hifi_[0-9a-ln-z]*' --tier 'UL=ont_*' \
  --threads 4 --polish-non-hpc
```

## Performance and memory

Measure sequential runs with identical Release compiler settings, input, model,
options and thread count. Repeat and alternate order for small differences;
filesystem cache and concurrent work can dominate short tests. Keep wall time,
CPU time and peak RSS distinct. macOS `time -l` reports RSS in bytes; Linux
`time -v` reports KiB. Caller `[memory]` lines are high-water marks, not live RSS,
and exclude samtools' concurrent allocation. Batch size is not an RSS limit.

The maintenance pass avoids reference-sized temporary substrings while applying
edits or writing draft FASTA, and transfers the draft buffer when no edits exist.
Insertion normalization tracks rotations and materializes the rotated allele
once. These changes preserve votes and sequence. Their effect on total RSS
depends on which earlier stage sets the peak; local savings need not lower the
process high-water mark.

Unexpected worker allocation failures now propagate after all workers join,
instead of being recorded as ordinary read-alignment failures. This is a
failure-handling change: resource exhaustion must fail the run, not lower QV by
silently removing evidence. Incomplete staging output remains available for
diagnosis and is not published as completed output.

Remaining major costs are native read storage, dense per-run arrays, sparse
multi-allelic/low-depth evidence, diagnostic output, and raw edlib realignment.
Replacing whole-read realignment with local windows or changing alignment ties
could affect QV and is not part of this behavior-preserving pass. Full-chromosome
performance and quality must still be verified on the target machine.
