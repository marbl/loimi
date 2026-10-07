# Project engineering expectations

These instructions apply throughout this repository.

## Maintainability and documentation

- Design for maintainability from the start: prefer focused functions, explicit
  ownership, clear interfaces, and descriptive names over dense or clever code.
- Document public APIs and non-obvious algorithms. Explain purpose, coordinate
  systems, units, preconditions, failure behavior, and important invariants.
- Use implementation comments to explain why an approach is needed, especially
  edge cases and tradeoffs. Avoid comments that merely restate the next line.
- Keep comments, CLI help, README, and tests synchronized with behavior changes.
  Document limitations honestly; distinguish demonstrated behavior from plans.
- Follow `.clang-format`; format C++ sources and tests after editing. Avoid
  unrelated behavioral refactors in documentation-only changes.

## Consensus-specific invariants

- Distinguish HPC/layout coordinates from raw/draft/final base coordinates.
  Internal intervals are zero-based, half-open; SAM/VCF output positions are
  one-based. Preserve orientation and clipping semantics explicitly.
- Separate consensus eligibility from alignment export. Preserve user-ordered
  tiers, placement-count uniqueness, primary/secondary rules, and deterministic
  mode/tie selection unless the user requests a change.
- Protect chromosome-scale memory behavior: share immutable read storage, keep
  per-read RLE/expanded traces temporary, preserve bounded alignment batches,
  and avoid silent deep copies or per-read-per-reference quadratic work.
- State remaining scaling costs. A batch-size limit is not a process RSS cap.
- Treat cnspack as a trusted, ABI-dependent format. Use native Verkko readers;
  do not invent an undocumented binary format or silently guess associations.

## Verification and change boundaries

- Run the CMake build and `ctest --test-dir build --output-on-failure` after
  code changes. Add regression coverage for new behavior and fixed edge cases.
- For behavior-preserving performance changes, compare consensus, diagnostics,
  and decoded alignments against the previous implementation where practical.
  Report measured memory/runtime separately from estimates and untested scales.
- Preserve unrelated user changes, source assemblies, and external Verkko
  checkouts. Keep test outputs separate; do not commit datasets or build products.
- This project's original code is MIT-0. Preserve third-party licenses and
  notices; the project license does not relicense linked dependencies.
