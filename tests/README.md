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

No E. coli test has been run yet. This checklist does not imply a benchmark
result or authorize changes to the supplied source data.
