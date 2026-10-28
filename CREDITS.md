# Credits and implementation provenance

This document separates scientific influences, code actually used, and tools
used to evaluate the caller. It describes the implementation audited on
2026-10-27, not every proposed backend in the original design document.
Attributions below identify specific dependencies, recorded design references,
and evaluation resources. Linked publications and repositories identify their
authors and contributors.

## Scientific and design influences

| Project | Contribution to this design | Relationship to the implementation |
| --- | --- | --- |
| **Verkko** | Graph-derived HPC templates, read layouts and partitioned packages; mature consensus supplied the accuracy baseline and motivated unique-read preference and coverage fallback. | Direct integration with native package/store code; independent consensus policy. [Paper](https://doi.org/10.1038/s41587-023-01662-6), [source](https://github.com/marbl/verkko). |
| **Canu / Celera Assembler** | Layout-oriented consensus, read clipping, placement multiplicity, and alignment bookkeeping. Canu's `utgcns` was the main existing caller examined during accuracy work. | Uses Canu's native `tgTig`/`sqStore` readers. Does not call or reproduce its read-tiling template builder or DAG consensus engine. [Canu paper](https://doi.org/10.1101/gr.215087.116), [source](https://github.com/marbl/canu). |
| **MBG** | The recorded design specifically cites `HPCConsensus` and `ConsensusMaker`: project reads onto an HPC graph/path and accumulate expansion evidence. | Architectural influence and upstream source of graph templates in Verkko, not a library called by this module. [Paper](https://doi.org/10.1093/bioinformatics/btab004), [source](https://github.com/maickrau/MBG). |
| **LJA / LJApolish** | The design cites homopolymer decompression using layout-constrained alignment and per-position run-length observations. | Conceptual influence; this module does not use LJA's median caller, ksw2 backend, or local SPOA machinery. [Paper](https://doi.org/10.1038/s41587-022-01220-6), [source](https://github.com/AntonBankevich/LJA). |
| **Shasta** | The original design explicitly cites empirical modal repeat-count calling. Shasta's modal and probabilistic run-length callers are also relevant prior art for separating base order from repeat counts. | Algorithmic precedent, not reused code or model files. Our optional estimator is a per-tier conditional-probability maximum-likelihood caller without a length prior; it is not a port of Shasta's Bayesian caller. [Paper](https://doi.org/10.1038/s41587-020-0503-6), [repeat-count documentation](https://paoloshasta.github.io/shasta/ComputationalMethods.html). |

## Code and libraries actually used

### Verkko, Canu and meryl-utility

The caller links the externally built `libverkko.a`. Its package adapter in
`src/main.cpp` calls `tgTig::loadFromBuffer` and
`sqStore::sqStore_loadReadFromBuffer`; native read decoding, serialization and
buffering are provided upstream. They have **not** been reimplemented here.
The adapter controls which records are retained and how their ownership is
shared, but the package ABI and encoding remain Canu's.

[meryl-utility](https://github.com/marbl/meryl-utility) supplies buffering,
sequence/bit utilities and the bundled alignment code. Its source records the
Canu, Celera Assembler and kmer-package lineage.

### Edlib

`src/core.cpp` and `src/offsets.cpp` directly call the **meryl-utility version of
Edlib**, through `merylutil::edlibAlign`. Edlib performs the unit-edit-distance
global, infix and prefix alignments, distance calculations and traceback.
The wrapper's interval selection, error limits, exact-match shortcut and
coordinate translation are local code; the alignment algorithm is Edlib's,
not a local reimplementation of Myers's bit-vector algorithm.

Cite Šošić M, Šikić M. *Edlib: a C/C++ library for fast, exact sequence alignment
using edit distance*. Bioinformatics 33:1394–1395 (2017).
[Paper](https://doi.org/10.1093/bioinformatics/btw753),
[upstream source](https://github.com/Martinsos/edlib).
The inspected bundled source carries the MIT notice, copyright 2014
Martin Šošić.

### SAMtools, HTSlib and HTSCodecs

The caller writes SAM text and computes its CIGAR/NM/MD fields locally. For
`--cram`, it invokes **SAMtools** for `faidx`, CRAM sorting/encoding, indexing
and `quickcheck`. CRAM compression and indexing are not implemented here.
There is a second dependency path: the inspected executable also contains
HTSlib and HTSCodecs symbols pulled in transitively from `libverkko.a`, even
though this project's export path uses the SAMtools command line rather than
direct HTSlib API calls. An archive's contents and what the linker includes can
vary with the Verkko build, so inspect the binary being distributed.

- Danecek et al. *Twelve years of SAMtools and BCFtools* (2021),
  [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008).
- Bonfield et al. *HTSlib: C library for reading/writing high-throughput
  sequencing data* (2021),
  [doi:10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007).
- [SAMtools publications and format credits](https://www.htslib.org/doc/#publications),
  [HTSlib source](https://github.com/samtools/htslib),
  [HTSCodecs source](https://github.com/samtools/htscodecs).

### Supporting libraries and build tools

The Makefile links the following in addition to `libverkko.a`. They supply
infrastructure rather than the consensus estimator:

| Dependency | Role |
| --- | --- |
| OpenSSL (`libssl`, `libcrypto`) | Inherited archive dependency. [Project](https://openssl-library.org/). |
| libcurl | Inherited I/O dependency. [Project and notice](https://curl.se/docs/copyright.html). |
| liblzma | XZ Utils compression support. [Project](https://tukaani.org/xz/). |
| zlib | Compression support. [Project](https://zlib.net/). |
| bzip2 | Compression support. [Project](https://sourceware.org/bzip2/). |
| OpenMP runtime | The runtime selected by the compiler, for example LLVM libomp or GCC libgomp. Required by the linked Verkko build. The caller's own worker scheduler uses `std::thread`. |
| C/C++ and operating-system libraries | Standard containers, threads, allocation, POSIX processes and `fnmatch`. The platform libraries vary by build. |
| GNU Make, pkg-config, compiler toolchain, Python | Build/test infrastructure. Python scripts use the standard library; Python is not required for a consensus invocation. |

The standalone Makefile follows Verkko's build-directory and compiler-flag
conventions. It does not include or vendor Verkko's boilermake implementation.

## Functionality implemented in this repository

The following are implemented in this repository: RLE and run/base coordinate
maps; layout-offset inference and rescue; ordered tier selection; sparse
run-length histograms and tie-breaking; conditional-probability table loading
and likelihood scoring; targeted HP re-estimation; linked read-local event
normalization; conservative allele/strand voting; end recovery; CIGAR rendering,
coordinate lifting and SAM tag generation; and bounded worker scheduling.

These implement common operations or project-specific policies; their presence
in this repository is not a claim of algorithmic originality. For example,
CIGAR/NM/MD generation follows SAM conventions and does not replace HTSlib's
file codec.

The run-length probability tables were calibrated in the project's experiments
using supplied reads and Q100 truth, rather than imported from Shasta or LJA.
Their training/validation manifests and data versions belong with the retained
experiment provenance. Do not infer a specific Q100 release from the filename
of a truth slice alone.

## Research context and evaluation tools

POA and affine-gap work discussed in `CONSENSUS_V2_DESIGN.md` supplied research
context. **POASTA, minipoa, Theseus, abPOA, SPOA, ksw2 and WFA/WFA2 are not
implemented alignment/consensus backends in this caller.** A proposed backend
is not a dependency and should not be credited as code actually used.

Development also benefited from tools and data that are not required for an
ordinary consensus run:

- **Merqury and Meryl** — k-mer-based accuracy and completeness evaluation.
  [Merqury paper](https://doi.org/10.1186/s13059-020-02134-9),
  [Meryl source](https://github.com/marbl/meryl).
- **MUMmer / nucmer / dnadiff** — assembly-to-truth
  alignments and difference analysis. [Project](https://github.com/mummer4/mummer).
- **FastANI** — early E. coli ANI comparisons.
  [Project](https://github.com/ParBLiSS/FastANI).
- **seqtk** — telomere checks during extension work.
  [Project](https://github.com/lh3/seqtk).
- **T2T-HG002 / Q100** — truth-guided calibration
  and validation; supplied mature Verkko assemblies were the comparison caller.
  [HG002 project and release information](https://github.com/marbl/HG002).

## Licensing and audit scope

The project's original code remains MIT-0. That does not relicense linked or
external components. Scientific citations and copyright/license notices serve
different purposes; preserve the notices for the exact dependency versions
included with a release.

Relevant upstream records include Verkko's `README.licenses`, Canu's
`README.licenses`, meryl-utility's `README.licenses`, Edlib's source-header MIT
notice, and the file-level notices in bundled HTSlib/HTSCodecs. Canu's repository
contains multiple licensing categories; a single blanket label is insufficient
for all files. The inspected `tgTig`/store headers identify their public-domain
lineage, while bundled Edlib, HTSlib and codecs retain their separate notices.
This credit inventory is not an exhaustive binary-distribution license manifest.

Audit baseline:

- Caller: `d6328f739948c7bdfbc9f5336c3286ba352f2e7d`.
- Verkko: `cb69f969eebfa2b2a74003a120edabf2b084248f`.
- Canu: `dacd2f07c40dfc408ac2fa6767dba7f640f30064`.
- meryl-utility: `f800fc4dada365a701122b5ac1ecb946eb402cb4`.
- MBG examined for design context: `2e1d44b7a001e35d2c40c31520d11f7c621a7506`.

Evidence comprised the current source/API calls, recorded design and project
history, upstream file notices and publications, and `nm -C`/`otool -L` on the
local macOS executable. Design references establish documented precedents,
not that their implementations were copied. This inventory is not an exhaustive
source-provenance or license audit, and does not establish originality.
