// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstddef>
#include <cstdint>
#include <map>
#include <istream>
#include <memory>
#include <string>
#include <vector>

namespace hpc {
// All intervals in this API are zero-based and half-open. RLE indexes count HPC
// runs; Alignment/Edit coordinates count uncompressed reference bases.
// offsets has one sentinel: offsets.back() equals the original raw length.
struct Rle {
  std::string bases;
  std::vector<uint32_t> lengths;
  std::vector<size_t> offsets;
};
// Accept uppercase ACGTN only; throw for invalid bases or uint32 run overflow.
Rle encode(const std::string &s);
std::string reverse_complement(std::string s);
// Each operation consumes one query and/or reference base: = X I D.
// Global unit-cost alignment, capped at ceil(max_error * max(sequence lengths)).
// Throws on failed alignment. Returned traces are expanded, not CIGAR strings.
std::string align(const std::string &query, const std::string &reference, double max_error);
// Centered local-quality gate in expanded HPC alignment-operation space.
// X/I/D all count as errors. Terminal windows are truncated to the trace.
// Empty means every matching operation is eligible (including exact traces).
// Otherwise one byte per operation: only locally acceptable '=' positions vote.
std::vector<uint8_t> local_hpc_votes(const std::string &ops, size_t radius, double max_error);
// Persist compact traces; expand only for the current read. Internal packed
// traces contain =/X/I/D, not soft clips or arbitrary externally supplied SAM.
std::string cigar(const std::string &ops);
std::string expand_cigar(const std::string &packed);
// Query bounds are HPC run indexes, translated to raw clips via Rle::offsets.
struct BoundaryAlignment {
  std::string ops;
  size_t query_begin = 0, query_end = 0;
};
// Consume the whole template overlap, with free read overhangs on flagged sides.
// Tied best endpoints are unresolved and throw, rather than choosing a repeat.
BoundaryAlignment align_boundary(const std::string &query, const std::string &reference, bool left,
                                 bool right, double max_error);
// Return the unique best full-query infix placement, or -1 if unresolved.
int64_t locate(const std::string &query, const std::string &reference, double max_error);
// Require a dominant (>=75%) tolerance-width cluster and return its lower median.
// Estimates are signed layout-start minus template-aligned-start coordinates.
int64_t offset_cluster(std::vector<int64_t> estimates, size_t required = 3, int64_t tolerance = 20);
struct Tier {
  std::string label, pattern;
};
int tier_of(const std::string &name, const std::vector<Tier> &tiers);
// tier == -1 denotes no evidence; length one is the explicit fallback, not a vote.
struct Call {
  uint32_t length = 1, depth = 0, peak = 0;
  int tier = -1;
};
using Histogram = std::map<uint32_t, uint32_t>;
// Rounded mean is an opt-in diagnostic estimator; mode remains the default.
Call call(const std::vector<Histogram> &tiers, uint32_t min_support, bool rounded_mean = false);
// Conditional run-length likelihoods. Each base has contiguous true-length
// rows beginning at one; the last observed bin includes all longer observations.
// Models are optional and never change tier eligibility or placement preference.
class RunModel {
public:
  static RunModel read(std::istream &input);
  uint32_t estimate(char base, const Histogram &observations, uint32_t fallback) const;
  uint32_t estimate_single(char base, uint32_t observed, uint32_t depth) const;
  // Experimental same-technology borrowing. An empty model uses modal support.
  // Unique calls with likelihood ratio >=100 (or modal support >=90%) are
  // protected. Borrow <=min(5, unique_depth/2) equivalent ambiguous votes only
  // if their candidate is compatible with unique evidence (LR<=20 / >=20%).
  bool uncertain(char base, const Histogram &unique) const;
  uint32_t borrow(char base, const Histogram &unique, const Histogram &ambiguous,
                  uint32_t current) const;

private:
  std::vector<double> likelihood_scores(char base, const Histogram &observations) const;
  size_t observed_max = 0;
  uint32_t minimum_depth = 5;
  std::map<char, std::vector<std::vector<double>>> log_probabilities;
  // Single-allele predictions are independent of depth once the minimum is met.
  std::map<char, std::vector<uint32_t>> single_allele_calls;
};
// Replace reference[begin:end] with alt. Empty intervals represent insertions.
struct Edit {
  size_t begin = 0, end = 0;
  std::string alt;
  bool operator<(const Edit &o) const;
};
// bases() returns the entire read oriented to the reference, including soft clips.
// begin/end exclude those clips. ops consumes only the aligned raw-base interval.
// shared_sequence owns immutable storage in production; sequence supports small
// independent alignments/tests. The shared representation takes precedence.
struct Alignment {
  uint32_t id = 0;
  std::string name, sequence, ops;
  size_t begin = 0, end = 0;
  size_t clip_left = 0, clip_right = 0;
  int tier = 0;
  bool reverse = false;
  bool packed_ops = false;
  std::shared_ptr<const std::string> shared_sequence;
  const std::string &bases() const {
    return shared_sequence ? *shared_sequence : sequence;
  }
};
// Only template-induced soft clips may be recovered; packaged raw clipping is
// excluded by the caller. Pointers borrow existing oriented read storage.
struct EndRead {
  Alignment *alignment;
  bool unique;
  size_t left, right;
  bool eligible = true;
};
struct EndExtension {
  size_t left = 0, right = 0;
  size_t left_candidates = 0, right_candidates = 0;
};
// Extend real template ends with anchored overhangs, requiring support at EVERY
// added run (no sparse fallback). Update coordinates, clips and compact traces.
// Internal reference sequence is unchanged. Left tails are reversed, not
// complemented. Memory depends on bounded terminal tails, not chromosome size.
EndExtension extend_ends(std::string &reference, const std::vector<EndRead> &reads, size_t tiers,
                         const std::map<int, RunModel> &models, unsigned support, size_t max_bases,
                         double max_error, bool ambiguous, bool left, bool right);
// Intersect an alignment with reference[begin:end], rebase reference coordinates,
// and soft-clip excluded query bases. Insertion at begin belongs to the crop;
// insertion at end does not. Returns false if no query/reference overlap remains.
bool crop_alignment(Alignment &a, size_t begin, size_t end);
// Extract read-supported replacements. When left-normalization would cross an
// earlier raw event, combine their original haplotype before normalizing. Never
// join pieces from different reads. Unlinked neighboring changes stay separate.
std::vector<Edit> events(const Alignment &a, const std::string &reference);
struct Accepted {
  Edit edit;
  int tier;
  uint32_t depth, support;
};
std::vector<Accepted> polish(const std::string &ref, const std::vector<Alignment> &alns,
                             size_t tiers, uint32_t support, double fraction, size_t max_indel,
                             size_t retain_unique_tiers = 0);
// Borrowed alignments must remain alive for this call. Callers prefilter consensus
// eligibility; returned edits are ordered, non-overlapping, and in draft space.
// Optional experimental retention: ranks [0, retain_unique_tiers) are unique.
// A qualified unique tier includes higher ranks' evidence in fraction/strand
// checks. Its own depth must still qualify. Ambiguous fallback is unchanged.
std::vector<Accepted> polish(const std::string &ref, const std::vector<const Alignment *> &alns,
                             size_t tiers, uint32_t support, double fraction, size_t max_indel,
                             size_t retain_unique_tiers = 0);
std::string apply(const std::string &ref, const std::vector<Accepted> &edits);
// Map a draft boundary through sorted edits. Points inside a replacement map to
// its start; a point at an insertion maps after inserted bases. Touching reads
// must subsequently be realigned, not merely assigned a shifted CIGAR.
size_t lift(size_t p, const std::vector<Accepted> &edits);
// Validate trace consumption/exact matches and calculate final-reference MD/NM.
std::string md(const Alignment &a, const std::string &ref, size_t &nm);
} // namespace hpc
