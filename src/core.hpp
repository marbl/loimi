// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstdint>
#include <map>
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
Call call(const std::vector<Histogram> &tiers, uint32_t min_support);
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
std::vector<Edit> events(const Alignment &a, const std::string &reference);
struct Accepted {
  Edit edit;
  int tier;
  uint32_t depth, support;
};
std::vector<Accepted> polish(const std::string &ref, const std::vector<Alignment> &alns,
                             size_t tiers, uint32_t support, double fraction, size_t max_indel);
// Borrowed alignments must remain alive for this call. Callers prefilter consensus
// eligibility; returned edits are ordered, non-overlapping, and in draft space.
std::vector<Accepted> polish(const std::string &ref, const std::vector<const Alignment *> &alns,
                             size_t tiers, uint32_t support, double fraction, size_t max_indel);
std::string apply(const std::string &ref, const std::vector<Accepted> &edits);
// Map a draft boundary through sorted edits. Points inside a replacement map to
// its start; a point at an insertion maps after inserted bases. Touching reads
// must subsequently be realigned, not merely assigned a shifted CIGAR.
size_t lift(size_t p, const std::vector<Accepted> &edits);
// Validate trace consumption/exact matches and calculate final-reference MD/NM.
std::string md(const Alignment &a, const std::string &ref, size_t &nm);
} // namespace hpc
