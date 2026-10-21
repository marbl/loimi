// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace hpc {
// Failure-only rescue anchors in HPC coordinates. Each exact chunk must occur
// once in the template; different chunks from one read are NOT independent votes.
std::vector<int64_t> internal_offsets(const std::string &query, const std::string &reference,
                                      int64_t layout_begin);
// Require >=75% of informative reads and reject a strong distant competitor.
// Empty/ambiguous evidence throws. Returns one-read-one-vote lower median.
int64_t internal_offset_cluster(const std::vector<std::vector<int64_t>> &reads, size_t required,
                                int64_t tolerance = 20);

// Offset discovery is separate from consensus eligibility. Each record is a
// known layout placement; repeated IDs never count as independent supporters.
struct OffsetEvidence {
  uint32_t id = 0;
  bool unique = false;
  std::vector<int64_t> offsets;
};
struct OffsetChoice {
  int64_t value = 0;
  size_t support = 0, informative = 0, runner_up = 0;
};
// Stable bounded sampling shared by whole-tig and region planning.
uint64_t offset_sample_key(uint32_t id);
// Whole shorter-sequence alignment first; tied alignments remain unresolved.
// Failed partial overlaps may use collinear, separated exact chunks, at most
// 64 sampled positions. Coordinates and distances are HPC bases.
std::vector<int64_t> discover_offsets(const std::string &query, const std::string &reference,
                                      int64_t layout_begin, double max_error);
// >=60% informative molecules, strictly >2x the distant runner-up, and the
// requested minimum independent support. Strong unique conflicts veto rescue.
OffsetChoice choose_offset(const std::vector<OffsetEvidence> &reads, size_t required,
                           bool unique_only = false);
bool validate_offset_overlap(const std::string &query, const std::string &reference,
                             int64_t layout_begin, int64_t layout_end, int64_t offset,
                             double max_error);
} // namespace hpc
