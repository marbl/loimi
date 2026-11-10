// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace hpc {
inline constexpr int64_t offset_local_span = 256;
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
  int tier = -1;            // first decisive tier; lower tiers cannot override it
  size_t repeat_period = 0; // nonzero only for verified periodic-template tie rescue
};
// Failure-only joint inference. Scores refer to a fixed predicted query window;
// an alignment cannot search another repeat copy. Invalid/nonoverlapping scores
// have overlap == 0. All coordinates and lengths are HPC bases.
struct OffsetScore {
  size_t overlap = 0, edits = 0;
};
struct OffsetRow {
  std::string molecule;
  int tier = 0;
  bool unique = false;
  std::vector<OffsetScore> scores;
};
std::string offset_molecule(const std::string &name, uint32_t id);
// Coalesce local coordinate uncertainty (256 HPC bases), never distant copies.
// Each locus is represented by the median proposals from its highest tier.
std::vector<int64_t> offset_candidates(const std::vector<OffsetEvidence> &reads,
                                       const std::vector<int> &tiers);
OffsetScore score_offset_overlap(const std::string &query, const std::string &reference,
                                 int64_t layout_begin, int64_t layout_end, int64_t offset,
                                 double max_error);
// Failure-only rescore allowing <=256 HPC bases of endpoint drift at the same
// candidate locus. Scores a trimmed reference core; requires >=1000 original
// overlap bases. Never searches outside the layout-predicted query overlap.
OffsetScore score_offset_overlap_local(const std::string &query, const std::string &reference,
                                       int64_t layout_begin, int64_t layout_end, int64_t offset,
                                       double max_error);
// Lower tiers corroborate higher-tier candidates. One molecule votes only in
// its highest usable tier; repeated placements and corrected/raw UUID aliases
// cannot increase independent support. Inconclusive tiers fall through; throws
// if no tier has sufficient corroboration, majority support and a resolved locus.
OffsetChoice choose_joint_offset(const std::vector<int64_t> &candidates,
                                 const std::vector<OffsetRow> &rows,
                                 const std::string &reference = {});
// Verify an exact period spanning at least two copies over the ENTIRE template.
// Bounded candidate search; returns zero when no period is proven.
size_t template_repeat_period(const std::string &reference);
struct BoundaryAlignment;
// Restrict boundary alignment to a known layout locus, with endpoint movement
// bounded by max_shift (0..256 HPC bases). Used after local-drift rescue.
BoundaryAlignment align_layout_boundary(const std::string &query, const std::string &reference,
                                        bool left, bool right, size_t predicted_begin,
                                        size_t max_shift, double max_error);
// Align an already placed repeat read near its predicted query start. The
// adjustment is capped below half a repeat to prevent jumps to adjacent copies.
BoundaryAlignment align_repeat_boundary(const std::string &query, const std::string &reference,
                                        bool left, bool right, size_t predicted_begin,
                                        size_t period, double max_error);
// Stable bounded sampling shared by whole-tig and region planning.
uint64_t offset_sample_key(uint32_t id);
// Whole shorter-sequence alignment first; distant tied loci remain unresolved.
// Optimal endpoints within a 256-base neighborhood are one uncertain locus.
// Failed partial overlaps may use collinear, separated exact chunks, at most
// 64 sampled positions. Coordinates and distances are HPC bases.
std::vector<int64_t> discover_offsets(const std::string &query, const std::string &reference,
                                      int64_t layout_begin, double max_error);
bool validate_offset_overlap(const std::string &query, const std::string &reference,
                             int64_t layout_begin, int64_t layout_end, int64_t offset,
                             double max_error);
} // namespace hpc
