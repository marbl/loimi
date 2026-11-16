// SPDX-License-Identifier: MIT-0
#pragma once
#include "core.hpp"
#include <functional>
namespace hpc {
// All intervals are zero-based half-open HPC run indexes. Raw strings are
// temporary and confined to a single candidate window.
struct RepeatConfig {
  size_t flank_bases = 1000;
  size_t max_window_bases = 6000;
  size_t max_reads = 256;
  size_t max_candidates = 8;
  size_t projection_cells = 16000000;
  size_t window_projection_cells = 64000000;
  uint32_t min_support = 5;
  bool ambiguous_fallback = true;
};
struct RepeatRead {
  const Alignment *alignment = nullptr; // Borrow immutable, reference-oriented bases.
  size_t begin = 0, end = 0;
  bool unique = false;
  unsigned placements = 0;
};
struct RepeatInterval {
  size_t begin = 0, end = 0;
};
struct RepeatUpdate {
  size_t run;
  bool model_changed, mixed;
};
struct RepeatDecision {
  RepeatInterval core, window;
  std::string status;
  size_t reads = 0, candidates = 0, changed_runs = 0;
  int tier = -1;
  bool pooled = false;
  int64_t before_score = 0, after_score = 0;
  std::vector<RepeatUpdate> updates;
};
// Linear scan for low-complexity or short-period HPC repeats, gated by discordant calls.
// Returned intervals are disjoint; no known-error coordinates are accepted.
std::vector<RepeatInterval> detect_repeats(const std::string &backbone,
                                           const std::vector<Call> &calls);
// Exact minimum-edit projection onto B1+ B2+ ... Bm+. Throws length_error at
// the explicit scratch bound; caller retains its original call on exhaustion.
std::string project_repeat_backbone(const std::string &sequence, const std::string &backbone,
                                    size_t max_cells);
// Process windows serially with bounded scratch; score all admitted molecules.
// Oversized/deep/unanchorable intervals are reported, never silently sampled.
// Existing call metadata is preserved unless a qualified new histogram wins.
std::vector<RepeatInterval>
refine_repeats(const std::string &backbone, std::vector<Call> &calls, std::vector<RepeatRead> reads,
               size_t tiers, const std::map<int, RunModel> &models, const RepeatConfig &config,
               const std::function<void(const RepeatDecision &)> &report);
} // namespace hpc
