// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstddef>
#include <string>
#include <utility>

namespace hpc::experimental {
// Diagnostic only: consecutive matching HPC operations immediately outside
// the focal run, capped to bound work per requested site. An indel or mismatch
// breaks the anchor. These are HPC runs, not raw bases or mapping uniqueness.
inline std::pair<std::size_t, std::size_t> matching_flanks(const std::string &ops, std::size_t pos,
                                                           std::size_t cap = 32) {
  if (pos >= ops.size() || ops[pos] != '=')
    return {0, 0};
  std::size_t left = 0, right = 0;
  while (left < cap && left < pos && ops[pos - left - 1] == '=')
    ++left;
  while (right < cap && right + 1 < ops.size() - pos && ops[pos + right + 1] == '=')
    ++right;
  return {left, right};
}
// Minimal first-stage filter, NOT a test of multi-base alignment quality.
// A length observation at an interior template run needs observed query runs
// on both sides. Natural/clipped read-terminal runs are otherwise censored.
// Preserve the legacy outward-facing template-end behavior; end extension and
// singleton passthrough have separate policies. All coordinates are HPC indexes
// in the oriented, clipped query and supplied template, respectively.
inline bool has_run_flanks(std::size_t query_run, std::size_t query_runs, std::size_t template_run,
                           std::size_t template_runs) {
  if (query_run >= query_runs || template_run >= template_runs)
    return false;
  return (query_run > 0 || template_run == 0) &&
         (query_run + 1 < query_runs || template_run + 1 == template_runs);
}

// Raw-base distances in the oriented, aligned (already clipped) query.
// Require adjacent matching HPC runs as local anchors, not a perfect long
// flank. Template endpoints waive only their outward-facing requirement.
// Boundary context counts only usable read sequence clipped by template
// alignment, never packaged/trimming clips. It supplies distance, not a waiver
// of the adjacent matching-run check. Non-boundary alignments pass zero context.
inline bool has_span_flanks(std::size_t begin, std::size_t end, std::size_t query_bases,
                            std::size_t minimum, bool left_match, bool right_match,
                            bool template_start, bool template_end, std::size_t left_context = 0,
                            std::size_t right_context = 0) {
  if (begin > end || end > query_bases)
    return false;
  // Subtraction avoids overflow when combining a long overhang with a span.
  bool left_span = begin >= minimum || left_context >= minimum - begin;
  auto right = query_bases - end;
  bool right_span = right >= minimum || right_context >= minimum - right;
  return (template_start || (left_match && left_span)) &&
         (template_end || (right_match && right_span));
}
} // namespace hpc::experimental
