// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstddef>

namespace hpc {
// Raw-base distances in the oriented, aligned (already clipped) query.
// Require adjacent matching HPC runs as local anchors, not a perfect long
// flank. Template endpoints waive only their outward-facing requirement.
// Boundary context counts only usable read sequence clipped by template
// alignment, never packaged/trimming clips. It supplies distance, not a waiver
// of the adjacent matching-run check. Non-boundary alignments pass zero
// context.
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
} // namespace hpc
