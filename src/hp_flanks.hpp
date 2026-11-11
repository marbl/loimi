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
// context. near_template_* marks the bounded biological-end neighborhood in
// template HPC coordinates; unlike exact endpoints, these flags waive distance
// only, retaining the adjacent-match test and requiring the full inward flank.
inline bool has_span_flanks(std::size_t begin, std::size_t end, std::size_t query_bases,
                            std::size_t minimum, bool left_match, bool right_match,
                            bool template_start, bool template_end, std::size_t left_context = 0,
                            std::size_t right_context = 0, bool near_template_start = false,
                            bool near_template_end = false) {
  if (begin > end || end > query_bases)
    return false;
  // Subtraction avoids overflow when combining a long overhang with a span.
  bool left_span = begin >= minimum || left_context >= minimum - begin;
  auto right = query_bases - end;
  bool right_span = right >= minimum || right_context >= minimum - right;
  // A biological end cannot supply an arbitrarily long outward flank. Rescue
  // only that distance requirement, and only with a full inward flank. Use the
  // ORIGINAL spans on both sides so a short read cannot waive both requirements.
  // Adjacent HPC matches remain mandatory except at the exact template edge.
  bool terminal_left = near_template_start && right_span;
  bool terminal_right = near_template_end && left_span;
  return (template_start || (left_match && (left_span || terminal_left))) &&
         (template_end || (right_match && (right_span || terminal_right)));
}
} // namespace hpc
