// SPDX-License-Identifier: MIT-0
// Exact master algorithm as of 2026-10-23, isolated from the branch default.
#include "experiments/master_events.hpp"
#include <algorithm>
#include <stdexcept>
namespace hpc::experimental {
std::vector<Edit> master_events(const Alignment &a, const std::string &ref) {
  std::string expanded;
  if (a.packed_ops)
    expanded = expand_cigar(a.ops);
  const auto &ops = a.packed_ops ? expanded : a.ops;
  std::vector<Edit> out;
  size_t q = a.clip_left, r = a.begin;
  for (size_t i = 0; i < ops.size();) {
    if (ops[i] == '=') {
      ++q;
      ++r;
      ++i;
      continue;
    }
    Edit e{r, r, {}};
    while (i < ops.size() && ops[i] != '=') {
      char op = ops[i++];
      if (op != 'D')
        e.alt += a.bases().at(q++);
      if (op != 'I')
        ++r;
    }
    e.end = r;
    if (e.begin == e.end) {
      size_t shift = 0;
      while (e.begin > 0 && !e.alt.empty() && ref[e.begin - 1] == e.alt[e.alt.size() - 1 - shift]) {
        shift = (shift + 1) % e.alt.size();
        --e.begin;
        --e.end;
      }
      if (shift)
        std::rotate(e.alt.begin(), e.alt.end() - shift, e.alt.end());
    } else if (e.alt.empty()) {
      while (e.begin > 0 && ref[e.begin - 1] == ref[e.end - 1]) {
        --e.begin;
        --e.end;
      }
    }
    out.push_back(e);
  }
  if (q + a.clip_right != a.bases().size() || r != a.end)
    throw std::runtime_error("Invalid alignment trace");
  return out;
}
} // namespace hpc::experimental
