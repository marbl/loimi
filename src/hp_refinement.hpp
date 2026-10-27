// SPDX-License-Identifier: MIT-0
#pragma once
#include "core.hpp"

namespace hpc {
struct HpSite {
  size_t begin, end; // complete run in raw draft coordinates
  char base;
};
// One-pass local realignment of uncertain complete runs. Borrowed alignments
// have unique ranks [0,tiers), ambiguous ranks [tiers,2*tiers). Evidence is
// rebuilt, not added to the original histogram. No new model or retention.
// Fixed 64-base context, 12 exact terminal matches, pure run between immediate
// flanking bases. Targets near reference/read ends are skipped. Per-read traces
// and query-coordinate mapping are temporary; no per-read chromosome arrays.
// query_flank requires raw aligned query bases on both sides of
// every observation, preventing polishing from reintroducing end-filtered
// votes.
std::vector<Accepted> realign_homopolymers(const std::string &reference,
                                           const std::vector<const Alignment *> &alignments,
                                           const std::vector<HpSite> &sites, size_t tiers,
                                           uint32_t support, const std::map<int, RunModel> &models,
                                           bool pooled_fallback, size_t max_indel,
                                           size_t query_flank = 0);
} // namespace hpc
