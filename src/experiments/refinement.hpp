// SPDX-License-Identifier: MIT-0
#pragma once
#include "core.hpp"
#include <optional>

namespace hpc::experimental {
struct UncertainCall {
  uint32_t first, second;
};
// Frozen corroboration rule: first tier only, unique depth qualified by caller,
// top-two log-likelihood gap <=3. Lower tiers keep their ordinary estimator.
std::optional<UncertainCall> corroboration_target(const RunModel &model, char base,
                                                  const Histogram &histogram);
// nullopt means insufficient model coverage/depth: try the next lower tier.
// Otherwise returns first/second and stops, never pooling LAO with raw UL.
std::optional<uint32_t> corroborate(const RunModel &model, char base, const Histogram &histogram,
                                    UncertainCall target);
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
// Optional corroboration recalculates the frozen first-tier unique decision on
// NEW histograms. query_flank requires raw aligned query bases on both sides of
// every observation, preventing polishing from reintroducing end-filtered votes.
std::vector<Accepted> realign_homopolymers(const std::string &reference,
                                           const std::vector<const Alignment *> &alignments,
                                           const std::vector<HpSite> &sites, size_t tiers,
                                           uint32_t support, const std::map<int, RunModel> &models,
                                           bool pooled_fallback, size_t max_indel,
                                           const std::map<int, RunModel> &corroboration_models = {},
                                           size_t query_flank = 0);
} // namespace hpc::experimental
