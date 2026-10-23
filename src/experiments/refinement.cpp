// SPDX-License-Identifier: MIT-0
#include "experiments/refinement.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace hpc::experimental {
std::optional<UncertainCall> corroboration_target(const RunModel &model, char base,
                                                  const Histogram &h) {
  auto scores = model.likelihood_scores(base, h);
  if (scores.size() < 2)
    return {};
  std::vector<size_t> order(scores.size());
  std::iota(order.begin(), order.end(), 0);
  std::partial_sort(order.begin(), order.begin() + 2, order.end(), [&](size_t a, size_t b) {
    return scores[a] != scores[b] ? scores[a] > scores[b] : a < b;
  });
  if (scores[order[0]] - scores[order[1]] > 3)
    return {};
  return UncertainCall{uint32_t(order[0] + 1), uint32_t(order[1] + 1)};
}
std::optional<uint32_t> corroborate(const RunModel &model, char base, const Histogram &h,
                                    UncertainCall target) {
  auto scores = model.likelihood_scores(base, h, false);
  auto depth = call({h}, 1).depth;
  if (depth < 5 || target.first > scores.size() || target.second > scores.size())
    return {};
  auto found = h.find(target.second);
  double fraction = found == h.end() ? 0 : double(found->second) / depth;
  return (scores[target.second - 1] - scores[target.first - 1]) / depth >= .5 && fraction >= .6
             ? target.second
             : target.first;
}

std::vector<Accepted>
realign_homopolymers(const std::string &ref, const std::vector<const Alignment *> &alns,
                     const std::vector<HpSite> &sites, size_t tiers, uint32_t support,
                     const std::map<int, RunModel> &models, bool pooled, size_t max_indel,
                     const std::map<int, RunModel> &corroboration_models, size_t query_flank) {
  constexpr size_t flank = 64, anchor = 12;
  // Sparse histogram maps allocate only at candidate sites, not every base.
  std::vector<std::map<int, Histogram>> evidence(sites.size());
  for (const auto *a : alns) {
    if (a->tier < 0 || size_t(a->tier) >= 2 * tiers)
      throw std::runtime_error("Invalid HP realignment tier");
    auto first = std::lower_bound(sites.begin(), sites.end(), a->begin + flank,
                                  [](const HpSite &s, size_t p) { return s.begin < p; });
    if (first == sites.end() || first->end + flank > a->end)
      continue;
    // Map each reference boundary to its query boundary once per read. Gap
    // placement at an outer boundary is harmless only if exact anchors pass.
    std::vector<size_t> query(a->end - a->begin + 1, SIZE_MAX);
    auto ops = a->packed_ops ? expand_cigar(a->ops) : a->ops;
    size_t q = a->clip_left, r = a->begin;
    for (char op : ops) {
      query.at(r - a->begin) = q;
      q += op != 'D';
      r += op != 'I';
    }
    query.at(r - a->begin) = q;
    for (auto site = first; site != sites.end() && site->end + flank <= a->end; ++site) {
      if (site->begin < flank || site->end + flank > ref.size() || site->end - site->begin > 1000)
        continue;
      size_t left = site->begin - flank, right = site->end + flank;
      size_t qb = query.at(left - a->begin), qe = query.at(right - a->begin);
      if (qe <= qb || qe - qb > 4096)
        continue;
      std::string sequence = a->bases().substr(qb, qe - qb), trace;
      try {
        trace = align(sequence, ref.substr(left, right - left), .15);
      } catch (const std::bad_alloc &) {
        throw;
      } catch (const std::runtime_error &) {
        continue;
      }
      if (trace.size() < 2 * anchor || trace.substr(0, anchor) != std::string(anchor, '=') ||
          trace.substr(trace.size() - anchor) != std::string(anchor, '='))
        continue;
      // Locate the two non-run bases bounding the complete allele. Observations
      // with substitutions/compound sequence between those anchors cannot vote.
      q = 0;
      r = left;
      size_t start = SIZE_MAX, end = SIZE_MAX;
      for (char op : trace) {
        if (op == '=' && r == site->begin - 1)
          start = q + 1;
        if (op == '=' && r == site->end)
          end = q;
        q += op != 'D';
        r += op != 'I';
      }
      if (start == SIZE_MAX || end == SIZE_MAX || end <= start ||
          !std::all_of(sequence.begin() + start, sequence.begin() + end,
                       [&](char b) { return b == site->base; }))
        continue;
      if (qb + start - a->clip_left < query_flank ||
          a->bases().size() - a->clip_right - (qb + end) < query_flank)
        continue;
      ++evidence[site - sites.begin()][a->tier][uint32_t(end - start)];
    }
  }
  std::vector<Accepted> edits;
  for (size_t i = 0; i < sites.size(); ++i) {
    const auto &s = sites[i];
    bool decided = false;
    for (size_t phase = 0; phase < (pooled ? 2u : 1u) && !decided; ++phase)
      for (size_t tier = 0; tier < tiers; ++tier) {
        Histogram h = evidence[i][int(tier)];
        if (phase)
          for (auto [length, count] : evidence[i][int(tier + tiers)])
            h[length] += count;
        auto c = call({h}, support);
        if (c.depth < support)
          continue;
        decided = true; // Never skip a qualifying tier because its allele loses.
        auto model = models.find(int(tier));
        if (model != models.end())
          c.length = model->second.estimate(s.base, h, c.length);
        // Sequential composition, not double counting: the first-stage vote
        // histogram is discarded. Only rebuilt UNIQUE evidence may corroborate.
        if (phase == 0 && tier == 0 && model != models.end() && !corroboration_models.empty()) {
          auto target = corroboration_target(model->second, s.base, h);
          if (target)
            for (size_t lower = 1; lower < tiers; ++lower) {
              auto calibration = corroboration_models.find(int(lower));
              if (calibration == corroboration_models.end())
                continue;
              auto result =
                  corroborate(calibration->second, s.base, evidence[i][int(lower)], *target);
              if (result) {
                c.length = *result;
                break;
              }
            }
        }
        size_t old = s.end - s.begin;
        if (c.length != old &&
            std::max(size_t(c.length), old) - std::min(size_t(c.length), old) <= max_indel)
          edits.push_back({{s.begin, s.end, std::string(c.length, s.base)},
                           int(tier + phase * tiers),
                           c.depth,
                           c.peak});
        break;
      }
  }
  return edits;
}
} // namespace hpc::experimental
