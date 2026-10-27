// SPDX-License-Identifier: MIT-0
#include "hp_refinement.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace hpc {
std::vector<Accepted> realign_homopolymers(const std::string &ref,
                                           const std::vector<const Alignment *> &alns,
                                           const std::vector<HpSite> &sites, size_t tiers,
                                           uint32_t support, const std::map<int, RunModel> &models,
                                           bool pooled, size_t max_indel, size_t query_flank) {
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
} // namespace hpc
