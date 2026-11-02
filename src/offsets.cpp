// SPDX-License-Identifier: MIT-0
#include "offsets.hpp"
#include "align/edlib.H"
#include "core.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace hpc {
namespace {
constexpr double missing_overlap_penalty = .01;
constexpr double indistinguishable_score_slack = .001;
constexpr size_t max_offset_loci = 128;
} // namespace
std::vector<int64_t> internal_offsets(const std::string &q, const std::string &r, int64_t begin) {
  const size_t width = std::min<size_t>(5000, r.size() / 4);
  std::vector<int64_t> result;
  if (width < 128 || q.size() < width)
    return result;
  for (size_t i = 0; i < 5; ++i) {
    size_t qb = (q.size() - width) * i / 4;
    auto chunk = q.substr(qb, width);
    size_t p = r.find(chunk);
    if (p != std::string::npos && r.find(chunk, p + 1) == std::string::npos)
      result.push_back(begin + int64_t(qb) - int64_t(p));
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

int64_t internal_offset_cluster(const std::vector<std::vector<int64_t>> &reads, size_t required,
                                int64_t tolerance) {
  std::vector<int64_t> starts;
  size_t informative = 0;
  for (const auto &r : reads) {
    informative += !r.empty();
    starts.insert(starts.end(), r.begin(), r.end());
  }
  std::sort(starts.begin(), starts.end());
  auto votes = [&](int64_t lo) {
    std::vector<int64_t> v;
    for (const auto &r : reads) {
      std::vector<int64_t> within;
      for (auto x : r)
        if (x >= lo && x - lo <= tolerance)
          within.push_back(x);
      if (!within.empty()) {
        std::sort(within.begin(), within.end());
        v.push_back(within[(within.size() - 1) / 2]);
      }
    }
    std::sort(v.begin(), v.end());
    return v;
  };
  std::vector<int64_t> best;
  int64_t start = 0;
  for (auto x : starts) {
    auto v = votes(x);
    if (v.size() > best.size()) {
      best = std::move(v);
      start = x;
    }
  }
  if (best.empty() || best.size() < required || best.size() * 4 < informative * 3)
    throw std::runtime_error("Insufficient independent internal offset anchors");
  for (auto x : starts)
    if (std::abs(x - start) > 2 * tolerance && votes(x).size() * 2 >= best.size())
      throw std::runtime_error("Competing internal offset anchors");
  return best[(best.size() - 1) / 2];
}
uint64_t offset_sample_key(uint32_t id) {
  uint64_t x = uint64_t(id) + 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

std::vector<int64_t> discover_offsets(const std::string &q, const std::string &r, int64_t begin,
                                      double rate) {
  // This rescue establishes an origin from sequence, not just layout geometry.
  // Tiny matches cannot provide a reliable anchor even when many reads agree.
  if (std::min(q.size(), r.size()) < 128 || q.size() > INT32_MAX || r.size() > INT32_MAX)
    return {};
  const bool enclosing = q.size() > r.size();
  const auto &small = enclosing ? r : q;
  const auto &large = enclosing ? q : r;
  auto a = merylutil::edlibAlign(small.data(), small.size(), large.data(), large.size(),
                                 merylutil::edlibNewAlignConfig(int(std::ceil(rate * small.size())),
                                                                merylutil::EDLIB_MODE_HW,
                                                                merylutil::EDLIB_TASK_LOC));
  // Do not turn tied whole-overlap evidence into a confident isolated seed.
  bool resolved = a.editDistance >= 0 && a.numLocations > 0 && a.startLocations;
  // Nearby optimal endpoints describe local alignment uncertainty, not a
  // different repeat copy. Preserve only a single local locus here.
  bool tied = false;
  if (resolved)
    for (int i = 1; i < a.numLocations; ++i)
      tied |= std::abs(a.startLocations[i] - a.startLocations[0]) > offset_local_span ||
              std::abs(a.endLocations[i] - a.endLocations[0]) > offset_local_span;
  resolved &= !tied;
  int64_t p = resolved ? a.startLocations[0] : -1;
  // A noisy reverse alignment pinned to a read boundary can be a truncated
  // template, rather than a truly enclosing read. Let internal chains test it.
  bool truncated =
      resolved && a.editDistance > 0 && (p == 0 || a.endLocations[0] + 1 == int(large.size()));
  merylutil::edlibFreeAlignResult(a);
  if (resolved && !truncated)
    return {enclosing ? begin + p : begin - p};
  if (tied)
    return {};

  // Sparse partial overlaps can contain several good shorter anchors but no
  // two 2-kb anchors separated by 2 kb. Try smaller exact chunks only when
  // the stronger chain fails. Every proposal still needs fixed-overlap scoring.
  for (size_t width = std::min<size_t>(2000, r.size() / 4); width >= 128; width /= 2) {
    if (q.size() < 2 * width)
      continue;
    size_t samples = std::min<size_t>(64, 1 + (q.size() - width + width / 2 - 1) / (width / 2));
    std::vector<std::pair<int64_t, size_t>> points;
    for (size_t i = 0; i < samples; ++i) {
      size_t qb = (q.size() - width) * i / (samples - 1);
      auto chunk = q.substr(qb, width);
      if (chunk.find('N') != std::string::npos)
        continue;
      size_t hit = r.find(chunk);
      if (hit != std::string::npos && r.find(chunk, hit + 1) == std::string::npos)
        points.emplace_back(begin + int64_t(qb) - int64_t(hit), qb);
    }
    std::sort(points.begin(), points.end());
    std::vector<int64_t> result;
    for (size_t i = 0; i < points.size(); ++i) {
      size_t j = i, lo = points[i].second, hi = lo;
      while (j < points.size() && points[j].first - points[i].first <= 20) {
        lo = std::min(lo, points[j].second);
        hi = std::max(hi, points[j].second);
        ++j;
      }
      if (hi - lo >= width)
        result.push_back(points[i + (j - i - 1) / 2].first);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    if (!result.empty())
      return result;
  }
  return {};
}

OffsetScore score_offset_overlap(const std::string &q, const std::string &r, int64_t begin,
                                 int64_t end, int64_t offset, double rate) {
  auto b = begin - offset, e = end - offset;
  auto lo = std::max<int64_t>(0, b), hi = std::min<int64_t>(r.size(), e);
  if (hi <= lo || hi - lo < 128 || hi - lo < int64_t(std::min<size_t>(1000, r.size())))
    return {};
  // Clip only the overhang predicted by the layout transform. Using an infix
  // alignment of the entire read here would silently relocate a candidate to
  // whichever repeat copy aligns best, invalidating the coordinate check.
  int64_t qb = std::max<int64_t>(0, -b);
  int64_t qe = e <= int64_t(r.size()) ? int64_t(q.size())
                                      : std::min<int64_t>(q.size(), int64_t(r.size()) - b);
  if (qb >= qe || qe - qb > INT32_MAX || hi - lo > INT32_MAX)
    return {};
  auto a = merylutil::edlibAlign(q.data() + qb, int(qe - qb), r.data() + lo, int(hi - lo),
                                 merylutil::edlibNewAlignConfig(int(std::ceil(rate * (hi - lo))),
                                                                merylutil::EDLIB_MODE_NW,
                                                                merylutil::EDLIB_TASK_DISTANCE));
  int edits = a.editDistance;
  merylutil::edlibFreeAlignResult(a);
  return edits < 0 ? OffsetScore{} : OffsetScore{size_t(hi - lo), size_t(edits)};
}

bool validate_offset_overlap(const std::string &q, const std::string &r, int64_t begin, int64_t end,
                             int64_t offset, double rate) {
  return score_offset_overlap(q, r, begin, end, offset, rate).overlap != 0;
}

std::string offset_molecule(const std::string &name, uint32_t id) {
  // Corrected ONT and its original UL record can have different packaged IDs.
  // Recognize only the explicit hifi_UUID / ont_UUID convention; do not guess
  // identities from arbitrary names or remove PacBio subread/CCS suffixes.
  size_t prefix = name.starts_with("hifi_") ? 5 : name.starts_with("ont_") ? 4 : 0;
  if (prefix && name.size() == prefix + 36) {
    auto uuid = name.substr(prefix);
    bool valid = true;
    for (size_t i = 0; i < uuid.size(); ++i) {
      bool dash = i == 8 || i == 13 || i == 18 || i == 23;
      valid &= dash ? uuid[i] == '-' : bool(std::isxdigit(static_cast<unsigned char>(uuid[i])));
      uuid[i] = std::tolower(static_cast<unsigned char>(uuid[i]));
    }
    if (valid)
      return "uuid:" + uuid;
  }
  return "id:" + std::to_string(id);
}

std::vector<int64_t> offset_candidates(const std::vector<OffsetEvidence> &reads,
                                       const std::vector<int> &tiers) {
  if (reads.size() != tiers.size())
    throw std::invalid_argument("Offset tier/evidence size mismatch");
  struct Proposal {
    int64_t value;
    int tier;
    uint32_t id;
  };
  std::vector<Proposal> proposals;
  for (size_t i = 0; i < reads.size(); ++i)
    for (auto x : reads[i].offsets)
      proposals.push_back({x, tiers[i], reads[i].id});
  std::sort(proposals.begin(), proposals.end(), [](auto a, auto b) { return a.value < b.value; });
  std::vector<int64_t> result;
  for (size_t i = 0; i < proposals.size();) {
    size_t end = i;
    int tier = proposals[i].tier;
    while (end < proposals.size() &&
           proposals[end].value - proposals[i].value <= offset_local_span) {
      tier = std::min(tier, proposals[end].tier);
      ++end;
    }
    std::vector<int64_t> preferred;
    std::map<uint32_t, std::vector<int64_t>> by_id;
    for (size_t j = i; j < end; ++j)
      if (proposals[j].tier == tier)
        by_id[proposals[j].id].push_back(proposals[j].value);
    for (auto &[id, values] : by_id) {
      values.erase(std::unique(values.begin(), values.end()), values.end());
      preferred.push_back(values[(values.size() - 1) / 2]);
    }
    std::sort(preferred.begin(), preferred.end());
    result.push_back(preferred[(preferred.size() - 1) / 2]);
    i = end;
  }
  // This is a failure-only bounded diagnostic, not an unbounded all-pairs
  // mapper. Excessive locus ambiguity must fail explicitly instead of dropping
  // an arbitrary competitor or growing chromosome-scale alignment work.
  if (result.size() > max_offset_loci)
    throw std::runtime_error("Offset rescue has more than 128 distinct candidate loci");
  return result;
}

size_t template_repeat_period(const std::string &reference) {
  if (reference.size() < 2 || reference.find('N') != std::string::npos)
    return 0;
  auto anchor = reference.substr(0, std::min<size_t>(256, reference.size() / 2));
  size_t pos = 0;
  for (size_t attempt = 0; attempt < 128; ++attempt) {
    pos = reference.find(anchor, pos + 1);
    if (pos == std::string::npos || pos > reference.size() / 2)
      break;
    if (std::equal(reference.begin() + pos, reference.end(), reference.begin()))
      return pos;
  }
  return 0;
}

BoundaryAlignment align_repeat_boundary(const std::string &q, const std::string &r, bool left,
                                        bool right, size_t predicted_begin, size_t period,
                                        double rate) {
  if (q.empty() || r.empty() || (!left && !right) || period < 2 || predicted_begin >= q.size() ||
      q.size() > INT32_MAX || r.size() > INT32_MAX)
    throw std::runtime_error("Invalid repeat boundary alignment");
  size_t radius = std::min<size_t>(offset_local_span, (period - 1) / 2);
  size_t qb = left ? predicted_begin : 0;
  size_t qe = right ? std::min(q.size(), qb + r.size()) : q.size();
  size_t lo = left && qb > radius ? qb - radius : 0;
  size_t hi = right ? std::min(q.size(), qe + radius) : q.size();
  auto query = q.substr(lo, hi - lo);
  auto target = r;
  bool reverse = left && !right;
  if (reverse) {
    std::reverse(query.begin(), query.end());
    std::reverse(target.begin(), target.end());
  }
  auto a = merylutil::edlibAlign(target.data(), target.size(), query.data(), query.size(),
                                 merylutil::edlibNewAlignConfig(int(std::ceil(rate * r.size())),
                                                                left && right
                                                                    ? merylutil::EDLIB_MODE_HW
                                                                    : merylutil::EDLIB_MODE_SHW,
                                                                merylutil::EDLIB_TASK_LOC));
  size_t best = SIZE_MAX, begin = 0, end = 0;
  for (int i = 0; a.editDistance >= 0 && i < a.numLocations; ++i) {
    size_t b = a.startLocations[i], e = a.endLocations[i] + 1;
    if (reverse) {
      auto old = b;
      b = query.size() - e;
      e = query.size() - old;
    }
    b += lo;
    e += lo;
    size_t db = b > qb ? b - qb : qb - b, de = e > qe ? e - qe : qe - e;
    if (db > (left ? radius : 0) || de > (right ? radius : 0) || e <= b)
      continue;
    if (db + de < best || (db + de == best && b < begin)) {
      best = db + de;
      begin = b;
      end = e;
    }
  }
  merylutil::edlibFreeAlignResult(a);
  if (best == SIZE_MAX)
    throw std::runtime_error("No repeat boundary alignment near layout prediction");
  return {align(q.substr(begin, end - begin), r, rate), begin, end};
}

OffsetChoice choose_joint_offset(const std::vector<int64_t> &candidates,
                                 const std::vector<OffsetRow> &rows, const std::string &reference) {
  const size_t n = candidates.size();
  if (!n || rows.empty())
    throw std::runtime_error("No joint offset candidates");
  std::map<std::string, std::vector<const OffsetRow *>> molecules;
  std::set<int> tiers;
  for (const auto &row : rows) {
    if (row.scores.size() != n || row.tier < 0)
      throw std::invalid_argument("Invalid joint offset score row");
    molecules[row.molecule].push_back(&row);
    tiers.insert(row.tier);
  }
  std::map<int, std::vector<size_t>> votes, unique_votes;
  std::map<int, size_t> informative;
  for (auto tier : tiers) {
    votes[tier].resize(n);
    unique_votes[tier].resize(n);
  }
  for (const auto &[molecule, placements] : molecules) {
    int tier = std::numeric_limits<int>::max();
    for (auto row : placements)
      if (std::any_of(row->scores.begin(), row->scores.end(),
                      [](auto s) { return s.overlap != 0; }))
        tier = std::min(tier, row->tier);
    if (tier == std::numeric_limits<int>::max())
      continue;
    size_t longest = 0;
    for (auto row : placements)
      if (row->tier == tier)
        for (auto s : row->scores)
          longest = std::max(longest, s.overlap);
    std::vector<double> cost(n, std::numeric_limits<double>::infinity());
    std::vector<bool> unique(n, false);
    for (auto row : placements) {
      if (row->tier != tier)
        continue;
      for (size_t i = 0; i < n; ++i) {
        auto s = row->scores[i];
        if (!s.overlap)
          continue;
        // Do not reward an incorrect origin merely for leaving a tiny perfect
        // overlap. A modest coverage penalty still permits genuine overhangs.
        double c = double(s.edits) / s.overlap +
                   missing_overlap_penalty * (1. - double(s.overlap) / longest);
        if (c < cost[i]) {
          cost[i] = c;
          unique[i] = row->unique;
        }
      }
    }
    double best = *std::min_element(cost.begin(), cost.end());
    ++informative[tier];
    for (size_t i = 0; i < n; ++i)
      // Near-equal alignments support both hypotheses, never an arbitrary
      // winner.
      if (cost[i] <= best + indistinguishable_score_slack) {
        ++votes[tier][i];
        unique_votes[tier][i] += unique[i];
      }
  }
  std::vector<size_t> support(n);
  for (auto tier : tiers)
    for (size_t i = 0; i < n; ++i)
      support[i] += votes[tier][i];
  // A three-molecule sample may contain one outlier, just as four/five-read
  // samples already can. Keep the existing minimum unchanged for every other
  // sample size; tier agreement and competing-candidate checks still apply.
  size_t required = std::min<size_t>(3, molecules.size() / 2 + 1);
  size_t total = 0;
  for (auto [tier, count] : informative)
    total += count;
  std::ostringstream rejected;
  // Decide within one tier at a time. A higher tier with split, tied, or
  // insufficiently corroborated evidence cannot veto a decisive lower tier.
  // Cross-tier support is used only for the independent-molecule minimum:
  // lower-tier abundance never overrides a qualifying higher-tier decision.
  // First preserve all evidence-based decisions, including a decisive later
  // tier. Only if every tier is unresolved may repeat equivalence break a tie.
  size_t period = 0;
  for (unsigned rescue = 0; rescue < 2; ++rescue) {
    if (rescue) {
      period = template_repeat_period(reference);
      if (!period)
        break;
    }
    for (auto tier : tiers) {
      if (!informative[tier])
        continue;
      auto compare = [&](size_t a, size_t b) {
        if (votes[tier][a] != votes[tier][b])
          return votes[tier][a] > votes[tier][b] ? 1 : -1;
        if (unique_votes[tier][a] != unique_votes[tier][b])
          return unique_votes[tier][a] > unique_votes[tier][b] ? 1 : -1;
        return 0;
      };
      size_t best = 0;
      for (size_t i = 1; i < n; ++i)
        if (compare(i, best) > 0 ||
            (rescue && compare(i, best) == 0 && candidates[i] < candidates[best]))
          best = i;
      size_t runner = 0;
      bool tied = false;
      for (size_t i = 0; i < n; ++i)
        if (i != best && std::abs(candidates[i] - candidates[best]) > 2 * offset_local_span) {
          runner = std::max(runner, support[i]);
          if (compare(i, best) == 0) {
            size_t distance = std::abs(candidates[i] - candidates[best]);
            size_t phase = period ? distance % period : 0;
            // Allow the same small coordinate drift as internal offset clusters,
            // not arbitrary phase shifts or partially repetitive templates.
            bool equivalent = rescue && distance >= period / 2 &&
                              std::min(phase, period - phase) <= std::min<size_t>(20, period / 4);
            tied |= !equivalent;
          }
        }
      bool majority = votes[tier][best] * 5 >= informative[tier] * 3;
      // Sparse tiers can be decisive even when other sampled molecules provide
      // no usable alignment. Keep the distant-candidate ambiguity veto: one read
      // supporting two repeat copies is not a resolved unanimous vote.
      bool unanimous = votes[tier][best] == informative[tier];
      if ((support[best] >= required || unanimous) && majority && !tied)
        return {candidates[best], support[best], total, runner, tier, period};
      rejected << " tier=" << tier << " candidate=" << candidates[best]
               << " tier_support=" << votes[tier][best] << "/" << informative[tier]
               << " support=" << support[best] << " required=" << required
               << " runner_up=" << runner << " tied=" << tied << " majority=" << majority << ";";
    }
  }
  throw std::runtime_error("Joint offset unresolved:" + rejected.str());
}
} // namespace hpc
