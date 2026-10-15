// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include <algorithm>
#include <set>
#include <stdexcept>

namespace hpc {
namespace {
// Align anchored prefixes, freeing only the excess of the longer sequence.
// Return query/reference consumption in their ORIGINAL roles.
struct Prefix {
  std::string ops;
  size_t query, reference;
};
Prefix prefix(const std::string &q, const std::string &r, double rate) {
  if (q.empty() || r.empty())
    throw std::runtime_error("Empty terminal alignment");
  if (q.size() >= r.size()) {
    auto a = align_boundary(q, r, false, true, rate);
    return {std::move(a.ops), a.query_end, r.size()};
  }
  auto a = align_boundary(r, q, false, true, rate);
  for (char &c : a.ops)
    if (c == 'I')
      c = 'D';
    else if (c == 'D')
      c = 'I';
  return {std::move(a.ops), q.size(), a.query_end};
}
struct Tail {
  Alignment *alignment;
  bool unique;
  std::string sequence;
};
std::string tail(const EndRead &r, bool left, size_t limit) {
  const auto &a = *r.alignment;
  size_t n = std::min(left ? r.left : r.right, limit);
  auto s = left ? a.bases().substr(a.clip_left - n, n)
                : a.bases().substr(a.bases().size() - a.clip_right, n);
  if (left)
    std::reverse(s.begin(), s.end());
  return s;
}
bool anchored(const Alignment &a, bool left, double rate) {
  auto ops = a.packed_ops ? expand_cigar(a.ops) : a.ops;
  size_t n = std::min<size_t>(100, ops.size());
  if (n < 20)
    return false;
  size_t errors = 0;
  for (size_t i = 0; i < n; ++i)
    errors += ops[left ? i : ops.size() - 1 - i] != '=';
  return errors <= rate * n;
}
std::string consensus(const std::vector<Tail> &tails, size_t nt,
                      const std::map<int, RunModel> &models, unsigned support, size_t limit,
                      double rate, bool ambiguous) {
  // Choose a path from the first adequately represented quality/placement
  // group. Never use an isolated longest/noisy read to define an end.
  const Tail *seed = nullptr;
  for (size_t phase = 0; phase < (ambiguous ? 2 : 1) && !seed; ++phase)
    for (size_t t = 0; t < nt && !seed; ++t) {
      size_t count = 0;
      const Tail *best = nullptr;
      for (const auto &x : tails)
        if (x.alignment->tier == int(t) && (phase || x.unique)) {
          ++count;
          if (!best || x.sequence.size() > best->sequence.size() ||
              (x.sequence.size() == best->sequence.size() && x.alignment->id < best->alignment->id))
            best = &x;
        }
      if (count >= support)
        seed = best;
    }
  if (!seed)
    return {};
  auto target = encode(seed->sequence);
  // Bounded by --max-end-extension and terminal coverage. Unlike chromosome
  // storage, this small local pileup can retain every tier simultaneously.
  std::vector<std::vector<Histogram>> votes(target.bases.size(), std::vector<Histogram>(2 * nt));
  for (const auto &x : tails) {
    if (!x.unique && !ambiguous)
      continue;
    auto q = encode(x.sequence);
    try {
      auto a = prefix(q.bases, target.bases, rate);
      auto mask = local_hpc_votes(a.ops, 100, rate);
      size_t qi = 0, ri = 0, opi = 0;
      for (char op : a.ops) {
        // A read's last run can be truncated by its end or the configured cap.
        if (op == '=' && qi + 1 < q.bases.size() && target.bases[ri] != 'N' &&
            (mask.empty() || mask[opi]))
          ++votes[ri][x.alignment->tier + (x.unique ? 0 : nt)][q.lengths[qi]];
        qi += op != 'D';
        ri += op != 'I';
        ++opi;
      }
    } catch (const std::runtime_error &) {
      // An unanchored or unresolved tail cannot establish extension coverage.
    }
  }
  std::string out;
  for (size_t i = 0; i < votes.size(); ++i) {
    Call chosen;
    Histogram histogram;
    for (size_t phase = 0; phase < (ambiguous ? 2 : 1) && chosen.tier < 0; ++phase)
      for (size_t t = 0; t < nt && chosen.tier < 0; ++t) {
        auto h = votes[i][t];
        if (phase)
          for (auto [length, n] : votes[i][nt + t])
            h[length] += n;
        auto c = call({h}, support);
        if (c.tier >= 0) {
          c.tier = t;
          chosen = c;
          histogram = std::move(h);
        }
      }
    if (chosen.tier < 0)
      break; // No gaps and no extrapolation on sparse evidence.
    auto model = models.find(chosen.tier);
    if (model != models.end())
      chosen.length = model->second.estimate(target.bases[i], histogram, chosen.length);
    if (chosen.length > limit - out.size())
      break;
    out.append(chosen.length, target.bases[i]);
  }
  return out;
}
} // namespace

EndExtension extend_ends(std::string &reference, const std::vector<EndRead> &reads, size_t tiers,
                         const std::map<int, RunModel> &models, unsigned support, size_t max_bases,
                         double rate, bool ambiguous, bool do_left, bool do_right) {
  if (!support || !max_bases || !tiers)
    throw std::runtime_error("Invalid end-extension settings");
  std::string extensions[2];
  EndExtension result;
  for (bool left : {true, false}) {
    if (left ? !do_left : !do_right)
      continue;
    std::vector<Tail> tails;
    std::set<uint32_t> seen;
    for (const auto &r : reads) {
      auto &a = *r.alignment;
      if ((left ? a.begin != 0 || !r.left : a.end != reference.size() || !r.right) || !r.eligible ||
          (!r.unique && !ambiguous) || !anchored(a, left, rate))
        continue;
      // Alternative placements of the same molecule are not independent depth.
      if (!seen.insert(a.id).second)
        continue;
      tails.push_back({&a, r.unique, tail(r, left, max_bases)});
    }
    (left ? result.left_candidates : result.right_candidates) = tails.size();
    extensions[left ? 0 : 1] = consensus(tails, tiers, models, support, max_bases, rate, ambiguous);
  }
  size_t old_size = reference.size(), shift = extensions[0].size();
  if (!shift && extensions[1].empty())
    return result; // Avoid copying a whole chromosome when neither end can extend.
  for (const auto &r : reads) {
    auto &a = *r.alignment;
    size_t old_begin = a.begin, old_end = a.end;
    a.begin += shift;
    a.end += shift;
    for (bool left : {true, false}) {
      const auto &extension = extensions[left ? 0 : 1];
      if (extension.empty() || (left ? old_begin != 0 || !r.left : old_end != old_size || !r.right))
        continue;
      try {
        auto x = prefix(tail(r, left, max_bases), extension, rate);
        auto ops = a.packed_ops ? expand_cigar(a.ops) : a.ops;
        if (left) {
          std::reverse(x.ops.begin(), x.ops.end());
          ops = x.ops + ops;
          a.begin -= x.reference;
          a.clip_left -= x.query;
        } else {
          ops += x.ops;
          a.end += x.reference;
          a.clip_right -= x.query;
        }
        a.ops = cigar(ops);
        a.packed_ops = true;
      } catch (const std::runtime_error &) {
        // Keep unresolved overhangs soft-clipped; never drop the interior read.
      }
    }
  }
  std::reverse(extensions[0].begin(), extensions[0].end());
  reference.reserve(old_size + shift + extensions[1].size());
  reference.insert(0, extensions[0]);
  reference += extensions[1];
  result.left = shift;
  result.right = extensions[1].size();
  return result;
}
} // namespace hpc
