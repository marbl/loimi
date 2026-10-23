// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include "experiments/master_events.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace hpc {
namespace {
// Canonicalize a complete read-supported replacement, not independent pieces
// of a haplotype. Trim shared sequence before left-normalizing a pure indel.
Edit normalize(Edit e, const std::string &ref) {
  while (e.end > e.begin && !e.alt.empty() && ref[e.end - 1] == e.alt.back()) {
    --e.end;
    e.alt.pop_back();
  }
  size_t shared = 0;
  while (e.begin + shared < e.end && shared < e.alt.size() &&
         ref[e.begin + shared] == e.alt[shared])
    ++shared;
  e.begin += shared;
  e.alt.erase(0, shared);
  if (e.begin == e.end && !e.alt.empty()) {
    size_t shift = 0;
    while (e.begin && ref[e.begin - 1] == e.alt[e.alt.size() - 1 - shift]) {
      shift = (shift + 1) % e.alt.size();
      --e.begin;
      --e.end;
    }
    if (shift)
      std::rotate(e.alt.begin(), e.alt.end() - shift, e.alt.end());
  } else if (e.alt.empty()) {
    while (e.begin && e.end > e.begin && ref[e.begin - 1] == ref[e.end - 1]) {
      --e.begin;
      --e.end;
    }
  }
  return e;
}
} // namespace
bool Edit::operator<(const Edit &o) const {
  return std::tie(begin, end, alt) < std::tie(o.begin, o.end, o.alt);
}
bool crop_alignment(Alignment &a, size_t begin, size_t end) {
  if (begin >= end || a.end <= begin || a.begin >= end)
    return false;
  auto ops = a.packed_ops ? expand_cigar(a.ops) : a.ops;
  size_t r = a.begin, left = 0, right = 0, query_kept = 0;
  std::string kept;
  for (char op : ops) {
    bool inside = r >= begin && r < end;
    if (inside) {
      kept += op;
      query_kept += op != 'D';
    } else if (op != 'D') {
      (r < begin ? left : right)++;
    }
    if (op != 'I')
      ++r;
  }
  if (kept.empty() || query_kept == 0)
    return false;
  a.begin = std::max(a.begin, begin) - begin;
  a.end = std::min(a.end, end) - begin;
  a.clip_left += left;
  a.clip_right += right;
  a.ops = cigar(kept);
  a.packed_ops = true;
  return true;
}
std::vector<Edit> events(const Alignment &a, const std::string &ref) {
  // Preserve raw alignment order while grouping linked changes. Independently
  // left-normalizing two insertions can move both to the same coordinate and
  // falsely turn one supported haplotype into competing alleles.
  std::string expanded;
  if (a.packed_ops)
    expanded = expand_cigar(a.ops);
  const auto &ops = a.packed_ops ? expanded : a.ops;
  struct Linked {
    Edit raw, normalized;
  };
  std::vector<Linked> linked;
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
    auto canonical = normalize(e, ref);
    // Merge only when normalization reaches/crosses an earlier raw event, not
    // merely because two independent variants are nearby. Include the original
    // intervening matches, retaining the read's actual haplotype and order.
    // Recheck the stack: normalization of a compound can reach another event.
    while (!linked.empty() && canonical.begin <= linked.back().raw.end) {
      auto previous = std::move(linked.back().raw);
      linked.pop_back();
      previous.alt.append(ref, previous.end, e.begin - previous.end);
      previous.alt += e.alt;
      previous.end = e.end;
      e = std::move(previous);
      canonical = normalize(e, ref);
    }
    linked.push_back({std::move(e), std::move(canonical)});
  }
  if (q + a.clip_right != a.bases().size() || r != a.end)
    throw std::runtime_error("Invalid alignment trace");
  std::vector<Edit> out;
  out.reserve(linked.size());
  for (auto &group : linked)
    if (group.normalized.begin != group.normalized.end || !group.normalized.alt.empty())
      out.push_back(std::move(group.normalized));
  return out;
}
std::vector<Accepted> polish(const std::string &ref, const std::vector<Alignment> &alns, size_t nt,
                             uint32_t support, double fraction, size_t max_indel,
                             size_t retain_unique_tiers, bool master_events) {
  std::vector<const Alignment *> ptrs;
  for (auto &a : alns)
    ptrs.push_back(&a);
  return polish(ref, ptrs, nt, support, fraction, max_indel, retain_unique_tiers, master_events);
}
std::vector<Accepted> polish(const std::string &ref, const std::vector<const Alignment *> &alns,
                             size_t nt, uint32_t support, double fraction, size_t max_indel,
                             size_t retain_unique_tiers, bool master_events) {
  if (retain_unique_tiers > nt)
    throw std::runtime_error("Unique tier count exceeds polishing tier count");
  // Filter before inserting into tree/set storage. Homopolymer noise must not
  // inflate the chromosome-wide candidate collection only to be discarded later.
  std::map<Edit, std::set<size_t>> candidates;
  auto eligible = [&](const Edit &e) {
    if (std::max(e.end - e.begin, e.alt.size()) > max_indel || e.begin == 0 ||
        e.end >= ref.size() || e.alt.find('N') != std::string::npos)
      return false;
    return encode(ref.substr(e.begin - 1, e.end - e.begin + 2)).bases !=
           encode(ref.substr(e.begin - 1, 1) + e.alt + ref.substr(e.end, 1)).bases;
  };
  for (size_t i = 0; i < alns.size(); ++i)
    for (const auto &e :
         master_events ? experimental::master_events(*alns[i], ref) : events(*alns[i], ref))
      if (eligible(e))
        candidates[e].insert(i);
  // Interval endpoint indexes give the exact same spanning depth without scanning all reads per
  // event.
  struct Coverage {
    std::vector<size_t> starts, ends;
  };
  std::vector<Coverage> coverage(2 * nt);
  std::vector<std::pair<size_t, const Alignment *>> by_end;
  for (auto a : alns)
    by_end.emplace_back(a->end, a);
  std::sort(by_end.begin(), by_end.end(), [](auto &a, auto &b) { return a.first < b.first; });
  for (auto a : alns) {
    auto &c = coverage[2 * a->tier + a->reverse];
    c.starts.push_back(a->begin);
    c.ends.push_back(a->end);
  }
  for (auto &c : coverage) {
    std::sort(c.starts.begin(), c.starts.end());
    std::sort(c.ends.begin(), c.ends.end());
  }
  std::vector<Accepted> accepted;
  for (const auto &[e, voters] : candidates) {
    // Eligibility was established before insertion; do not repeat allocating
    // flanked RLE comparisons for the same immutable candidate/reference.
    std::vector<uint32_t> depth(nt), votes(nt), forward(nt), reverse(nt), vf(nt), vr(nt);
    for (size_t t = 0; t < nt; ++t)
      for (size_t strand = 0; strand < 2; ++strand) {
        const auto &c = coverage[2 * t + strand];
        // Reads starting strictly before begin minus reads ending at/before begin
        // gives left-flank coverage. The endpoint-range pass below removes those
        // that do not also reach beyond the event's right boundary.
        auto started =
            std::lower_bound(c.starts.begin(), c.starts.end(), e.begin) - c.starts.begin();
        auto ended = std::upper_bound(c.ends.begin(), c.ends.end(), e.begin) - c.ends.begin();
        (strand ? reverse[t] : forward[t]) = started - ended;
      }
    // Remove reads ending inside this event. Only these short endpoint ranges need inspection.
    auto end_it = std::upper_bound(by_end.begin(), by_end.end(), e.begin,
                                   [](size_t p, const auto &a) { return p < a.first; });
    for (; end_it != by_end.end() && end_it->first <= e.end; ++end_it) {
      const auto &a = *end_it->second;
      if (a.begin < e.begin)
        --(a.reverse ? reverse[a.tier] : forward[a.tier]);
    }
    for (size_t t = 0; t < nt; ++t)
      depth[t] = forward[t] + reverse[t];
    for (size_t i : voters) {
      const auto &a = *alns[i];
      // Require aligned bases on both sides; clipped/terminal events never vote.
      if (a.begin >= e.begin || a.end <= e.end)
        continue;
      ++votes[a.tier];
      ++(a.reverse ? vr[a.tier] : vf[a.tier]);
    }
    for (size_t t = 0; t < nt; ++t)
      if (depth[t] >= support) {
        // Select with this tier's depth first; never qualify by accumulating
        // weak tiers. Retained higher unique evidence can support OR veto.
        if (t < retain_unique_tiers)
          for (size_t higher = 0; higher < t; ++higher) {
            depth[t] += depth[higher];
            votes[t] += votes[higher];
            forward[t] += forward[higher];
            reverse[t] += reverse[higher];
            vf[t] += vf[higher];
            vr[t] += vr[higher];
          }
        // Once a tier supplies sufficient depth, disagreement there vetoes the
        // candidate; do not shop lower tiers for a more favorable allele fraction.
        if (votes[t] >= support && double(votes[t]) / depth[t] >= fraction &&
            (!forward[t] || !reverse[t] || (vf[t] && vr[t])))
          accepted.push_back({e, static_cast<int>(t), depth[t], votes[t]});
        break;
      }
  }
  // Conflicting or adjacent edits are left unchanged; do not choose by input order.
  std::vector<Accepted> result;
  size_t prior_end = 0;
  // Candidates are already ordered by Edit. A prefix maximum catches preceding
  // containing intervals; the next start catches any following conflict.
  for (size_t i = 0; i < accepted.size(); ++i) {
    bool conflict = (i > 0 && prior_end >= accepted[i].edit.begin) ||
                    (i + 1 < accepted.size() && accepted[i + 1].edit.begin <= accepted[i].edit.end);
    prior_end = std::max(prior_end, accepted[i].edit.end);
    if (!conflict)
      result.push_back(accepted[i]);
  }
  return result;
}
std::string apply(const std::string &ref, const std::vector<Accepted> &edits) {
  // All edits refer to the original draft, so apply in one forward pass without
  // mutating reference coordinates after each individual change.
  std::string out;
  // Reserve the final size once and append reference slices directly. substr()
  // would allocate chromosome-sized temporary strings between sparse edits.
  size_t length = ref.size();
  size_t previous_end = 0;
  for (const auto &a : edits) {
    const auto &e = a.edit;
    if (e.begin < previous_end || e.begin > e.end || e.end > ref.size())
      throw std::runtime_error("Edits must be ordered, non-overlapping reference intervals");
    length -= e.end - e.begin;
    if (e.alt.size() > out.max_size() - length)
      throw std::length_error("Edited reference exceeds string limit");
    length += e.alt.size();
    previous_end = e.end;
  }
  out.reserve(length);
  size_t p = 0;
  for (const auto &a : edits) {
    out.append(ref, p, a.edit.begin - p);
    out += a.edit.alt;
    p = a.edit.end;
  }
  out.append(ref, p, ref.size() - p);
  return out;
}
size_t lift(size_t p, const std::vector<Accepted> &edits) {
  int64_t shift = 0;
  for (const auto &a : edits) {
    const auto &e = a.edit;
    if (p < e.begin)
      break;
    if (p < e.end)
      return static_cast<size_t>(int64_t(e.begin) + shift);
    shift += int64_t(e.alt.size()) - int64_t(e.end - e.begin);
  }
  return static_cast<size_t>(int64_t(p) + shift);
}
std::string md(const Alignment &a, const std::string &ref, size_t &nm) {
  // MD omits insertions and includes deleted reference bases after '^'. NM
  // counts each inserted/deleted/mismatched base, not the number of gap events.
  std::string expanded;
  if (a.packed_ops)
    expanded = expand_cigar(a.ops);
  const auto &ops = a.packed_ops ? expanded : a.ops;
  size_t r = a.begin, q = a.clip_left, matches = 0;
  nm = 0;
  std::string out;
  for (size_t i = 0; i < ops.size(); ++i) {
    char op = ops[i];
    if (op == '=') {
      if (a.bases().at(q++) != ref.at(r++))
        throw std::runtime_error("Invalid exact match");
      ++matches;
    } else if (op == 'I') {
      ++q;
      ++nm;
    } else if (op == 'X') {
      out += std::to_string(matches);
      matches = 0;
      out += ref.at(r++);
      ++q;
      ++nm;
    } else if (op == 'D') {
      out += std::to_string(matches);
      matches = 0;
      out += '^';
      do {
        out += ref.at(r++);
        ++nm;
        ++i;
      } while (i < ops.size() && ops[i] == 'D');
      --i;
    } else
      throw std::runtime_error("Invalid CIGAR operation");
  }
  if (r != a.end || q + a.clip_right != a.bases().size())
    throw std::runtime_error("CIGAR consumption mismatch");
  return out + std::to_string(matches);
}
} // namespace hpc
