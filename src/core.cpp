// SPDX-License-Identifier: MIT-0
// Pure sequence/consensus operations. Package I/O and pipeline ownership live in main.cpp.
#include "core.hpp"
#include "align/edlib.H"
#include <algorithm>
#include <cmath>
#include <fnmatch.h>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace hpc {
RunModel RunModel::read(std::istream &input) {
  RunModel model;
  std::string line;
  bool header = false;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream row(line);
    if (!header) {
      std::string magic, extra;
      if (!(row >> magic >> model.observed_max >> model.minimum_depth) ||
          magic != "HPC_RUN_MODEL_V1" || model.observed_max < 1 || model.observed_max > 1000 ||
          model.minimum_depth < 1 || model.minimum_depth > 1000000 || (row >> extra))
        throw std::runtime_error("Invalid run model header");
      header = true;
      continue;
    }
    char base;
    size_t truth, loci;
    if (!(row >> base >> truth >> loci) || std::string("ACGT").find(base) == std::string::npos ||
        truth < 1 || truth > 1000 || loci < 1)
      throw std::runtime_error("Invalid run model row");
    auto &rows = model.log_probabilities[base];
    if (truth != rows.size() + 1)
      throw std::runtime_error("Run model lengths must be contiguous from one per base");
    std::vector<double> probabilities(model.observed_max + 1);
    double sum = 0;
    for (double &p : probabilities) {
      if (!(row >> p) || !std::isfinite(p) || p <= 0 || p > 1)
        throw std::runtime_error("Run model probabilities must be finite and positive");
      sum += p;
      p = std::log(p);
    }
    std::string extra;
    if ((row >> extra) || std::abs(sum - 1) > 1e-6)
      throw std::runtime_error("Invalid run model probability row size or normalization");
    rows.push_back(std::move(probabilities));
  }
  if (!input.eof() || !header || model.log_probabilities.empty())
    throw std::runtime_error("Empty or unreadable run model");
  for (const auto &[base, rows] : model.log_probabilities) {
    auto &calls = model.single_allele_calls[base];
    calls.resize(model.observed_max + 1);
    for (size_t observed = 0; observed <= model.observed_max; ++observed) {
      double best = -std::numeric_limits<double>::infinity();
      for (size_t length = 0; length < rows.size(); ++length)
        if (rows[length][observed] > best) {
          best = rows[length][observed];
          calls[observed] = length + 1;
        }
    }
  }
  return model;
}

uint32_t RunModel::estimate_single(char base, uint32_t observed, uint32_t depth) const {
  // Avoid constructing an allocating map for the common single-allele run.
  auto found = log_probabilities.find(base);
  if (depth < minimum_depth || found == log_probabilities.end() || observed < 1 ||
      observed > found->second.size())
    return observed;
  return single_allele_calls.at(base)[std::min<size_t>(observed, observed_max)];
}

uint32_t RunModel::estimate(char base, const Histogram &observations, uint32_t fallback) const {
  auto found = log_probabilities.find(base);
  if (found == log_probabilities.end() || fallback < 1 || fallback > found->second.size())
    return fallback;
  uint64_t depth = 0;
  for (auto [length, count] : observations)
    depth += count;
  if (depth < minimum_depth)
    return fallback;
  if (observations.size() == 1)
    return single_allele_calls.at(
        base)[std::min<size_t>(observations.begin()->first, observed_max)];
  double best = -std::numeric_limits<double>::infinity();
  uint32_t result = fallback;
  for (size_t length = 0; length < found->second.size(); ++length) {
    double score = 0;
    for (auto [observed, count] : observations)
      score += count * found->second[length][std::min<size_t>(observed, observed_max)];
    // Rows are ascending: exact score ties deterministically favor shorter runs.
    if (score > best) {
      best = score;
      result = length + 1;
    }
  }
  return result;
}

Rle encode(const std::string &s) {
  // Keep offsets temporary: eight bytes per read run becomes prohibitive if all
  // reads of a chromosome retain this representation at once.
  Rle r;
  for (size_t i = 0; i < s.size(); ++i) {
    char b = s[i];
    if (b != 'A' && b != 'C' && b != 'G' && b != 'T' && b != 'N')
      throw std::runtime_error("Non-ACGTN sequence");
    if (r.bases.empty() || r.bases.back() != b) {
      r.bases += b;
      r.lengths.push_back(1);
      r.offsets.push_back(i);
    } else {
      if (r.lengths.back() == UINT32_MAX)
        throw std::runtime_error("Run length overflow");
      ++r.lengths.back();
    }
  }
  r.offsets.push_back(s.size());
  return r;
}
std::string reverse_complement(std::string s) {
  std::reverse(s.begin(), s.end());
  for (char &b : s) {
    switch (b) {
    case 'A':
      b = 'T';
      break;
    case 'T':
      b = 'A';
      break;
    case 'C':
      b = 'G';
      break;
    case 'G':
      b = 'C';
      break;
    case 'N':
      break;
    default:
      throw std::runtime_error("Invalid base");
    }
  }
  return s;
}
std::string align(const std::string &q, const std::string &r, double rate) {
  // Most HiFi HPC reads match their supplied interval exactly; avoid edlib then.
  if (q == r)
    return std::string(q.size(), '=');
  if (std::max(q.size(), r.size()) > INT32_MAX)
    throw std::runtime_error("Alignment exceeds edlib limit");
  int k = static_cast<int>(std::ceil(rate * std::max(q.size(), r.size())));
  auto result = merylutil::edlibAlign(
      q.data(), q.size(), r.data(), r.size(),
      merylutil::edlibNewAlignConfig(k, merylutil::EDLIB_MODE_NW, merylutil::EDLIB_TASK_PATH));
  if (result.editDistance < 0) {
    merylutil::edlibFreeAlignResult(result);
    throw std::runtime_error("Alignment exceeds maximum error rate");
  }
  std::string ops;
  ops.reserve(result.alignmentLength);
  for (int i = 0; i < result.alignmentLength; ++i)
    ops += "=IDX"[result.alignment[i]];
  merylutil::edlibFreeAlignResult(result);
  return ops;
}
std::string cigar(const std::string &ops) {
  std::string out;
  for (size_t i = 0; i < ops.size();) {
    size_t j = i + 1;
    while (j < ops.size() && ops[j] == ops[i])
      ++j;
    out += std::to_string(j - i);
    out += ops[i];
    i = j;
  }
  return out.empty() ? "*" : out;
}
std::string expand_cigar(const std::string &packed) {
  std::string out;
  size_t n = 0;
  for (char c : packed) {
    if (c >= '0' && c <= '9')
      n = n * 10 + size_t(c - '0');
    else {
      if (n == 0 || std::string("=XID").find(c) == std::string::npos)
        throw std::runtime_error("Invalid packed trace");
      out.append(n, c);
      n = 0;
    }
  }
  if (n)
    throw std::runtime_error("Truncated packed trace");
  return out;
}
BoundaryAlignment align_boundary(const std::string &query, const std::string &reference, bool left,
                                 bool right, double rate) {
  if (query.empty() || reference.empty() || (!left && !right) || query.size() > INT32_MAX ||
      reference.size() > INT32_MAX)
    throw std::runtime_error("Invalid boundary alignment");
  // Swap roles: consume the template overlap, allowing excess read sequence free.
  auto q = query, r = reference;
  bool reverse = left && !right;
  // Prefix mode frees a right overhang. Reverse both strings (not complement)
  // to apply that same operation to a left overhang, then restore coordinates.
  if (reverse) {
    std::reverse(q.begin(), q.end());
    std::reverse(r.begin(), r.end());
  }
  auto mode = left && right ? merylutil::EDLIB_MODE_HW : merylutil::EDLIB_MODE_SHW;
  auto a = merylutil::edlibAlign(r.data(), r.size(), q.data(), q.size(),
                                 merylutil::edlibNewAlignConfig(int(std::ceil(rate * r.size())),
                                                                mode, merylutil::EDLIB_TASK_PATH));
  if (a.editDistance < 0 || a.numLocations != 1 || !a.startLocations || !a.alignment) {
    merylutil::edlibFreeAlignResult(a);
    throw std::runtime_error("Unresolved boundary alignment");
  }
  BoundaryAlignment out;
  out.query_begin = a.startLocations[0];
  out.query_end = a.endLocations[0] + 1;
  // Swapping template/read roles also swaps insertion/deletion meanings.
  for (int i = 0; i < a.alignmentLength; ++i)
    out.ops += "=DIX"[a.alignment[i]];
  merylutil::edlibFreeAlignResult(a);
  if (reverse) {
    std::reverse(out.ops.begin(), out.ops.end());
    auto b = out.query_begin;
    out.query_begin = query.size() - out.query_end;
    out.query_end = query.size() - b;
  }
  return out;
}
int64_t locate(const std::string &q, const std::string &r, double rate) {
  // LOC finds endpoints without retaining a traceback for the whole template.
  // Reject tied placements: a convenient repeat hit is not offset evidence.
  if (q.empty() || q.size() > INT32_MAX || r.size() > INT32_MAX)
    return -1;
  auto a = merylutil::edlibAlign(q.data(), q.size(), r.data(), r.size(),
                                 merylutil::edlibNewAlignConfig(int(std::ceil(rate * q.size())),
                                                                merylutil::EDLIB_MODE_HW,
                                                                merylutil::EDLIB_TASK_LOC));
  int64_t pos = -1;
  if (a.editDistance >= 0 && a.numLocations == 1 && a.startLocations && a.startLocations[0] >= 0)
    pos = a.startLocations[0];
  merylutil::edlibFreeAlignResult(a);
  return pos;
}
int64_t offset_cluster(std::vector<int64_t> x, size_t required, int64_t tolerance) {
  std::sort(x.begin(), x.end());
  size_t best = 0, begin = 0;
  // A sliding window finds the largest bounded-width cluster in sorted order.
  // Strict improvement preserves the first (lowest-coordinate) cluster on ties.
  for (size_t i = 0, j = 0; i < x.size(); ++i) {
    while (j < x.size() && x[j] - x[i] <= tolerance)
      ++j;
    if (j - i > best) {
      best = j - i;
      begin = i;
    }
  }
  if (best < required || best * 4 < x.size() * 3)
    throw std::runtime_error("Insufficient consistent offset evidence");
  return x[begin + (best - 1) / 2];
}
int tier_of(const std::string &name, const std::vector<Tier> &tiers) {
  // First-match semantics implement explicit user priority; names imply no type.
  for (size_t i = 0; i < tiers.size(); ++i)
    if (fnmatch(tiers[i].pattern.c_str(), name.c_str(), 0) == 0)
      return i;
  return -1;
}
Call call(const std::vector<Histogram> &tiers, uint32_t min_support, bool rounded_mean) {
  // Qualify tiers by total observations, not by the modal count. A lower tier
  // cannot outvote a qualifying higher one, even with much greater coverage.
  for (size_t t = 0; t < tiers.size(); ++t) {
    uint64_t n = 0;
    uint32_t peak = 0;
    for (auto [k, c] : tiers[t]) {
      n += c;
      peak = std::max(peak, c);
    }
    if (n < min_support)
      continue;
    if (n > UINT32_MAX)
      throw std::runtime_error("Coverage overflow");
    if (rounded_mean) {
      long double total = 0;
      for (auto [length, count] : tiers[t])
        total += static_cast<long double>(length) * count;
      return {static_cast<uint32_t>(total / n + 0.5L), static_cast<uint32_t>(n), peak,
              static_cast<int>(t)};
    }
    uint32_t best = 0;
    uint64_t distance = UINT64_MAX;
    // Tied modes minimize total absolute deviation; sorted keys and strict '<'
    // make the final tie deterministic in favor of the shorter run.
    for (auto [k, c] : tiers[t])
      if (c == peak) {
        uint64_t d = 0;
        for (auto [v, count] : tiers[t])
          d += uint64_t(k > v ? k - v : v - k) * count;
        if (d < distance) {
          best = k;
          distance = d;
        }
      }
    return {best, static_cast<uint32_t>(n), peak, static_cast<int>(t)};
  }
  return {};
}
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
  // Coalesce consecutive non-matches into a single replacement allele. This
  // allows substitutions, indels, and combined events to use one voting policy.
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
    // Normalize pure indels by rotating their repeat sequence to the left.
    if (e.begin == e.end) {
      while (e.begin > 0 && !e.alt.empty() && ref[e.begin - 1] == e.alt.back()) {
        e.alt = e.alt.back() + e.alt.substr(0, e.alt.size() - 1);
        --e.begin;
        --e.end;
      }
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
std::vector<Accepted> polish(const std::string &ref, const std::vector<Alignment> &alns, size_t nt,
                             uint32_t support, double fraction, size_t max_indel) {
  std::vector<const Alignment *> ptrs;
  for (auto &a : alns)
    ptrs.push_back(&a);
  return polish(ref, ptrs, nt, support, fraction, max_indel);
}
std::vector<Accepted> polish(const std::string &ref, const std::vector<const Alignment *> &alns,
                             size_t nt, uint32_t support, double fraction, size_t max_indel) {
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
    for (const auto &e : events(*alns[i], ref))
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
    if (std::max(e.end - e.begin, e.alt.size()) > max_indel || e.begin == 0 ||
        e.end >= ref.size() || e.alt.find('N') != std::string::npos)
      continue;
    // A flanked compression comparison excludes pure run-length changes.
    std::string before = ref.substr(e.begin - 1, e.end - e.begin + 2);
    std::string after = ref.substr(e.begin - 1, 1) + e.alt + ref.substr(e.end, 1);
    if (encode(before).bases == encode(after).bases)
      continue;
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
  size_t p = 0;
  for (const auto &a : edits) {
    out += ref.substr(p, a.edit.begin - p);
    out += a.edit.alt;
    p = a.edit.end;
  }
  out += ref.substr(p);
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
