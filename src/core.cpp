// SPDX-License-Identifier: MIT-0
// Sequence, alignment and empirical run calling. See models.cpp and polish.cpp for\n// calibration
// and raw-base correction; package ownership lives in main.cpp.
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
std::string align(const std::string &q, const std::string &r, double rate, size_t distance_bound) {
  // Most HiFi HPC reads match their supplied interval exactly; avoid edlib then.
  if (q == r)
    return std::string(q.size(), '=');
  if (std::max(q.size(), r.size()) > INT32_MAX)
    throw std::runtime_error("Alignment exceeds edlib limit");
  int k = static_cast<int>(std::ceil(rate * std::max(q.size(), r.size())));
  k = static_cast<int>(std::min<size_t>(k, distance_bound));
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
std::vector<uint8_t> local_hpc_votes(const std::string &ops, size_t radius, double max_error) {
  if (radius == 0 || !std::isfinite(max_error) || max_error < 0 || max_error > 1)
    throw std::runtime_error("Invalid local HPC filter parameters");
  if (ops.find_first_not_of("=XID") != std::string::npos)
    throw std::runtime_error("Invalid expanded HPC trace");
  if (ops.find_first_not_of('=') == std::string::npos)
    return {};
  // A rolling sum avoids a per-base integer prefix array. The caller releases
  // this temporary mask after voting for one read, never retaining it per tig.
  std::vector<uint8_t> mask(ops.size());
  size_t left = 0, right = 0, errors = 0;
  for (size_t i = 0; i < ops.size(); ++i) {
    size_t begin = i > radius ? i - radius : 0;
    size_t end = i + std::min(radius, ops.size() - i - 1) + 1;
    while (right < end)
      errors += ops[right++] != '=';
    while (left < begin)
      errors -= ops[left++] != '=';
    mask[i] = ops[i] == '=' && double(errors) <= max_error * double(end - begin);
  }
  return mask;
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
  // In a one-sided alignment the inward start is fixed. Multiple optimal
  // endpoints describe clipping uncertainty, not alternative genomic loci.
  // Edlib traces the first (shortest) optimal endpoint; retain that conservative
  // choice. Two-sided infix alignment still needs an unambiguous placement.
  if (a.editDistance < 0 || a.numLocations < 1 || (left && right && a.numLocations != 1) ||
      !a.startLocations || !a.alignment) {
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
} // namespace hpc
