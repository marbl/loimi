// SPDX-License-Identifier: MIT-0
#include "offsets.hpp"
#include "core.hpp"
#include "align/edlib.H"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace hpc {
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
  bool tied = a.editDistance >= 0 && a.numLocations != 1;
  bool resolved = a.editDistance >= 0 && a.numLocations == 1 && a.startLocations;
  int64_t p = resolved ? a.startLocations[0] : -1;
  // A noisy reverse alignment pinned to a read boundary can be a truncated
  // template, rather than a truly enclosing read. Let internal chains test it.
  bool truncated = resolved && enclosing && a.editDistance > 0 &&
                   (p == 0 || a.endLocations[0] + 1 == int(q.size()));
  merylutil::edlibFreeAlignResult(a);
  if (resolved && !truncated)
    return {enclosing ? begin + p : begin - p};
  if (tied)
    return {};

  size_t width = std::min<size_t>(2000, r.size() / 4);
  if (width < 128 || q.size() < 2 * width)
    return {};
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
  return result;
}

OffsetChoice choose_offset(const std::vector<OffsetEvidence> &reads, size_t required,
                           bool unique_only) {
  std::set<uint32_t> informative;
  std::vector<int64_t> starts;
  for (const auto &read : reads) {
    if (unique_only && !read.unique)
      continue;
    if (!read.offsets.empty())
      informative.insert(read.id);
    starts.insert(starts.end(), read.offsets.begin(), read.offsets.end());
  }
  std::sort(starts.begin(), starts.end());
  starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
  auto votes = [&](int64_t start, bool only_unique) {
    std::map<uint32_t, std::vector<int64_t>> by_id;
    for (const auto &read : reads) {
      if (only_unique && !read.unique)
        continue;
      for (auto x : read.offsets)
        if (x >= start && x - start <= 20)
          by_id[read.id].push_back(x);
    }
    std::vector<int64_t> v;
    for (auto &[id, offsets] : by_id) {
      std::sort(offsets.begin(), offsets.end());
      v.push_back(offsets[(offsets.size() - 1) / 2]);
    }
    std::sort(v.begin(), v.end());
    return v;
  };
  std::vector<int64_t> best;
  int64_t start = 0;
  for (auto x : starts) {
    auto v = votes(x, unique_only);
    if (v.size() > best.size()) {
      start = x;
      best = std::move(v);
    }
  }
  size_t runner = 0, unique_runner = 0;
  for (auto x : starts)
    if (std::abs(x - start) > 40) {
      runner = std::max(runner, votes(x, unique_only).size());
      unique_runner = std::max(unique_runner, votes(x, true).size());
    }
  if (best.size() < required || best.size() * 5 < informative.size() * 3 ||
      runner * 2 >= best.size() ||
      (!unique_only && unique_runner >= 3 && unique_runner >= votes(start, true).size())) {
    std::ostringstream error;
    error << "offset candidate=" << start << " support=" << best.size()
          << " informative=" << informative.size() << " runner_up=" << runner
          << " unique_runner=" << unique_runner << " required=" << required;
    throw std::runtime_error(error.str());
  }
  return {best[(best.size() - 1) / 2], best.size(), informative.size(), runner};
}

bool validate_offset_overlap(const std::string &q, const std::string &r, int64_t begin, int64_t end,
                             int64_t offset, double rate) {
  auto b = begin - offset, e = end - offset;
  auto lo = std::max<int64_t>(0, b), hi = std::min<int64_t>(r.size(), e);
  if (hi <= lo || hi - lo < 128 || hi - lo < int64_t(std::min<size_t>(1000, r.size())))
    return false;
  try {
    auto target = r.substr(size_t(lo), size_t(hi - lo));
    if (b < 0 || e > int64_t(r.size()))
      (void)hpc::align_boundary(q, target, b<0, e> int64_t(r.size()), rate);
    else
      (void)hpc::align(q, target, rate);
    return true;
  } catch (const std::runtime_error &) {
    return false;
  }
}
} // namespace hpc
