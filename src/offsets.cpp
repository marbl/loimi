// SPDX-License-Identifier: MIT-0
#include "offsets.hpp"
#include <algorithm>
#include <cstdlib>
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
} // namespace hpc
