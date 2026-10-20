// SPDX-License-Identifier: MIT-0
// Opt-in same-tier borrowing policy. No model loading or default estimation lives here.
#include "core.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hpc {
std::vector<double> RunModel::likelihood_scores(char base, const Histogram &h) const {
  auto found = log_probabilities.find(base);
  auto mode = call({h}, 1);
  if (found == log_probabilities.end() || mode.depth < minimum_depth ||
      mode.length > found->second.size())
    return {};
  std::vector<double> scores(found->second.size());
  for (size_t i = 0; i < scores.size(); ++i)
    for (auto [length, count] : h)
      scores[i] += count * found->second[i][std::min<size_t>(length, observed_max)];
  return scores;
}
bool RunModel::uncertain(char base, const Histogram &h) const {
  if (h.empty())
    return false;
  auto scores = likelihood_scores(base, h);
  if (scores.empty()) {
    auto mode = call({h}, 1);
    return double(mode.peak) / mode.depth < .9;
  }
  if (scores.size() < 2)
    return false;
  std::partial_sort(scores.begin(), scores.begin() + 2, scores.end(), std::greater<double>());
  return scores[0] - scores[1] < std::log(100.0);
}
uint32_t RunModel::borrow(char base, const Histogram &u, const Histogram &a,
                          uint32_t current) const {
  auto ua = call({u}, 1), aa = call({a}, 1);
  if (aa.depth < 5 || !uncertain(base, u))
    return current;
  uint32_t alt = estimate(base, a, aa.length);
  auto scores = likelihood_scores(base, u);
  if (!scores.empty()) {
    if (alt < 1 || alt > scores.size() ||
        *std::max_element(scores.begin(), scores.end()) - scores[alt - 1] > std::log(20.0))
      return current;
  } else {
    auto found = u.find(alt);
    if (found == u.end() || double(found->second) / ua.depth < .2)
      return current;
  }
  double weight = std::min(5.0, .5 * ua.depth) / aa.depth;
  std::map<uint32_t, double> mixed;
  for (auto [length, count] : u)
    mixed[length] += count;
  for (auto [length, count] : a)
    mixed[length] += weight * count;
  // Reproduce the ordinary mode's absolute-deviation then shorter tie break.
  uint32_t modal = 0;
  double peak = -1, deviation = 0;
  for (auto [length, count] : mixed) {
    double dev = 0;
    for (auto [other, n] : mixed)
      dev += std::abs(double(length) - other) * n;
    if (count > peak || (count == peak && dev < deviation)) {
      modal = length;
      peak = count;
      deviation = dev;
    }
  }
  auto found = log_probabilities.find(base);
  if (found == log_probabilities.end() || modal > found->second.size() ||
      ua.depth + weight * aa.depth < minimum_depth)
    return modal;
  double best = -std::numeric_limits<double>::infinity();
  uint32_t result = modal;
  for (size_t i = 0; i < found->second.size(); ++i) {
    double score = 0;
    for (auto [length, count] : mixed)
      score += count * found->second[i][std::min<size_t>(length, observed_max)];
    if (score > best) {
      best = score;
      result = i + 1;
    }
  }
  return result;
}
} // namespace hpc
