// SPDX-License-Identifier: MIT-0
#include "repeat_refinement.hpp"
#include "align/edlib.H"
#include <algorithm>
#include <array>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
namespace hpc {
namespace {
int code(char b) {
  return b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : 4;
}
size_t unique_at(const std::string &s, const std::string &word) {
  size_t p = s.find(word);
  return p != std::string::npos && s.find(word, p + 1) == std::string::npos ? p : std::string::npos;
}
struct Molecule {
  uint32_t id;
  int tier;
  bool unique;
  std::string sequence;
};
int distance(const std::string &q, const std::string &r) {
  auto a = merylutil::edlibAlign(
      q.data(), q.size(), r.data(), r.size(),
      merylutil::edlibNewAlignConfig(-1, merylutil::EDLIB_MODE_NW, merylutil::EDLIB_TASK_DISTANCE));
  int d = a.editDistance;
  merylutil::edlibFreeAlignResult(a);
  if (d < 0)
    throw std::runtime_error("Repeat distance failed");
  return d;
}
bool informative_anchor(const std::string &s) {
  std::array<size_t, 5> counts{};
  for (char b : s)
    ++counts[code(b)];
  if (counts[4])
    return false;
  std::sort(counts.begin(), counts.end(), std::greater<size_t>());
  return counts[0] <= 30 && counts[0] + counts[1] <= 45;
}
// Raw majority proposal. Unlike independent HP votes, the complete raw trace
// retains the run-length pattern that distinguishes neighboring repeat units.
std::string propose(const std::string &seed, const std::vector<const Molecule *> &reads) {
  std::vector<std::array<unsigned, 5>> bases(seed.size());
  std::map<size_t, std::map<std::string, unsigned>> insertions;
  for (auto read : reads) {
    auto ops = align(read->sequence, seed, 1.0);
    size_t q = 0, r = 0;
    std::string insertion;
    auto flush = [&]() {
      if (!insertion.empty()) {
        ++insertions[r][insertion];
        insertion.clear();
      }
    };
    for (char op : ops) {
      if (op == 'I') {
        insertion += read->sequence[q++];
        continue;
      }
      flush();
      if (op == 'D') {
        ++bases.at(r)[4];
        ++r;
      } else {
        ++bases.at(r)[code(read->sequence[q])];
        ++q;
        ++r;
      }
    }
    flush();
  }
  std::string result;
  for (size_t i = 0; i <= seed.size(); ++i) {
    auto it = insertions.find(i);
    if (it != insertions.end())
      for (const auto &[s, n] : it->second)
        if (2 * n > reads.size()) {
          result += s;
          break;
        }
    if (i == seed.size())
      break;
    auto best = std::max_element(bases[i].begin(), bases[i].end());
    size_t b = best - bases[i].begin();
    if (2 * (*best) <= reads.size())
      result += seed[i];
    else if (b < 4)
      result += "ACGT"[b];
  }
  return result;
}
} // namespace
std::vector<RepeatInterval> detect_repeats(const std::string &s, const std::vector<Call> &calls) {
  if (s.size() != calls.size())
    throw std::invalid_argument("Repeat call/backbone size mismatch");
  constexpr size_t width = 32, join = 32;
  std::vector<RepeatInterval> regions;
  std::array<size_t, 9> periodic_matches{};
  std::array<size_t, 5> counts{};
  for (size_t i = 0; i < s.size(); ++i) {
    ++counts[code(s[i])];
    for (size_t period = 2; period <= 8; ++period) {
      if (i >= period)
        periodic_matches[period] += s[i] == s[i - period];
      if (i >= width) {
        size_t j = i - (width - period);
        periodic_matches[period] -= s[j] == s[j - period];
      }
    }
    if (i >= width)
      --counts[code(s[i - width])];
    if (i + 1 < width || (i + 1 - width) % 8)
      continue;
    auto sorted = counts;
    std::sort(sorted.begin(), sorted.begin() + 4, std::greater<size_t>());
    bool periodic = false;
    for (size_t period = 2; period <= 8; ++period)
      periodic |= periodic_matches[period] * 10 >= (width - period) * 9;
    if (counts[4] || (sorted[0] + sorted[1] < 30 && !periodic))
      continue;
    size_t b = i + 1 - width, e = i + 1;
    if (!regions.empty() && b <= regions.back().end + join)
      regions.back().end = e;
    else
      regions.push_back({b, e});
  }
  regions.erase(std::remove_if(regions.begin(), regions.end(),
                               [&](const RepeatInterval &r) {
                                 size_t uncertain = 0;
                                 for (size_t i = r.begin; i < r.end; ++i)
                                   uncertain +=
                                       calls[i].depth >= 3 &&
                                       uint64_t(calls[i].peak) * 10 < uint64_t(calls[i].depth) * 9;
                                 return uncertain < 2;
                               }),
                regions.end());
  return regions;
}
std::string project_repeat_backbone(const std::string &q, const std::string &b, size_t max_cells) {
  if (b.empty())
    throw std::invalid_argument("Empty repeat backbone");
  size_t stride = b.size() + 1;
  if (stride > max_cells || q.size() + 1 > max_cells / stride)
    throw std::length_error("Repeat projection budget");
  std::vector<uint8_t> trace((q.size() + 1) * stride);
  std::vector<int> previous(stride), current(stride);
  for (size_t k = 1; k < stride; ++k) {
    previous[k] = int(k);
    trace[k] = 'D';
  }
  for (size_t i = 1; i <= q.size(); ++i) {
    current[0] = int(i);
    trace[i * stride] = 'I';
    for (size_t k = 1; k < stride; ++k) {
      int mismatch = q[i - 1] != b[k - 1], score = previous[k - 1] + mismatch;
      char op = 'A';
      if (previous[k] + mismatch < score) {
        score = previous[k] + mismatch;
        op = 'R';
      }
      if (previous[k] + 1 < score) {
        score = previous[k] + 1;
        op = 'I';
      }
      if (current[k - 1] + 1 < score) {
        score = current[k - 1] + 1;
        op = 'D';
      }
      current[k] = score;
      trace[i * stride + k] = op;
    }
    previous.swap(current);
  }
  size_t i = q.size(), k = b.size();
  std::string out;
  while (i || k) {
    char op = trace[i * stride + k];
    if (op == 'A') {
      out += b[k - 1];
      --i;
      --k;
    } else if (op == 'R') {
      out += b[k - 1];
      --i;
    } else if (op == 'D') {
      out += b[k - 1];
      --k;
    } else if (op == 'I')
      --i;
    else
      throw std::logic_error("Invalid repeat projection trace");
  }
  std::reverse(out.begin(), out.end());
  return out;
}
std::vector<RepeatInterval>
refine_repeats(const std::string &backbone, std::vector<Call> &calls, std::vector<RepeatRead> reads,
               size_t tiers, const std::map<int, RunModel> &models, const RepeatConfig &cfg,
               const std::function<void(const RepeatDecision &)> &report) {
  auto regions = detect_repeats(backbone, calls);
  std::vector<RepeatInterval> changed;
  std::sort(reads.begin(), reads.end(), [](auto &a, auto &b) { return a.begin < b.begin; });
  // Avoid a chromosome-wide raw draft/prefix allocation here. Detection and
  // context use run lengths; only an admitted window is expanded.
  size_t last_end = 0, cursor = 0;
  std::multimap<size_t, const RepeatRead *> active;
  for (auto core : regions) {
    RepeatDecision info;
    info.core = core;
    info.window = core;
    auto done = [&](const char *status) {
      info.status = status;
      report(info);
    };
    size_t b = core.begin, e = core.end, left = 0, right = 0;
    while (b && left < cfg.flank_bases) {
      --b;
      left += calls[b].length;
    }
    while (e < backbone.size() && right < cfg.flank_bases) {
      right += calls[e].length;
      ++e;
    }
    info.window = {b, e};
    if (left < cfg.flank_bases || right < cfg.flank_bases) {
      done("TERMINAL");
      continue;
    }
    if (b < last_end) {
      done("OVERLAPPING_WINDOW");
      continue;
    }
    size_t length = 0;
    for (size_t j = b; j < e; ++j) {
      length += calls[j].length;
      if (length > cfg.max_window_bases)
        break;
    }
    if (length > cfg.max_window_bases) {
      done("WINDOW_BUDGET");
      continue;
    }
    std::vector<const RepeatRead *> spanning;
    // Read intervals are sorted. Stop at the left boundary; do not inspect
    // sequence storage for nonoverlapping placements.
    while (cursor < reads.size() && reads[cursor].begin <= b) {
      active.emplace(reads[cursor].end, &reads[cursor]);
      ++cursor;
    }
    while (!active.empty() && active.begin()->first < core.end)
      active.erase(active.begin());
    for (auto it = active.lower_bound(e); it != active.end(); ++it)
      if (it->second->begin <= b)
        spanning.push_back(it->second);
    if (spanning.size() > cfg.max_reads) {
      done("READ_BUDGET");
      continue;
    }
    if (spanning.size() < cfg.min_support) {
      done("LOW_SPANNING_DEPTH");
      continue;
    }
    std::string draft;
    std::vector<size_t> offsets{0};
    for (size_t j = b; j < e; ++j) {
      draft.append(calls[j].length, backbone[j]);
      offsets.push_back(draft.size());
    }
    // Restrict anchors to run boundaries. Then the refined interior maps back
    // to exactly the original run indexes; no partial run can be overwritten.
    auto anchor = [&](bool first) -> std::pair<size_t, std::string> {
      std::vector<size_t> best_counts(2 * tiers, 0);
      size_t best = std::string::npos, next_position = 0;
      for (size_t p : offsets) {
        if (first ? (p > left / 2 || p + 51 >= draft.size())
                  : (p < draft.size() - right / 2 || p < 51 || p == draft.size()))
          continue;
        if (p < next_position)
          continue;
        next_position = p + 64;
        size_t start = first ? p : p - 51;
        auto word = draft.substr(start, 51);
        if (!informative_anchor(word) || unique_at(draft, word) == std::string::npos)
          continue;
        std::vector<size_t> support(2 * tiers, 0);
        for (auto r : spanning)
          if (unique_at(r->alignment->bases(), word) != std::string::npos)
            ++support.at(size_t(r->alignment->tier) + (r->unique ? 0 : tiers));
        if (best == std::string::npos || support > best_counts) {
          best = p;
          best_counts = std::move(support);
        }
        // Fixed stride in raw space bounds anchor scans, independently of run lengths.
        if (first && p >= 256)
          break;
      }
      return best == std::string::npos ? std::make_pair(best, std::string())
                                       : std::make_pair(best, first ? draft.substr(best, 51)
                                                                    : draft.substr(best - 51, 51));
    };
    auto [lb, la] = anchor(true);
    auto [re, ra] = anchor(false);
    if (lb == std::string::npos || re == std::string::npos || lb + 51 >= re - 51) {
      done("NO_ANCHORS");
      continue;
    }
    // Include both anchors in the pattern. They are trusted context and are
    // preserved by the same positive-run backbone invariant as the interior.
    auto ib = std::lower_bound(offsets.begin(), offsets.end(), lb) - offsets.begin();
    auto ie = std::lower_bound(offsets.begin(), offsets.end(), re) - offsets.begin();
    std::string initial = draft.substr(lb, re - lb), hpc = backbone.substr(b + ib, ie - ib);
    std::vector<Molecule> molecules;
    std::set<uint32_t> used;
    for (auto r : spanning) {
      auto &a = *r->alignment;
      auto &s = a.bases();
      size_t x = unique_at(s, la), y = unique_at(s, ra);
      if (x == std::string::npos || y == std::string::npos || x + 51 > y ||
          y + 51 - x > cfg.max_window_bases)
        continue;
      if (!used.insert(a.id).second)
        continue;
      molecules.push_back({a.id, a.tier, r->unique, s.substr(x, y + 51 - x)});
    }
    info.reads = molecules.size();
    std::vector<const Molecule *> scoring;
    for (size_t phase = 0; phase < (cfg.ambiguous_fallback ? 2u : 1u) && scoring.empty(); ++phase)
      for (size_t t = 0; t < tiers; ++t) {
        std::vector<const Molecule *> group;
        for (auto &m : molecules)
          if (size_t(m.tier) == t && (m.unique || phase))
            group.push_back(&m);
        if (group.size() >= cfg.min_support) {
          scoring = std::move(group);
          info.tier = int(t);
          info.pooled = phase;
          break;
        }
      }
    if (scoring.empty()) {
      done("LOW_ANCHORED_DEPTH");
      continue;
    }
    std::vector<std::string> candidates{initial};
    auto add = [&](const std::string &s) {
      if (s.size() <= cfg.max_window_bases && encode(s).bases == hpc &&
          std::find(candidates.begin(), candidates.end(), s) == candidates.end() &&
          candidates.size() < cfg.max_candidates)
        candidates.push_back(s);
    };
    for (size_t phase = 0; phase < (cfg.ambiguous_fallback ? 2u : 1u); ++phase)
      for (size_t t = 0; t < tiers; ++t) {
        std::map<std::string, size_t> counts;
        for (auto &m : molecules)
          if (size_t(m.tier) == t && (m.unique || phase))
            ++counts[m.sequence];
        std::vector<std::pair<size_t, std::string>> order;
        for (auto &[s, n] : counts)
          if (encode(s).bases == hpc)
            order.emplace_back(n, s);
        std::sort(order.begin(), order.end(), [](auto &x, auto &y) {
          return x.first != y.first ? x.first > y.first : x.second < y.second;
        });
        for (size_t j = 0; j < std::min<size_t>(2, order.size()); ++j)
          add(order[j].second);
      }
    auto score = [&](const std::string &s) {
      int64_t total = 0;
      for (auto r : scoring)
        total += distance(r->sequence, s);
      return total;
    };
    int64_t baseline = score(initial), best = baseline;
    size_t winner = 0;
    auto choose = [&]() {
      for (size_t j = 1; j < candidates.size(); ++j) {
        auto value = score(candidates[j]);
        if (value < best) {
          best = value;
          winner = j;
        }
      }
    };
    choose();
    bool projection_exhausted = false;
    if (!winner) {
      size_t budget = cfg.window_projection_cells;
      auto project = [&](const std::string &s) {
        size_t stride = hpc.size() + 1;
        if (stride > cfg.projection_cells || s.size() + 1 > cfg.projection_cells / stride ||
            (s.size() + 1) * stride > budget)
          throw std::length_error("Repeat projection budget");
        budget -= (s.size() + 1) * stride;
        return project_repeat_backbone(s, hpc, cfg.projection_cells);
      };
      // Small fixed number of starts/rounds; every candidate is still scored
      // against every admitted scoring-tier molecule, never a read subsample.
      try {
        for (size_t start = 0; start < std::min<size_t>(3, scoring.size() + 1); ++start) {
          std::string seed = start ? project(scoring[start - 1]->sequence) : initial;
          add(seed);
          for (size_t round = 0; round < 2 && candidates.size() < cfg.max_candidates; ++round) {
            auto next = project(propose(seed, scoring));
            add(next);
            if (next == seed)
              break;
            seed = std::move(next);
          }
        }
      } catch (const std::length_error &) {
        projection_exhausted = true;
      }
      choose();
    }
    info.candidates = candidates.size();
    info.before_score = baseline;
    info.after_score = best;
    if (!winner || baseline - best < int64_t(scoring.size())) {
      done(projection_exhausted ? "PROJECTION_BUDGET" : "NO_SCORE_GAIN");
      continue;
    }
    size_t improved = 0;
    for (auto r : scoring)
      improved += distance(r->sequence, candidates[winner]) < distance(r->sequence, initial);
    if (2 * improved <= scoring.size()) {
      done("NO_MOLECULE_MAJORITY");
      continue;
    }
    auto seed = encode(candidates[winner]);
    std::vector<std::map<int, Histogram>> hist(hpc.size());
    for (auto &m : molecules) {
      auto forward = align(m.sequence, candidates[winner], 1.0);
      std::string q(m.sequence.rbegin(), m.sequence.rend()),
          r(candidates[winner].rbegin(), candidates[winner].rend());
      auto backward = align(q, r, 1.0);
      std::reverse(backward.begin(), backward.end());
      std::array<std::vector<size_t>, 2> maps;
      for (size_t direction = 0; direction < 2; ++direction) {
        auto &map = maps[direction];
        map.assign(candidates[winner].size(), SIZE_MAX);
        size_t x = 0, y = 0;
        for (char op : direction ? backward : forward) {
          if (op == '=')
            map[y] = x;
          x += op != 'D';
          y += op != 'I';
        }
      }
      for (size_t j = 1; j + 1 < hpc.size(); ++j) {
        size_t x = seed.offsets[j] - 1, y = seed.offsets[j + 1];
        if (maps[0][x] == SIZE_MAX || maps[0][y] == SIZE_MAX || maps[0][x] != maps[1][x] ||
            maps[0][y] != maps[1][y])
          continue;
        size_t qb = maps[0][x] + 1, qe = maps[0][y];
        if (qe <= qb || !std::all_of(m.sequence.begin() + qb, m.sequence.begin() + qe,
                                     [&](char c) { return c == hpc[j]; }))
          continue;
        ++hist[j][m.tier + (m.unique ? 0 : int(tiers))][uint32_t(qe - qb)];
      }
    }
    for (size_t j = 1; j + 1 < hpc.size(); ++j) {
      // Do not modify flanking context outside the detected repeat itself.
      size_t run = b + ib + j;
      if (run < core.begin || run >= core.end)
        continue;
      bool decided = false;
      for (size_t phase = 0; phase < (cfg.ambiguous_fallback ? 2u : 1u) && !decided; ++phase)
        for (size_t t = 0; t < tiers; ++t) {
          auto h = hist[j][int(t)];
          if (phase)
            for (auto [l, n] : hist[j][int(t + tiers)])
              h[l] += n;
          auto c = call({h}, cfg.min_support);
          if (c.depth < cfg.min_support)
            continue;
          decided = true;
          c.tier = int(t + (phase ? tiers : 0));
          auto model = models.find(int(t));
          auto uncalibrated = c.length;
          if (model != models.end())
            c.length = model->second.estimate(hpc[j], h, c.length);
          if (c.length != calls[run].length) {
            calls[run] = c;
            info.updates.push_back(
                {run, c.length != uncalibrated, bool(phase && !hist[j][int(t)].empty())});
            ++info.changed_runs;
          }
          break;
        }
    }
    if (info.changed_runs) {
      changed.push_back(core);
      last_end = e;
      done("ACCEPTED");
    } else
      done("NO_LENGTH_CHANGE");
  }
  return changed;
}
} // namespace hpc
