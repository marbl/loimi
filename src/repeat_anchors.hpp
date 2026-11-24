// SPDX-License-Identifier: MIT-0
#pragma once
#include <array>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace hpc::repeat_detail {
// Bounded, window-local Aho-Corasick matcher for exact 51-base DNA anchors.
// Search the ENTIRE oriented read, including clips: narrowing the search to
// the expected placement would change whether an anchor is globally unique.
// No chromosome/read-wide index is retained between candidate windows.
class UniqueAnchors {
public:
  static constexpr size_t width = 51;
  explicit UniqueAnchors(const std::vector<std::string> &words) : count_(words.size()) {
    nodes_.emplace_back();
    for (size_t id = 0; id < words.size(); ++id) {
      if (words[id].size() != width)
        throw std::invalid_argument("Repeat anchors must have length 51");
      size_t state = 0;
      for (unsigned char c : words[id]) {
        auto base = alphabet()[c];
        if (base == 4)
          throw std::invalid_argument("Repeat anchors must contain only ACGT");
        int child = nodes_[state].next[base];
        if (child < 0) {
          child = int(nodes_.size());
          nodes_[state].next[base] = child;
          nodes_.emplace_back();
        }
        state = size_t(child);
      }
      // Duplicate queries share a state but retain their separate result slots.
      nodes_[state].outputs.push_back(id);
    }
    std::queue<size_t> pending;
    for (size_t base = 0; base < 4; ++base) {
      int child = nodes_[0].next[base];
      if (child < 0)
        nodes_[0].next[base] = 0;
      else
        pending.push(size_t(child));
    }
    while (!pending.empty()) {
      size_t state = pending.front();
      pending.pop();
      for (size_t base = 0; base < 4; ++base) {
        int child = nodes_[state].next[base];
        int fallback = nodes_[nodes_[state].failure].next[base];
        if (child < 0)
          nodes_[state].next[base] = fallback;
        else {
          nodes_[child].failure = size_t(fallback);
          pending.push(size_t(child));
        }
      }
    }
  }

  // Return the zero-based start only for exactly one occurrence, otherwise
  // npos. Count overlapping matches, just like find(word, first + 1).
  std::vector<size_t> find_unique(std::string_view text) const {
    std::vector<size_t> positions(count_, std::string::npos);
    if (!count_)
      return positions;
    std::vector<uint8_t> seen(count_);
    size_t state = 0;
    const auto &letters = alphabet();
    for (size_t i = 0; i < text.size(); ++i) {
      auto base = letters[static_cast<unsigned char>(text[i])];
      state = base == 4 ? 0 : size_t(nodes_[state].next[base]);
      // All queries have identical length. A proper failure suffix is shorter
      // than 51, so it cannot contribute an inherited terminal output.
      for (size_t id : nodes_[state].outputs) {
        if (!seen[id]) {
          positions[id] = i + 1 - width;
          seen[id] = 1;
        } else {
          positions[id] = std::string::npos;
          seen[id] = 2; // Saturate: a third match must never restore uniqueness.
        }
      }
    }
    return positions;
  }

private:
  struct Node {
    std::array<int, 4> next{{-1, -1, -1, -1}};
    size_t failure = 0;
    std::vector<size_t> outputs;
  };
  static const std::array<uint8_t, 256> &alphabet() {
    static const auto letters = [] {
      std::array<uint8_t, 256> result{};
      result.fill(4);
      result['A'] = 0;
      result['C'] = 1;
      result['G'] = 2;
      result['T'] = 3;
      return result;
    }();
    return letters;
  }
  size_t count_;
  std::vector<Node> nodes_;
};
} // namespace hpc::repeat_detail
