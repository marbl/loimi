// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace hpc {
// Failure-only rescue anchors in HPC coordinates. Each exact chunk must occur
// once in the template; different chunks from one read are NOT independent votes.
std::vector<int64_t> internal_offsets(const std::string &query, const std::string &reference,
                                      int64_t layout_begin);
// Require >=75% of informative reads and reject a strong distant competitor.
// Empty/ambiguous evidence throws. Returns one-read-one-vote lower median.
int64_t internal_offset_cluster(const std::vector<std::vector<int64_t>> &reads, size_t required,
                                int64_t tolerance = 20);
} // namespace hpc
