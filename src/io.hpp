// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstdint>
#include <map>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace hpc {
// Plain FASTA, uppercase bases, first header token as ID. Reject duplicate IDs.
// An empty optional path returns an empty map.
std::map<std::string, std::string> fasta(const std::string &path);
// Spawn literal argv (never a shell), inheriting PATH/environment and streams.
// Avoids forking the large consensus address space; waits and throws on failure.
void run(std::vector<std::string> args);
// Wrapped at 80 bases; sequence is borrowed, not copied.
void write_fasta(std::ostream &output, const std::string &name, std::string_view sequence);
} // namespace hpc
