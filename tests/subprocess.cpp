// SPDX-License-Identifier: MIT-0
// Re-execute this fixture to test process boundaries without invoking a shell.
#include "io.hpp"
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>

static void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}

int main(int argc, char **argv) {
  if (argc > 1) {
    const std::string mode = argv[1];
    if (mode == "arguments") {
      check(argc == 5 && std::string(argv[2]) == "space in filename" &&
                std::string(argv[3]) == "literal;$(command)*'\"" && std::string(argv[4]).empty(),
            "Spawn must preserve literal arguments, including an empty argument");
      const char *value = std::getenv("HPC_CONSENSUS_SPAWN_TEST");
      check(value && std::string(value) == "inherited", "Spawn must inherit environment");
      return 0;
    }
    if (mode == "failure")
      return 23;
    if (mode == "signal") {
      std::raise(SIGTERM);
      return 99;
    }
    return 98;
  }
  const std::string self = std::filesystem::absolute(argv[0]).string();
  check(setenv("HPC_CONSENSUS_SPAWN_TEST", "inherited", 1) == 0, "Set fixture environment");
  hpc::run({self, "arguments", "space in filename", "literal;$(command)*'\"", ""});
  // PATH lookup remains available to callers such as the samtools launcher.
  const std::string parent = std::filesystem::path(self).parent_path().string();
  const char *path = std::getenv("PATH");
  const std::string search = parent + ":" + (path ? path : "");
  check(setenv("PATH", search.c_str(), 1) == 0, "Set fixture PATH");
  hpc::run({std::filesystem::path(self).filename().string(), "arguments", "space in filename",
            "literal;$(command)*'\"", ""});

  bool missing = false;
  const std::string absent = self + "/missing-command";
  try {
    hpc::run({absent});
  } catch (const std::system_error &e) {
    // A path through a regular file fails with ENOTDIR on POSIX.
    missing =
        e.code().value() == ENOTDIR && std::string(e.what()).find(absent) != std::string::npos;
  }
  check(missing, "Launch failures must include the OS error and command");
  for (const auto &mode : {"failure", "signal"}) {
    bool failed = false;
    try {
      hpc::run({self, mode});
    } catch (const std::runtime_error &e) {
      const std::string message = e.what();
      const std::string detail =
          std::string(mode) == "failure" ? "exit status 23" : "signal " + std::to_string(SIGTERM);
      failed = message.find(self) != std::string::npos && message.find(detail) != std::string::npos;
    }
    check(failed, "Child failures must identify the command and exit status or signal");
  }
}
