// SPDX-License-Identifier: MIT-0
// Diagnostic fixture extraction using native readers to discover record bounds.
// Copy serialized reads verbatim: do not rebuild names, clipping, or sequences.
#include "sqStore.H"
#include "tgTig.H"
#include "files.H"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <iostream>

static void check(bool ok, const char *msg) {
  if (!ok)
    throw std::runtime_error(msg);
}
struct Region {
  std::string prefix;
  int64_t begin, end;
};
int main(int argc, char **argv) {
  try {
    check(argc == 5, "Usage: extract-regions input.cnspack tigID template_offset "
                     "regions.tsv\nRows: output_prefix HPC_begin HPC_end (zero-based half-open)");
    uint32_t target = std::stoul(argv[2]);
    int64_t offset = std::stoll(argv[3]);
    std::ifstream regions(argv[4]);
    check(bool(regions), "Cannot open region list");
    std::vector<Region> windows;
    Region window;
    while (regions >> window.prefix >> window.begin >> window.end) {
      check(window.begin >= 0 && window.end > window.begin, "Invalid region");
      for (const auto &suffix : {".cnspack", ".placements.tsv", ".span.tsv"})
        check(!std::filesystem::exists(window.prefix + suffix), "Output already exists");
      windows.push_back(window);
    }
    check(regions.eof() && !windows.empty(), "Malformed/empty region list");
    std::vector<tgPosition> positions;
    std::map<uint32_t, unsigned> occurrences;
    std::map<uint32_t, unsigned> declared;
    std::map<uint32_t, std::pair<uint64_t, uint64_t>> records;
    std::set<uint32_t> wanted;
    std::map<uint32_t, std::string> names;
    readBuffer input(argv[1]);
    sqRead_setDefaultVersion(sqRead_raw | sqRead_normal);
    bool found = false;
    while (!input.eof()) {
      tgTig tig;
      if (!tig.loadFromBuffer(&input))
        break;
      std::set<std::tuple<uint32_t, int32_t, int32_t, bool, int32_t, int32_t>> seen;
      for (uint32_t i = 0; i < tig.numberOfChildren(); ++i) {
        auto p = tig.getChild(i);
        if (!seen.emplace(p->ident(), p->min(), p->max(), p->isReverse(), p->_askip, p->_bskip)
                 .second)
          continue;
        ++occurrences[p->ident()];
        declared[p->ident()] = std::max(declared[p->ident()], unsigned(p->_numPlacement));
        if (tig.tigID() == target) {
          positions.push_back(*p);
          for (const auto &w : windows)
            if (int64_t(p->min()) - offset < w.end && int64_t(p->max()) - offset > w.begin)
              wanted.insert(p->ident());
        }
      }
      if (tig.tigID() == target) {
        check(!found, "Repeated target tig");
        found = true;
      }
      sqRead read;
      for (uint32_t i = 0; i <= tig.numberOfChildren(); ++i) {
        auto begin = input.tell();
        sqStore::sqStore_loadReadFromBuffer(&input, &read);
        if (tig.tigID() == target && wanted.count(read.sqRead_readID())) {
          records.emplace(read.sqRead_readID(), std::make_pair(begin, input.tell() - begin));
          names.emplace(read.sqRead_readID(), read.sqRead_name());
        }
      }
    }
    check(found, "Target tig absent");
    std::ifstream source(argv[1], std::ios::binary);
    std::vector<char> buffer(1024 * 1024);
    for (const auto &w : windows) {
      std::vector<tgPosition> selected;
      int64_t begin = w.begin, end = w.end;
      for (const auto &p : positions)
        if (int64_t(p.min()) - offset < w.end && int64_t(p.max()) - offset > w.begin) {
          selected.push_back(p);
          begin = std::min(begin, std::max<int64_t>(0, int64_t(p.min()) - offset));
          end = std::max(end, int64_t(p.max()) - offset);
        }
      check(!selected.empty(), "No overlapping placements");
      tgTig subset;
      subset._tigID = target;
      std::ofstream meta(w.prefix + ".placements.tsv");
      meta << "read_id\tname\toriginal_hpc_begin\toriginal_hpc_end\tplacements\tskipConsensus\n";
      for (auto p : selected) {
        // Retain ambiguity even when alternative placements fall outside the
        // fixture. The native field saturates at 31; preserve unknown (zero).
        if (p._numPlacement)
          p._numPlacement = std::min(31u, std::max(declared[p.ident()], occurrences[p.ident()]));
        meta << p.ident() << '\t' << names.at(p.ident()) << '\t' << int64_t(p.min()) - offset
             << '\t' << int64_t(p.max()) - offset << '\t' << p._numPlacement << '\t'
             << p.skipConsensus() << '\n';
        p.setMinMax(p.min() - begin, p.max() - begin);
        check(p.min() >= 0, "Negative rebased layout coordinate");
        subset._layoutLen = std::max(subset._layoutLen, uint32_t(p.max()));
        *subset.addChild() = p;
      }
      writeBuffer out((w.prefix + ".cnspack").c_str(), "w");
      subset.saveToBuffer(&out);
      auto copy = [&](uint32_t id) {
        auto [pos, len] = records.at(id);
        source.clear();
        source.seekg(pos);
        while (len) {
          size_t n = std::min<uint64_t>(len, buffer.size());
          source.read(buffer.data(), n);
          check(size_t(source.gcount()) == n, "Truncated native read record");
          out.write(buffer.data(), n);
          len -= n;
        }
      };
      copy(subset.getChild(0)->ident()); // Native redundant leading read record.
      for (uint32_t i = 0; i < subset.numberOfChildren(); ++i)
        copy(subset.getChild(i)->ident());
      std::ofstream span(w.prefix + ".span.tsv");
      span << "template_begin\ttemplate_end\trequested_begin\trequested_end\toffset\n"
           << begin << '\t' << end << '\t' << w.begin << '\t' << w.end << '\t' << offset << '\n';
      std::cerr << w.prefix << ": " << selected.size() << " placements; HPC " << begin << "-" << end
                << '\n';
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
