// SPDX-License-Identifier: MIT-0
// Local diagnostic only: same read chunks, HPC evidence, raw evidence, and
// Verkko's actual alignment graph. Not the full utgcns template/placement path.
#include "core.hpp"
#include "Alignment.H"
#include "AlnGraphBoost.H"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace hpc;
struct Read {
  std::string name, tier, sequence, query, target;
  int hpc_length = -1;
  std::string raw_allele;
};

// Represent mismatches as insertion/deletion alternatives for the DAG, then
// move equivalent gaps toward the right. This implements the normalization
// convention of unitigConsensus::alignEdLib, using independently held strings.
static void normalize(Read &read, const std::string &ref) {
  auto ops = align(read.sequence, ref, 0.5);
  size_t q = 0, r = 0;
  for (char op : ops) {
    if (op == 'X') {
      read.query += read.sequence[q++];
      read.target += '-';
      read.query += '-';
      read.target += ref[r++];
    } else {
      read.query += op == 'D' ? '-' : read.sequence[q++];
      read.target += op == 'I' ? '-' : ref[r++];
    }
  }
  for (size_t i = 0; i + 1 < read.query.size(); ++i)
    for (bool query_gap : {false, true}) {
      auto &gapped = query_gap ? read.query : read.target;
      const auto &opposite = query_gap ? read.target : read.query;
      if (gapped[i] != '-')
        continue;
      auto next = gapped.find_first_not_of('-', i + 1);
      if (next != std::string::npos && gapped[next] == opposite[i])
        std::swap(gapped[i], gapped[next]);
    }
  size_t out = 0;
  for (size_t i = 0; i < read.query.size(); ++i)
    if (read.query[i] != '-' || read.target[i] != '-') {
      read.query[out] = read.query[i];
      read.target[out++] = read.target[i];
    }
  read.query.resize(out);
  read.target.resize(out);
  auto ungap = [](std::string text) {
    text.erase(std::remove(text.begin(), text.end(), '-'), text.end());
    return text;
  };
  if (ungap(read.query) != read.sequence || ungap(read.target) != ref)
    throw std::runtime_error("Normalization changed an input sequence");
}

static std::string graph(const std::string &ref, const std::vector<const Read *> &reads) {
  AlnGraphBoost graph(ref);
  for (auto read : reads) {
    dagAlignment a;
    a.start = 1;
    a.end = ref.size();
    a.length = read->query.size();
    a.qstr = new char[a.length + 1]{};
    a.tstr = new char[a.length + 1]{};
    std::copy(read->query.begin(), read->query.end(), a.qstr);
    std::copy(read->target.begin(), read->target.end(), a.tstr);
    graph.addAln(a);
  }
  graph.mergeNodes();
  std::vector<uint32_t> mapping(ref.size() + 1);
  return graph.consensusNoSplit(0, mapping.data(), ref.size());
}

int main(int argc, char **argv) {
  try {
    if (argc != 5)
      throw std::runtime_error("Usage: hp-graph-probe reference.txt reads.tsv run_begin run_end");
    std::ifstream ref_in(argv[1]), read_in(argv[2]);
    std::string ref;
    ref_in >> ref;
    size_t begin = std::stoull(argv[3]), end = std::stoull(argv[4]);
    if (begin == 0 || begin >= end || end >= ref.size())
      throw std::runtime_error("Invalid run");
    auto reference = encode(ref);
    auto run = std::lower_bound(reference.offsets.begin(), reference.offsets.end(), begin) -
               reference.offsets.begin();
    if (reference.offsets.at(run) != begin || reference.offsets.at(run + 1) != end)
      throw std::runtime_error("Not a complete reference run");
    std::vector<Read> reads;
    std::set<std::string> names;
    Read read;
    while (read_in >> read.name >> read.tier >> read.sequence) {
      if (!names.insert(read.name).second ||
          (read.tier != "LA" && read.tier != "UL" && read.tier != "LAO"))
        throw std::runtime_error("Expected unique read names and LA/UL/LAO diagnostic tiers");
      try {
        auto query = encode(read.sequence);
        auto ops = align(query.bases, reference.bases, 0.5);
        auto eligible = local_hpc_votes(ops, 500, 0.15);
        size_t q = 0, r = 0;
        for (size_t i = 0; i < ops.size(); ++i) {
          if (r == size_t(run) && ops[i] == '=' && (eligible.empty() || eligible[i]))
            read.hpc_length = query.lengths[q];
          q += ops[i] != 'D';
          r += ops[i] != 'I';
        }
        normalize(read, ref);
        // Count query bases between the immediate reference flanks. Include
        // interrupted alleles in the report; never silently call them pure runs.
        r = 0;
        bool inside = false, left = false, right = false;
        for (size_t i = 0; i < read.target.size(); ++i) {
          if (read.target[i] != '-' && r == end) {
            right = read.query[i] == ref[end];
            break;
          }
          if (inside && read.query[i] != '-')
            read.raw_allele += read.query[i];
          if (read.target[i] != '-' && r == begin - 1) {
            left = read.query[i] == ref[begin - 1];
            inside = true;
          }
          r += read.target[i] != '-';
        }
        if (!left || !right)
          read.raw_allele = "UNANCHORED";
        std::cout << "READ\t" << read.name << '\t' << read.tier << '\t' << read.hpc_length << '\t'
                  << (read.raw_allele.empty() ? "EMPTY" : read.raw_allele) << '\n';
        reads.push_back(std::move(read));
      } catch (const std::runtime_error &error) {
        // Bad local chunks are excluded identically from both estimators.
        // Keep an explicit audit record; do not suppress unrelated failures.
        if (std::string(error.what()) != "Alignment exceeds maximum error rate")
          throw;
        std::cout << "REJECT\t" << read.name << '\t' << error.what() << '\n';
      }
      read = {};
    }
    // An all-unique control tests evidence mixing, without pretending to
    // reproduce utgcns's whole-read promotion decisions on this local slice.
    std::vector<const Read *> all;
    for (const auto &r : reads)
      if (r.hpc_length >= 0) {
        all.push_back(&r);
        std::cout << "SELECT\tall_unique\t" << r.name << '\n';
      }
    if (!all.empty())
      std::cout << "GRAPH\tall_unique\t" << graph(ref, all) << '\n';
    std::vector<const Read *> without_ul;
    for (auto r : all)
      if (r->tier != "UL") {
        without_ul.push_back(r);
        std::cout << "SELECT\tLA_LAO_unique\t" << r->name << '\n';
      }
    if (!without_ul.empty())
      std::cout << "GRAPH\tLA_LAO_unique\t" << graph(ref, without_ul) << '\n';
    for (const auto &order : {std::vector<std::string>{"LA", "UL", "LAO"},
                              std::vector<std::string>{"LA", "LAO", "UL"}}) {
      std::map<std::string, size_t> depth;
      for (const auto &r : reads)
        if (r.hpc_length >= 0)
          ++depth[r.tier];
      size_t cutoff = order.size() - 1;
      // Match the caller's sparse fallback: first nonempty tier if none
      // individually reaches five. This diagnostic uses unique reads only.
      for (size_t t = 0; t < order.size(); ++t)
        if (depth[order[t]]) {
          cutoff = t;
          break;
        }
      for (size_t t = 0; t < order.size(); ++t)
        if (depth[order[t]] >= 5) {
          cutoff = t;
          break;
        }
      for (bool retain : {false, true}) {
        std::vector<const Read *> chosen;
        std::string label = order[1] + (retain ? "_additive" : "_exclusive");
        for (const auto &r : reads) {
          auto rank = std::find(order.begin(), order.end(), r.tier) - order.begin();
          if (r.hpc_length >= 0 && (retain ? rank <= cutoff : rank == cutoff)) {
            chosen.push_back(&r);
            std::cout << "SELECT\t" << label << '\t' << r.name << '\n';
          }
        }
        if (!chosen.empty())
          std::cout << "GRAPH\t" << label << '\t' << graph(ref, chosen) << '\n';
      }
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
