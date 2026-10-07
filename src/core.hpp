#pragma once
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace hpc {
struct Rle { std::string bases; std::vector<uint32_t> lengths; std::vector<size_t> offsets; };
Rle encode(const std::string &s);
std::string reverse_complement(std::string s);
// Each operation consumes one query and/or reference base: = X I D.
std::string align(const std::string &query, const std::string &reference, double max_error);
std::string cigar(const std::string &ops);
std::string expand_cigar(const std::string &packed);
struct BoundaryAlignment { std::string ops; size_t query_begin=0,query_end=0; };
BoundaryAlignment align_boundary(const std::string &query,const std::string &reference,
                                 bool left,bool right,double max_error);
// Return the unique best full-query infix placement, or -1 if unresolved.
int64_t locate(const std::string &query, const std::string &reference, double max_error);
int64_t offset_cluster(std::vector<int64_t> estimates, size_t required=3, int64_t tolerance=20);
struct Tier { std::string label, pattern; };
int tier_of(const std::string &name, const std::vector<Tier> &tiers);
struct Call { uint32_t length=1, depth=0, peak=0; int tier=-1; };
using Histogram = std::map<uint32_t,uint32_t>;
Call call(const std::vector<Histogram> &tiers, uint32_t min_support);
struct Edit { size_t begin=0, end=0; std::string alt; bool operator<(const Edit &o) const; };
struct Alignment {
  uint32_t id=0; std::string name, sequence, ops; size_t begin=0, end=0;
  size_t clip_left=0, clip_right=0; int tier=0; bool reverse=false;
  bool packed_ops=false;
  std::shared_ptr<const std::string> shared_sequence;
  const std::string &bases() const {return shared_sequence?*shared_sequence:sequence;}
};
std::vector<Edit> events(const Alignment &a, const std::string &reference);
struct Accepted { Edit edit; int tier; uint32_t depth, support; };
std::vector<Accepted> polish(const std::string &ref, const std::vector<Alignment> &alns,
                            size_t tiers, uint32_t support, double fraction, size_t max_indel);
std::vector<Accepted> polish(const std::string &ref, const std::vector<const Alignment*> &alns,
                            size_t tiers, uint32_t support, double fraction, size_t max_indel);
std::string apply(const std::string &ref, const std::vector<Accepted> &edits);
size_t lift(size_t p, const std::vector<Accepted> &edits);
std::string md(const Alignment &a, const std::string &ref, size_t &nm);
}
