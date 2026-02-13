#pragma once
#include <unordered_map>
#include <cstdint>
#include "network/link_graph_routing.hpp"

namespace FMM::MM {

using FMM::NETWORK::EdgeIndex;
using FMM::ROUTING::Path;
using FMM::ROUTING::LinkGraph;
using FMM::ROUTING::DijkstraState;
using FMM::ROUTING::IndexedMinHeap;
using FMM::ROUTING::shortest_edge_to_edges;

// Pack two 32-bit EdgeIndex values into one 64-bit key.
// If EdgeIndex is already 32-bit-ish, this is fast and simple.
static inline uint64_t pack_key(EdgeIndex s, EdgeIndex t) {
  return (uint64_t(uint32_t(s)) << 32) | uint64_t(uint32_t(t));
}

class RoutingCache {
public:
  void clear() { cache_.clear(); }
  void reserve(size_t n) { cache_.reserve(n); }

  // Get routing results for a single source to many targets.
  // Only runs Dijkstra for "missing" targets.
  std::unordered_map<EdgeIndex, Path> get_or_compute(
      const LinkGraph& G,
      DijkstraState& state,
      IndexedMinHeap& heap,
      EdgeIndex source,
      const std::vector<EdgeIndex>& targets
  ) {
    std::unordered_map<EdgeIndex, Path> out;
    out.reserve(targets.size());

    std::vector<EdgeIndex> missing;
    missing.reserve(targets.size());

    for (EdgeIndex t : targets) {
      const uint64_t k = pack_key(source, t);
      auto it = cache_.find(k);
      if (it != cache_.end()) {
        out.emplace(t, it->second);
      } else {
        // placeholder; we’ll fill after routing
        out.emplace(t, Path{});
        missing.push_back(t);
      }
    }

    if (!missing.empty()) {
      auto fresh = shortest_edge_to_edges(G, state, heap, source, missing);

      // shortest_edge_to_edges returns an entry for each goal (found or not)
      for (auto& [t, p] : fresh) {
        const uint64_t k = pack_key(source, t);
        cache_[k] = p;
        out[t] = p;
      }
    }
    return out;
  }

  // Convenience: single target
  Path get_or_compute_one(
      const LinkGraph& G,
      DijkstraState& state,
      IndexedMinHeap& heap,
      EdgeIndex source,
      EdgeIndex target)
  {
    const uint64_t k = pack_key(source, target);
    auto it = cache_.find(k);
    if (it != cache_.end()) return it->second;

    std::vector<EdgeIndex> goals{target};
    auto fresh = shortest_edge_to_edges(G, state, heap, source, goals);
    Path p = fresh[target];
    cache_[k] = p;
    return p;
  }

private:
  std::unordered_map<uint64_t, Path> cache_;
};

} // namespace FMM::MM