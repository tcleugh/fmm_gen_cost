#pragma once

#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "network/type.hpp"
#include "network/network.hpp"

namespace FMM {
namespace ROUTING {

using FMM::NETWORK::Edge;
using FMM::NETWORK::EdgeIndex;
using FMM::NETWORK::NodeIndex;

inline constexpr EdgeIndex kNoEdge = std::numeric_limits<EdgeIndex>::max();

// ----------------------
// Public structs/classes
// ----------------------

struct Path {
  double total_cost = std::numeric_limits<double>::infinity();
  std::vector<EdgeIndex> edges;
  bool found = false;
};

// ------------------------------------
// LinkGraph: vertices are EdgeIndex-es
// ------------------------------------
class LinkGraph {
public:
  struct Arc { EdgeIndex to; double w; };

  explicit LinkGraph(const FMM::NETWORK::Network& net);

  const std::vector<Arc>& neighbors(EdgeIndex in_e) const;

  const FMM::NETWORK::Network& net() const { return net_; }

private:
  const FMM::NETWORK::Network& net_;

  std::vector<std::vector<EdgeIndex>> out_by_node_;
  std::vector<std::vector<Arc>> adj_;

  void buildNeighbors_(EdgeIndex in_e);
};


// ---------------------------
// Indexed min-heap (decrease-key)
// ---------------------------
class IndexedMinHeap {
public:
  static constexpr size_t npos = std::numeric_limits<size_t>::max();

  void ensure_size(size_t n);

  bool empty() const { return heap_.empty(); }
  size_t size() const { return heap_.size(); }

  void push_or_decrease(EdgeIndex v, double new_key);
  EdgeIndex pop_min();

  double key(EdgeIndex v) const { return key_[static_cast<size_t>(v)]; }

  void clear();

private:
  std::vector<EdgeIndex> heap_;
  std::vector<size_t>    pos_;
  std::vector<double>    key_;

  bool less_(EdgeIndex a, EdgeIndex b) const;
  void swap_nodes_(size_t i, size_t j);
  void sift_up_(size_t i);
  void sift_down_(size_t i);
};


// ---------------------------
// Sparse-reset Dijkstra state (epoch stamps)
// ---------------------------
struct DijkstraState {
  std::vector<double>    dist;
  std::vector<EdgeIndex> parent;
  std::vector<uint32_t>  seen;
  std::vector<uint32_t>  settled;
  uint32_t epoch = 1;

  // Goal tracking (epoch-stamped to avoid per-call allocation)
  std::vector<uint32_t>  goal_marker;
  uint32_t goal_epoch = 0;

  // Instrumentation counters (accumulated across calls, caller resets)
  size_t dijkstra_calls = 0;
  size_t nodes_explored = 0;
  std::vector<size_t> per_call_nodes;  // nodes explored per call, for distribution stats

  void ensure_size(size_t E);

  bool has_dist(EdgeIndex v) const;
  double get_dist(EdgeIndex v) const;

  void set_dist(EdgeIndex v, double d, EdgeIndex p);

  bool is_settled(EdgeIndex v) const;
  void mark_settled(EdgeIndex v);

  bool is_goal(EdgeIndex v) const;
  void mark_goal(EdgeIndex v);
  void clear_goal(EdgeIndex v);
  void next_goal_epoch();

  void next_epoch();
};


// ----------------------
// API
// ----------------------
Path reconstruct_path(const std::vector<EdgeIndex>& parent, EdgeIndex end_e);

/**
 * Dijkstra from start_e to all goal_edges.
 * Results are written into out[i] corresponding to goal_edges[i].
 * out must be pre-sized to goal_edges.size().
 *
 * upper_bound_factor: if > 0, after finding min_found_for_bound goals,
 * sets an upper bound of max_found_cost * factor. Nodes beyond this
 * bound are not explored; unfound goals get infinity.
 * min_found_for_bound: how many goals must be found before activating
 * the bound. Default: max(1, K/2).
 */
void shortest_edge_to_edges(
  const LinkGraph& G,
  DijkstraState& state,
  IndexedMinHeap& heap,
  EdgeIndex start_e,
  const std::vector<EdgeIndex>& goal_edges,
  std::vector<Path>& out,
  double upper_bound_factor = 0.0
);

} // namespace ROUTING
} // namespace FMM