#include "network/link_graph_routing.hpp"

#include <algorithm>   // std::reverse
#include <cassert>

// #include "util/debug.hpp" // optional logging

namespace FMM {
namespace ROUTING {

constexpr size_t IndexedMinHeap::npos;

// ----------------------
// LinkGraph
// ----------------------
LinkGraph::LinkGraph(const FMM::NETWORK::Network& net)
  : net_(net)
{
  const auto n_nodes = static_cast<size_t>(net_.get_node_count());
  out_by_node_.assign(n_nodes, {});

  const auto& edges = net_.get_edges();
  for (const auto& e : edges) {
    out_by_node_[e.source].push_back(e.index);
  }

  for (const auto& e : edges) {
    buildNeighbors_(e.index);
  }
}

const std::vector<LinkGraph::Arc>& LinkGraph::neighbors(EdgeIndex in_e) const {
  return adj_[in_e];
}

void LinkGraph::buildNeighbors_(EdgeIndex in_e) {
  const Edge& in = net_.get_edge(in_e);
  const NodeIndex join = in.target;

  std::vector<Arc> nbrs;
  for (EdgeIndex out_e : out_by_node_[join]) {
    if (net_.is_turn_banned(in_e, out_e)) continue;
    const Edge& out = net_.get_edge(out_e);
    nbrs.push_back(Arc{out_e, out.cost});
  }

  if (in_e >= adj_.size()) adj_.resize(net_.get_edge_count());
  adj_[in_e] = std::move(nbrs);
}


// ----------------------
// IndexedMinHeap
// ----------------------
void IndexedMinHeap::ensure_size(size_t n) {
  if (pos_.size() < n) {
    pos_.assign(n, npos);
    key_.assign(n, std::numeric_limits<double>::infinity());
  }
}

void IndexedMinHeap::push_or_decrease(EdgeIndex v, double new_key) {
  const size_t vi = static_cast<size_t>(v);
  if (pos_[vi] == npos) {
    key_[vi] = new_key;
    pos_[vi] = heap_.size();
    heap_.push_back(v);
    sift_up_(pos_[vi]);
  } else if (new_key < key_[vi]) {
    key_[vi] = new_key;
    sift_up_(pos_[vi]);
  }
}

EdgeIndex IndexedMinHeap::pop_min() {
  assert(!heap_.empty());
  EdgeIndex min_v = heap_.front();
  EdgeIndex last  = heap_.back();
  heap_.pop_back();

  pos_[static_cast<size_t>(min_v)] = npos;

  if (!heap_.empty()) {
    heap_[0] = last;
    pos_[static_cast<size_t>(last)] = 0;
    sift_down_(0);
  }
  return min_v;
}

void IndexedMinHeap::clear() {
  for (EdgeIndex v : heap_) pos_[static_cast<size_t>(v)] = npos;
  heap_.clear();
}

bool IndexedMinHeap::less_(EdgeIndex a, EdgeIndex b) const {
  return key_[static_cast<size_t>(a)] < key_[static_cast<size_t>(b)];
}

void IndexedMinHeap::swap_nodes_(size_t i, size_t j) {
  EdgeIndex vi = heap_[i];
  EdgeIndex vj = heap_[j];
  heap_[i] = vj; heap_[j] = vi;
  pos_[static_cast<size_t>(vi)] = j;
  pos_[static_cast<size_t>(vj)] = i;
}

void IndexedMinHeap::sift_up_(size_t i) {
  while (i > 0) {
    size_t p = (i - 1) >> 1;
    if (!less_(heap_[i], heap_[p])) break;
    swap_nodes_(i, p);
    i = p;
  }
}

void IndexedMinHeap::sift_down_(size_t i) {
  const size_t n = heap_.size();
  while (true) {
    size_t l = (i << 1) + 1;
    if (l >= n) break;
    size_t r = l + 1;
    size_t m = (r < n && less_(heap_[r], heap_[l])) ? r : l;
    if (!less_(heap_[m], heap_[i])) break;
    swap_nodes_(i, m);
    i = m;
  }
}


// ----------------------
// DijkstraState
// ----------------------
void DijkstraState::ensure_size(size_t E) {
  if (dist.size() < E) {
    dist.resize(E);
    parent.resize(E);
    seen.assign(E, 0);
    settled.assign(E, 0);
    goal_marker.assign(E, 0);
    epoch = 1;
    goal_epoch = 0;
  }
}

bool DijkstraState::has_dist(EdgeIndex v) const {
  return seen[v] == epoch;
}

double DijkstraState::get_dist(EdgeIndex v) const {
  return has_dist(v) ? dist[v] : std::numeric_limits<double>::infinity();
}

void DijkstraState::set_dist(EdgeIndex v, double d, EdgeIndex p) {
  if (seen[v] != epoch) {
    seen[v] = epoch;
    // settled has its own epoch marking
  }
  dist[v] = d;
  parent[v] = p;
}

bool DijkstraState::is_settled(EdgeIndex v) const {
  return settled[v] == epoch;
}

void DijkstraState::mark_settled(EdgeIndex v) {
  settled[v] = epoch;
}

bool DijkstraState::is_goal(EdgeIndex v) const {
  return goal_marker[v] == goal_epoch;
}

void DijkstraState::mark_goal(EdgeIndex v) {
  goal_marker[v] = goal_epoch;
}

void DijkstraState::clear_goal(EdgeIndex v) {
  goal_marker[v] = 0;  // any value != goal_epoch
}

void DijkstraState::next_goal_epoch() {
  ++goal_epoch;
  if (goal_epoch == 0) {
    std::fill(goal_marker.begin(), goal_marker.end(), 0);
    goal_epoch = 1;
  }
}

void DijkstraState::next_epoch() {
  ++epoch;
  if (epoch == 0) {
    std::fill(seen.begin(), seen.end(), 0);
    std::fill(settled.begin(), settled.end(), 0);
    epoch = 1;
  }
}


// ----------------------
// Path reconstruction
// ----------------------
Path reconstruct_path(const std::vector<EdgeIndex>& parent, EdgeIndex end_e) {
  Path P;
  if (parent[end_e] == kNoEdge) return P;

  std::vector<EdgeIndex> seq;
  EdgeIndex cur = end_e;
  while (true) {
    seq.push_back(cur);
    EdgeIndex par = parent[cur];
    if (par == kNoEdge) break;
    cur = par;
  }
  std::reverse(seq.begin(), seq.end());
  P.edges = std::move(seq);
  P.found = true;
  return P;
}

void shortest_edge_to_edges(
  const LinkGraph& G,
  DijkstraState& state,
  IndexedMinHeap& heap,
  EdgeIndex start_e,
  const std::vector<EdgeIndex>& goal_edges,
  std::vector<Path>& out,
  double upper_bound_factor
) {
  const size_t E = static_cast<size_t>(G.net().get_edge_count());
  const size_t K = goal_edges.size();
  state.ensure_size(E);
  state.next_epoch();
  heap.ensure_size(E);
  heap.clear();
  ++state.dijkstra_calls;
  size_t call_nodes = 0;

  // Reset output slots
  for (size_t i = 0; i < K; ++i) {
    out[i] = Path{};
  }

  if (K == 0) return;

  // Epoch-stamped goal tracking: O(1) lookup, no per-call allocation
  state.next_goal_epoch();
  size_t remaining = 0;
  size_t unique_goals = 0;

  for (size_t i = 0; i < K; ++i) {
    EdgeIndex g = goal_edges[i];
    if (!state.is_goal(g)) {
      state.mark_goal(g);
      ++remaining;
      ++unique_goals;
    }
  }

  // Upper bound state
  const bool use_bound = upper_bound_factor > 0.0;
  const size_t min_found_for_bound =
    use_bound ? static_cast<size_t>(std::max(1, static_cast<int>(unique_goals) / 2)) : 0;
  size_t goals_found = 0;
  double max_found_cost = 0.0;
  double bound = std::numeric_limits<double>::infinity();

  if (state.is_goal(start_e)) {
    // Fill all output slots that match start_e
    for (size_t i = 0; i < K; ++i) {
      if (goal_edges[i] == start_e) {
        out[i].found = true;
        out[i].total_cost = 0.0;
      }
    }
    state.clear_goal(start_e);
    --remaining;
    // Don't count cost-0 start toward bound (would make bound = 0)
    if (remaining == 0) {
      state.per_call_nodes.push_back(0);
      return;
    }
  }

  state.set_dist(start_e, 0.0, kNoEdge);
  heap.push_or_decrease(start_e, 0.0);

  constexpr double eps = 1e-15;

  while (!heap.empty() && remaining > 0) {
    const EdgeIndex u = heap.pop_min();
    if (state.is_settled(u)) continue;
    state.mark_settled(u);
    ++state.nodes_explored;
    ++call_nodes;

    const double du = state.get_dist(u);

    // Upper bound check: stop if we've exceeded the bound
    if (use_bound && du > bound) break;

    if (state.is_goal(u)) {
      Path p = reconstruct_path(state.parent, u);
      p.total_cost = du;
      p.found = true;
      // Fill all output slots that match u
      for (size_t i = 0; i < K; ++i) {
        if (goal_edges[i] == u) {
          out[i] = p;
        }
      }
      state.clear_goal(u);
      --remaining;

      // Update adaptive bound
      if (use_bound && du > 0.0) {
        ++goals_found;
        if (du > max_found_cost) max_found_cost = du;
        if (goals_found >= min_found_for_bound) {
          bound = max_found_cost * upper_bound_factor;
        }
      }
    }

    const auto& nbrs = G.neighbors(u);
    for (size_t i = 0; i < nbrs.size(); ++i) {
      const EdgeIndex v = nbrs[i].to;
      const double nd = du + nbrs[i].w;

      if (nd + eps < state.get_dist(v)) {
        state.set_dist(v, nd, u);
        heap.push_or_decrease(v, nd);
      }
    }
  }

  state.per_call_nodes.push_back(call_nodes);
}

} // namespace ROUTING
} // namespace FMM