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
    epoch = 1;
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

std::unordered_map<EdgeIndex, Path> shortest_edge_to_edges(
  const LinkGraph& G,
  DijkstraState& state,
  IndexedMinHeap& heap,
  EdgeIndex start_e,
  const std::vector<EdgeIndex>& goal_edges
) {
  const size_t E = static_cast<size_t>(G.net().get_edge_count());
  state.ensure_size(E);
  state.next_epoch();
  heap.ensure_size(E);
  ++state.dijkstra_calls;

  std::unordered_map<EdgeIndex, Path> out;

  if (goal_edges.empty()) return out;

  std::vector<uint8_t> is_goal(E, 0);
  size_t remaining = 0;

  for (EdgeIndex g : goal_edges) {
    if (!is_goal[g]) {
      is_goal[g] = 1;
      ++remaining;
      out.emplace(g, Path{});
    }
  }

  if (is_goal[start_e]) {
    Path P;
    P.found = true;
    P.total_cost = 0.0;
    out[start_e] = std::move(P);

    is_goal[start_e] = 0;
    --remaining;
    if (remaining == 0) return out;
  }

  state.set_dist(start_e, 0.0, kNoEdge);
  heap.push_or_decrease(start_e, 0.0);

  constexpr double eps = 1e-15;

  while (!heap.empty() && remaining > 0) {
    const EdgeIndex u = heap.pop_min();
    if (state.is_settled(u)) continue;
    state.mark_settled(u);
    ++state.nodes_explored;

    const double du = state.get_dist(u);

    if (is_goal[u]) {
      Path& slot = out[u];
      if (!slot.found || du < slot.total_cost) {
        slot = reconstruct_path(state.parent, u);
        slot.total_cost = du;
        slot.found = true;
      }
      is_goal[u] = 0;
      --remaining;
    }

    const auto& nbrs = G.neighbors(u);
    for (size_t i = 0; i < nbrs.size(); ++i) {
      const EdgeIndex v = nbrs[i].to;
      const double nd = du + nbrs[i].w;

      if (nd + eps < state.get_dist(v)) {
        state.set_dist(v, nd, u);
        heap.push_or_decrease(v, nd);  // IMPORTANT: (v, key)
      }
    }
  }

  return out;
}

} // namespace ROUTING
} // namespace FMM