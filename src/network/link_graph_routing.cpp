
#pragma once
#include <vector>
#include <queue>
#include <limits>
#include <functional>
#include <utility>
#include <cassert>
#include <cstdint>
#include <unordered_map>

#include "network/type.hpp"
#include "network/network.hpp"
#include "util/debug.hpp"

namespace FMM {
namespace ROUTING {

using FMM::NETWORK::Edge;
using FMM::NETWORK::EdgeIndex;
using FMM::NETWORK::NodeIndex;

static constexpr EdgeIndex kNoEdge = std::numeric_limits<EdgeIndex>::max();

// ------------------------------
// Configurable edge cost functor
// ------------------------------
struct LengthCost {
  double operator()(const Edge& e) const noexcept { return e.length; }
};

struct WeightedCost {
  double operator()(const Edge& e) const noexcept { return e.length * e.weight; }
};


// ------------------------------------
// LinkGraph: vertices are EdgeIndex-es
// ------------------------------------
template<class CostFn = WeightedCost>
class LinkGraph {
public:
  struct Arc { EdgeIndex to; double w; };

  struct BuildOptions {
    bool allow_u_turns = true;    // drop turns where M.source==L.target && M.target==L.source
  };

  explicit LinkGraph(const FMM::NETWORK::Network& net,
                     BuildOptions opt = {},
                     CostFn cost_fn = CostFn())
    : net_(net), opt_(opt), cost_fn_(cost_fn)
  {
    initOutgoingByNode_();
  }

  // Adjacency on demand (built once per in-edge)
  const std::vector<Arc>& neighbors(EdgeIndex in_e) {
    if (in_e >= adj_.size()) adj_.resize(net_.get_edge_count());
    if (!built_.count(in_e)) buildNeighbors_(in_e);
    return adj_[in_e];
  }

  const FMM::NETWORK::Network& net() const { return net_; }

private:
  const FMM::NETWORK::Network& net_;
  BuildOptions opt_;
  CostFn cost_fn_;

  // node -> list of outgoing edges (EdgeIndex)
  std::vector<std::vector<EdgeIndex>> out_by_node_;

  // cache: in-edge -> arcs
  std::vector<std::vector<Arc>> adj_;
  std::unordered_set<EdgeIndex> built_;

  void initOutgoingByNode_() {
    const auto n_nodes = static_cast<size_t>(net_.get_node_count());
    out_by_node_.assign(n_nodes, {});
    const auto& edges = net_.get_edges();
    out_by_node_.shrink_to_fit();

    // Map: NodeIndex is already dense [0, num_vertices)
    for (const auto& e : edges) {
      if (e.source >= out_by_node_.size())
        throw std::runtime_error("Bad NodeIndex in edge.source");
      out_by_node_[e.source].push_back(e.index);
    }
  }

  void buildNeighbors_(EdgeIndex in_e) {
    const Edge& in = net_.get_edge(in_e);
    const NodeIndex join = in.target; // turns occur at the end node of in-edge
    std::vector<Arc> nbrs;

    // enumerate all links that start at the join node
    if (join < out_by_node_.size()) {
      for (EdgeIndex out_e : out_by_node_[join]) {
        const Edge& out = net_.get_edge(out_e);

        // optional U-turn filter
        bool u_turn = (in.source == out.target) && (in.target == out.source);
        if (!opt_.allow_u_turns && u_turn) continue;

        if (net_.is_turn_banned(in_e, out_e))
            SPDLOG_TRACE("Turn banned {} -> {}", in_e, out_e);

        // turn bans
        if (net_.is_turn_banned(in_e, out_e)) continue;

        // weight = cost(out-link)
        nbrs.push_back(Arc{out_e, cost_fn_(out)});
      }
    }

    if (in_e >= adj_.size()) adj_.resize(net_.get_edge_count());
    adj_[in_e] = std::move(nbrs);
    built_.insert(in_e);
  }
};

// -------------------------------------
// Dijkstra over the link-as-vertex graph
// -------------------------------------
struct DijResult {
  std::vector<double> dist;        // size = |E|, INF if unreachable
  std::vector<EdgeIndex> parent;   // predecessor edge (UINT_MAX sentinel)
};

template<class CostFn = WeightedCost>
DijResult dijkstra(LinkGraph<CostFn>& G,
                   const std::vector<std::pair<EdgeIndex, double>>& sources, // (edge, initCost)
                   const std::function<bool(EdgeIndex)>& isTarget = nullptr,
                   EdgeIndex stopAt = kNoEdge) {

  const auto M = static_cast<size_t>(G.net().get_edge_count());
  DijResult R;
  R.dist.assign(M, std::numeric_limits<double>::infinity());
  R.parent.assign(M, std::numeric_limits<EdgeIndex>::max());

  using QItem = std::pair<double, EdgeIndex>; // (dist, edge)
  std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;

  // seed
  for (size_t i = 0; i < sources.size(); ++i) {
    const EdgeIndex e = sources[i].first;
    const double d0   = sources[i].second;
    R.dist[e] = d0;
    R.parent[e] = std::numeric_limits<EdgeIndex>::max();
    pq.emplace(d0, e);
  }  
    // ---- SPECIAL CASE: trivial target found among seeds ----
    if (isTarget) {
        for (size_t i = 0; i < sources.size(); ++i) {
            const EdgeIndex e = sources[i].first;
            const double d0   = sources[i].second;

            if (isTarget(e)) {
                DijResult R2;
                R2.dist.assign(M, std::numeric_limits<double>::infinity());
                R2.parent.assign(M, std::numeric_limits<EdgeIndex>::max());
                R2.dist[e] = d0;   // typically d0 = 0
                return R2;         // no edges explored → preserves turn bans
            }
        }
    }

  std::vector<char> settled(M, 0);

  while (!pq.empty()) {
    const QItem top = pq.top(); pq.pop();
    const double    du = top.first;
    const EdgeIndex u  = top.second;
    if (settled[u]) continue;
    settled[u] = 1;
    
    // Early stop iff:
    //  - caller provided a concrete stop edge (stopAt != kNoEdge), AND
    //  - a target predicate exists (isTarget), AND
    //  - we popped that exact edge from the queue (optimal now), AND
    //  - it satisfies the target predicate
    if (stopAt != kNoEdge && isTarget && isTarget(u) && u == stopAt) break;

    for (const auto& arc : G.neighbors(u)) {
      EdgeIndex v = arc.to;
      double nd = du + arc.w;
      if (nd + 1e-15 < R.dist[v]) {
        R.dist[v] = nd;
        R.parent[v] = u;
        pq.emplace(nd, v);
      }
    }
  }
  return R;
}

// ----------------------
// Path reconstruction
// ----------------------
struct Path {
  double total_cost = std::numeric_limits<double>::infinity();
  std::vector<EdgeIndex> edges; // edge sequence
  std::vector<NodeIndex> nodes; // optional (filled if requested)
  bool found = false;
};

inline Path reconstruct_path(const FMM::NETWORK::Network& net,
                             const std::vector<EdgeIndex>& parent,
                             EdgeIndex end_e,
                             bool with_nodes = true) {
  Path P;
  if (end_e >= parent.size()) return P;
  if (parent[end_e] == std::numeric_limits<EdgeIndex>::max()) {
    // could still be a seed if the end_e was seeded; treat as a 1-edge path
    // will be validated by caller via dist
  }
  std::vector<EdgeIndex> seq;
  EdgeIndex cur = end_e;
  // climb parents; stop when parent is UINT_MAX
  while (true) {
    seq.push_back(cur);
    EdgeIndex par = parent[cur];
    if (par == std::numeric_limits<EdgeIndex>::max()) break;
    cur = par;
  }
  std::reverse(seq.begin(), seq.end());
  P.edges = seq;
  P.found = true;

  if (with_nodes && !seq.empty()) {
    const auto& E0 = net.get_edge(seq.front());
    P.nodes.push_back(E0.source);
    for (auto e : seq) P.nodes.push_back(net.get_edge(e).target);
  }
  return P;
}

// ---------------------------------------
// Convenience wrappers for common tasks
// ---------------------------------------

template<class CostFn = WeightedCost>
Path shortest_edge_to_edge(const FMM::NETWORK::Network& net,
                           EdgeIndex start_e,
                           EdgeIndex goal_e,
                           bool allow_u_turns = true,
                           bool include_start_edge_cost = false,
                           CostFn cost_fn = CostFn()) {
  typename LinkGraph<CostFn>::BuildOptions opt;
  opt.allow_u_turns = allow_u_turns;
  LinkGraph<CostFn> G(net, opt, cost_fn);

  double seed_cost = include_start_edge_cost ? cost_fn(net.get_edge(start_e)) : 0.0;

  auto R = dijkstra(G,
                    /*sources*/ {{start_e, seed_cost}},
                    /*isTarget*/ [goal_e](EdgeIndex e){ return e == goal_e; },
                    /*stopAt*/ goal_e);

  Path P = reconstruct_path(net, R.parent, goal_e, /*with_nodes=*/true);
  if (P.found) P.total_cost = R.dist[goal_e];
  return P;
}

template<class CostFn = WeightedCost>
std::unordered_map<EdgeIndex, Path>
shortest_edge_to_edges(const FMM::NETWORK::Network& net,
                       EdgeIndex start_e,
                       const std::vector<EdgeIndex>& goal_edges,
                       bool allow_u_turns = true,
                       bool include_start_edge_cost = false,
                       CostFn cost_fn = CostFn())
{
  SPDLOG_TRACE("Starting path search from edge {} to edges {}", start_e, goal_edges);
  std::unordered_map<EdgeIndex, Path> out;

  // --- Handle no goals ---
  if (goal_edges.empty()) {
    SPDLOG_TRACE("No goal edges, exiting early");
    return out;
  }
  const size_t E = static_cast<size_t>(net.get_edge_count());
  std::vector<char> is_goal(E, 0);
  size_t remaining = 0;

  // Mark unique goal edges
  for (auto g : goal_edges) {
    if (g < is_goal.size() && !is_goal[g]) {
      is_goal[g] = 1;
      ++remaining;
      out.emplace(g, Path{});
    }
  }

  // --- Trivial case: start edge is a target edge ---
  if (is_goal[start_e]) {
    Path P;
    P.found = true;
    P.total_cost = include_start_edge_cost ? cost_fn(net.get_edge(start_e)) : 0.0;

    // Empty edge list
    P.edges.clear();

    // Node path has at least source->target of the start edge
    const auto& es = net.get_edge(start_e);
    P.nodes.push_back(es.source);
    P.nodes.push_back(es.target);

    out[start_e] = std::move(P);
    is_goal[start_e] = 0;
    --remaining;

    // If that was the only goal, we are done
    if (remaining == 0) {
      SPDLOG_TRACE("No goals remaining, exiting early");
      return out;
    }
  }

  // ---- Build link graph ----
  typename LinkGraph<CostFn>::BuildOptions opt;
  opt.allow_u_turns = allow_u_turns;
  LinkGraph<CostFn> G(net, opt, cost_fn);

  // Seed: cost of entering the start edge
  double seed_cost = include_start_edge_cost ? cost_fn(net.get_edge(start_e)) : 0.0;
  std::vector<std::pair<EdgeIndex,double>> seeds;
  seeds.emplace_back(start_e, seed_cost);

  // --- Dijkstra state ---
  std::vector<double> dist(E, std::numeric_limits<double>::infinity());
  std::vector<EdgeIndex> parent(E, std::numeric_limits<EdgeIndex>::max());
  using QItem = std::pair<double, EdgeIndex>;
  std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;
  std::vector<char> settled(E, 0);

  // Seed PQ
  dist[start_e]   = seed_cost;
  parent[start_e] = std::numeric_limits<EdgeIndex>::max();
  pq.emplace(seed_cost, start_e);

  SPDLOG_TRACE("Starting main dijkstra loop {} goals to find", remaining);
  // ---- Multi-goal Dijkstra ----
  while (!pq.empty() && remaining > 0) {
    const QItem top = pq.top(); pq.pop();
    const double    du = top.first;
    const EdgeIndex u  = top.second;

    if (settled[u]) continue;
    settled[u] = 1;

    // If this edge is a goal, we have its optimal dist
    if (u < is_goal.size() && is_goal[u]) {
      Path& slot = out[u];
      if (!slot.found || du < slot.total_cost) {
        slot = reconstruct_path(net, parent, u, /*with_nodes=*/true);
        slot.total_cost = du;
        slot.found = true;
      }
      // Mark goal as satisfied
      is_goal[u] = 0;
      --remaining;
    }

    // Relax neighbors
    const auto& nbrs = G.neighbors(u);
    for (size_t i = 0; i < nbrs.size(); ++i) {
      const EdgeIndex v = nbrs[i].to;
      const double nd   = du + nbrs[i].w;

      if (nd + 1e-15 < dist[v]) {
        dist[v]   = nd;
        parent[v] = u;
        pq.emplace(nd, v);
      }
    }
  }

  return out;
}

} // namespace ROUTING
} // namespace FMM
