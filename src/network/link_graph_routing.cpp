
#pragma once
#include <vector>
#include <queue>
#include <limits>
#include <functional>
#include <utility>
#include <cassert>
#include <cstdint>
#include <unordered_map>

#include "network/type.hpp"   // for Edge, EdgeIndex, NodeIndex, etc.
#include "network/network.hpp" // your Network class header (path as per your tree)

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

// ------------------------------------
// LinkGraph: vertices are EdgeIndex-es
// ------------------------------------
template<class CostFn = LengthCost>
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

template<class CostFn = LengthCost>
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

template<class CostFn = LengthCost>
Path shortest_edge_to_edge(const FMM::NETWORK::Network& net,
                           EdgeIndex start_e,
                           EdgeIndex goal_e,
                           bool allow_u_turns = true,
                           bool include_start_edge_cost = true,
                           CostFn cost_fn = CostFn()) {
  LinkGraph<CostFn> G(net, typename LinkGraph<CostFn>::BuildOptions{allow_u_turns}, cost_fn);

  double seed_cost = include_start_edge_cost ? cost_fn(net.get_edge(start_e)) : 0.0;

  auto R = dijkstra(G,
                    /*sources*/ {{start_e, seed_cost}},
                    /*isTarget*/ [goal_e](EdgeIndex e){ return e == goal_e; },
                    /*stopAt*/ goal_e);

  Path P = reconstruct_path(net, R.parent, goal_e, /*with_nodes=*/true);
  if (P.found) P.total_cost = R.dist[goal_e];
  return P;
}

template<class CostFn = LengthCost>
Path shortest_node_to_node(const FMM::NETWORK::Network& net,
                           NodeIndex start_n,
                           NodeIndex goal_n,
                           bool allow_u_turns = true,
                           CostFn cost_fn = CostFn()) {
  typename LinkGraph<CostFn>::BuildOptions opt;
  opt.allow_u_turns = allow_u_turns;                       
  LinkGraph<CostFn> G(net, opt, cost_fn);

  // seed with all edges leaving start node, cost = cost(out-edge)
  std::vector<std::pair<EdgeIndex,double>> seeds;
  for (const auto& e : net.get_edges()) {
    if (e.source == start_n) {
      seeds.emplace_back(e.index, cost_fn(e));
    }
  }
  if (seeds.empty()) return {};

  // target: any edge whose *target* is goal_n
  auto isTarget = [&net, goal_n](EdgeIndex e) {
    return net.get_edge(e).target == goal_n;
  };

  auto R = dijkstra(G, seeds, isTarget, /*stopAt*/ kNoEdge);

  // choose best among edges that end at goal
  double best = std::numeric_limits<double>::infinity();
  EdgeIndex best_e = std::numeric_limits<EdgeIndex>::max();
  for (size_t e = 0; e < R.dist.size(); ++e) {
    if (!std::isfinite(R.dist[e])) continue;
    if (net.get_edge(static_cast<EdgeIndex>(e)).target == goal_n && R.dist[e] < best) {
      best = R.dist[e];
      best_e = static_cast<EdgeIndex>(e);
    }
  }
  if (best_e == std::numeric_limits<EdgeIndex>::max()) return {};

  Path P = reconstruct_path(net, R.parent, best_e, /*with_nodes=*/true);
  P.total_cost = best;
  return P;
}


template<class CostFn = LengthCost>
std::unordered_map<NodeIndex, Path>
shortest_node_to_each_of(const FMM::NETWORK::Network& net,
                         NodeIndex start_n,
                         const std::vector<NodeIndex>& goal_nodes,
                         bool allow_u_turns = true,
                         CostFn cost_fn = CostFn())
{
  std::string line = "";
  for (auto g : goal_nodes) {
    line = line + std::to_string(g) + ",";
  }
  SPDLOG_TRACE("Starting path search from {} to {}", start_n, line);
  std::unordered_map<NodeIndex, Path> out;

  if (goal_nodes.empty()) {
    SPDLOG_TRACE("No goal nodes, exiting early");
    return out;
  }

  std::vector<char> is_goal(net.get_node_count(), 0);
  size_t remaining = 0;
  for (auto g : goal_nodes) {
    if (g < is_goal.size() && !is_goal[g]) {
      is_goal[g] = 1;
      ++remaining;
      out.emplace(g, Path{}); // create placeholders
    }
  }

//   // handle trivial direct matches
//   for (NodeIndex g : goal_nodes) {
//     if (g == start_n) {
//       Path P;
//       P.found = true;
//       P.total_cost = 0.0;
//       P.nodes = { start_n };
//       out[g] = P;
//       --remaining;
//       is_goal[g] = 0;
//     }
//   }
  if (remaining == 0) {
    SPDLOG_TRACE("No goals remaining, exiting early");
    return out;
  }

  typename LinkGraph<CostFn>::BuildOptions opt;
  opt.allow_u_turns = allow_u_turns;
  LinkGraph<CostFn> G(net, opt, cost_fn);

  // Seeds from start node
  std::vector<std::pair<EdgeIndex,double>> seeds;
  for (const auto& e : net.get_edges()) {
    if (e.source == start_n) seeds.emplace_back(e.index, cost_fn(e));
  }
  if (seeds.empty()) {
    SPDLOG_TRACE("No seed edges, exiting early");
    return out;
  }

  // Custom Dijkstra loop here so we can early‑stop after all goals are settled
  const size_t M = static_cast<size_t>(net.get_edge_count());
  std::vector<double> dist(M, std::numeric_limits<double>::infinity());
  std::vector<EdgeIndex> parent(M, std::numeric_limits<EdgeIndex>::max());
  using QItem = std::pair<double, EdgeIndex>;
  std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;
  std::vector<char> settled(M, 0);

  for (size_t i = 0; i < seeds.size(); ++i) {
    const EdgeIndex e = seeds[i].first;
    const double d0   = seeds[i].second;
    dist[e] = d0; parent[e] = std::numeric_limits<EdgeIndex>::max();
    pq.emplace(d0, e);
  }

  while (!pq.empty() && remaining > 0) {
    const QItem top = pq.top(); pq.pop();
    const double    du = top.first;
    const EdgeIndex u  = top.second;
    if (settled[u]) continue;
    settled[u] = 1;

    const auto& ue = net.get_edge(u);
    if (is_goal[ue.target]) {
      // produce a path for this goal if not already set (or if this is better)
      auto &slot = out[ue.target];
      if (!slot.found || du < slot.total_cost) {
        slot = reconstruct_path(net, parent, u, /*with_nodes=*/true);
        slot.total_cost = du;
        slot.found = true;
      }
      // We can mark this goal as satisfied once its first settling arrives (which is optimal).
      if (slot.total_cost == du) {
        is_goal[ue.target] = 0;
        --remaining;
      }
      // Continue; there might be other goals still pending.
    }

    for (const auto& arc : G.neighbors(u)) {
      EdgeIndex v = arc.to;
      double nd = du + arc.w;
      if (nd + 1e-15 < dist[v]) {
        dist[v] = nd; parent[v] = u;
        pq.emplace(nd, v);
      }
    }
  }

  // Any remaining goals stay with default (found=false)
  return out;
}

} // namespace ROUTING
} // namespace FMM
