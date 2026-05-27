#include "network/poly_link_graph.hpp"

#include <algorithm>
#include <boost/geometry.hpp>
#include <cassert>
#include <limits>
#include <unordered_set>

#include "util/debug.hpp"

namespace FMM {
namespace ROUTING {

namespace bg = boost::geometry;

PolyLinkGraph::PolyLinkGraph(const Network &network,
                             const LinkGraph &link_graph,
                             const PolygonLayer &polygon_layer,
                             const AccessPointLayer &ap_layer,
                             double through_penalty_factor)
    : n_edges_(static_cast<size_t>(network.get_edge_count())),
      n_polygons_(polygon_layer.size()),
      through_penalty_factor_(through_penalty_factor) {
  // Build sub-vertex bookkeeping (sizes + maps) before allocating adjacency.
  build_polygon_sub_vertices(polygon_layer, ap_layer);

  adjacency_.assign(n_edges_ + n_sub_vertices_, {});
  build_link_arcs(link_graph);
  build_polygon_arcs(network, polygon_layer, ap_layer);
}

void PolyLinkGraph::build_polygon_sub_vertices(
    const PolygonLayer &polygon_layer, const AccessPointLayer &ap_layer) {
  polygon_ap_to_subvertex_.assign(n_polygons_, {});
  polygon_local_ap_index_.assign(n_polygons_, {});
  sub_vertex_polygon_.clear();
  sub_vertex_ap_.clear();

  for (PolygonIndex p = 0; p < n_polygons_; ++p) {
    const auto &aps = ap_layer.aps_for_polygon(p);
    for (size_t k = 0; k < aps.size(); ++k) {
      AccessPointIndex ap = aps[k];
      uint32_t sv = static_cast<uint32_t>(n_edges_ + sub_vertex_polygon_.size());
      sub_vertex_polygon_.push_back(p);
      sub_vertex_ap_.push_back(ap);
      polygon_ap_to_subvertex_[p][ap] = sv;
      polygon_local_ap_index_[p][ap] = static_cast<uint16_t>(k);
    }
  }
  n_sub_vertices_ = sub_vertex_polygon_.size();

  // Pre-compute factor-independent through-cost table per polygon for the
  // introspection API. The factor-baked weight goes into arcs in
  // build_polygon_arcs; this table mirrors the original R11 design so a
  // caller can read weight * dist(a, b) without re-running the search.
  through_cost_tables_.assign(n_polygons_, {});
  for (PolygonIndex p = 0; p < n_polygons_; ++p) {
    const auto &aps = ap_layer.aps_for_polygon(p);
    size_t n = aps.size();
    through_cost_tables_[p].assign(n * n, 0.0);
    if (n == 0) continue;
    const auto &polys = polygon_layer.polygons();
    double w = polys[p].weight;
    for (size_t i = 0; i < n; ++i) {
      const auto &api = ap_layer.access_points()[aps[i]];
      for (size_t j = 0; j < n; ++j) {
        if (i == j) continue;
        const auto &apj = ap_layer.access_points()[aps[j]];
        through_cost_tables_[p][i * n + j] =
            w * bg::distance(api.point, apj.point);
      }
    }
  }
}

void PolyLinkGraph::build_link_arcs(const LinkGraph &link_graph) {
  for (size_t v = 0; v < n_edges_; ++v) {
    EdgeIndex e = static_cast<EdgeIndex>(v);
    const auto &nbrs = link_graph.neighbors(e);
    auto &out = adjacency_[v];
    out.reserve(nbrs.size());
    for (const auto &arc : nbrs) {
      out.push_back(Arc{static_cast<uint32_t>(arc.to), arc.w});
    }
  }
}

void PolyLinkGraph::build_polygon_arcs(const Network &network,
                                       const PolygonLayer &polygon_layer,
                                       const AccessPointLayer &ap_layer) {
  const auto &aps = ap_layer.access_points();
  const auto &edges = network.get_edges();

  // 1) link <-> polygon-sub-vertex arcs. We mirror LinkGraph's "arc weight
  //    is the cost-to-traverse-the-target" convention:
  //      edge E -> sub-vertex S: weight 0  (S is at a node, no cost to arrive
  //        having already paid for E earlier).
  //      sub-vertex S -> edge E: weight E.cost  (we must traverse E to leave
  //        the polygon onto E). Without this charge the polygon would "absorb"
  //        adjacent edges and break link-only routing in the matcher.
  //    For each AP with a network-node attachment, connect every incident
  //    edge to the sub-vertex that represents this AP within each polygon
  //    it references.
  for (const auto &ap : aps) {
    if (!ap.attached_node.has_value()) continue;
    for (PolygonIndex p_idx : ap.polygons) {
      auto it = polygon_ap_to_subvertex_[p_idx].find(ap.index);
      if (it == polygon_ap_to_subvertex_[p_idx].end()) continue;
      uint32_t sv = it->second;
      for (EdgeIndex e_idx : ap.attached_edges) {
        adjacency_[e_idx].push_back(Arc{sv, 0.0});
        adjacency_[sv].push_back(
            Arc{static_cast<uint32_t>(e_idx), edges[e_idx].cost});
      }
    }
  }

  // 2) Intra-polygon through-cost arcs sub-vertex (P, a) -> (P, b).
  //    The factor is baked in here so Dijkstra needs no runtime computation.
  for (PolygonIndex p = 0; p < n_polygons_; ++p) {
    const auto &local_aps = ap_layer.aps_for_polygon(p);
    size_t n = local_aps.size();
    if (n < 2) continue;
    double w_poly = polygon_layer.polygons()[p].weight;
    for (size_t i = 0; i < n; ++i) {
      uint32_t sv_i = polygon_ap_to_subvertex_[p][local_aps[i]];
      const auto &api = aps[local_aps[i]];
      for (size_t j = 0; j < n; ++j) {
        if (i == j) continue;
        uint32_t sv_j = polygon_ap_to_subvertex_[p][local_aps[j]];
        const auto &apj = aps[local_aps[j]];
        double w = w_poly * bg::distance(api.point, apj.point) *
                   through_penalty_factor_;
        adjacency_[sv_i].push_back(Arc{sv_j, w});
      }
    }
  }

  // 3) Cross-polygon shared-AP arcs: for every AP s shared between polygons
  //    P and Q, connect sub-vertex (P, s) <-> sub-vertex (Q, s) with weight 0.
  //    The same AP point in two different polygons is geometrically identical,
  //    so the transit between polygons is free (the polygons' own through-
  //    costs are paid by the intra-polygon arcs above).
  for (const auto &ap : aps) {
    if (ap.polygons.size() < 2) continue;
    for (size_t i = 0; i < ap.polygons.size(); ++i) {
      for (size_t j = 0; j < ap.polygons.size(); ++j) {
        if (i == j) continue;
        PolygonIndex pi = ap.polygons[i];
        PolygonIndex pj = ap.polygons[j];
        auto it_i = polygon_ap_to_subvertex_[pi].find(ap.index);
        auto it_j = polygon_ap_to_subvertex_[pj].find(ap.index);
        if (it_i == polygon_ap_to_subvertex_[pi].end() ||
            it_j == polygon_ap_to_subvertex_[pj].end()) {
          continue;
        }
        adjacency_[it_i->second].push_back(Arc{it_j->second, 0.0});
      }
    }
  }
}

uint32_t PolyLinkGraph::sub_vertex(PolygonIndex p, AccessPointIndex ap) const {
  if (p >= n_polygons_) return static_cast<uint32_t>(n_vertices());
  const auto &m = polygon_ap_to_subvertex_[p];
  auto it = m.find(ap);
  if (it == m.end()) return static_cast<uint32_t>(n_vertices());
  return it->second;
}

double PolyLinkGraph::through_cost_raw(PolygonIndex p, AccessPointIndex from,
                                       AccessPointIndex to) const {
  if (p >= n_polygons_) return std::numeric_limits<double>::infinity();
  const auto &map = polygon_local_ap_index_[p];
  auto it_f = map.find(from);
  auto it_t = map.find(to);
  if (it_f == map.end() || it_t == map.end()) {
    return std::numeric_limits<double>::infinity();
  }
  size_t n = map.size();
  return through_cost_tables_[p][it_f->second * n + it_t->second];
}

void shortest_polylink_to_polylinks(
    const PolyLinkGraph &G,
    DijkstraState &state,
    IndexedMinHeap &heap,
    uint32_t start_v,
    const std::vector<uint32_t> &goal_vs,
    std::vector<Path> &out,
    double upper_bound_factor) {
  const size_t V = G.n_vertices();
  const size_t K = goal_vs.size();
  state.ensure_size(V);
  state.next_epoch();
  heap.ensure_size(V);
  heap.clear();

  for (size_t i = 0; i < K; ++i) out[i] = Path{};
  if (K == 0) return;

  state.next_goal_epoch();
  size_t remaining = 0;
  size_t unique_goals = 0;
  for (size_t i = 0; i < K; ++i) {
    if (!state.is_goal(goal_vs[i])) {
      state.mark_goal(goal_vs[i]);
      ++remaining;
      ++unique_goals;
    }
  }

  const bool use_bound = upper_bound_factor > 0.0;
  const size_t min_found_for_bound =
      use_bound
          ? static_cast<size_t>(std::max(1, static_cast<int>(unique_goals) / 2))
          : 0;
  size_t goals_found = 0;
  double max_found_cost = 0.0;
  double bound = std::numeric_limits<double>::infinity();

  if (state.is_goal(start_v)) {
    for (size_t i = 0; i < K; ++i) {
      if (goal_vs[i] == start_v) {
        out[i].found = true;
        out[i].total_cost = 0.0;
      }
    }
    state.clear_goal(start_v);
    --remaining;
    if (remaining == 0) return;
  }

  state.set_dist(start_v, 0.0, kNoEdge);
  heap.push_or_decrease(start_v, 0.0);
  constexpr double eps = 1e-15;

  while (!heap.empty() && remaining > 0) {
    const uint32_t u = heap.pop_min();
    if (state.is_settled(u)) continue;
    state.mark_settled(u);
    const double du = state.get_dist(u);
    if (use_bound && du > bound) break;

    if (state.is_goal(u)) {
      Path p = reconstruct_path(state.parent, u);
      p.total_cost = du;
      p.found = true;
      for (size_t i = 0; i < K; ++i) {
        if (goal_vs[i] == u) out[i] = p;
      }
      state.clear_goal(u);
      --remaining;
      if (use_bound && du > 0.0) {
        ++goals_found;
        if (du > max_found_cost) max_found_cost = du;
        if (goals_found >= min_found_for_bound) {
          bound = max_found_cost * upper_bound_factor;
        }
      }
    }

    const auto &nbrs = G.neighbors(u);
    for (const auto &arc : nbrs) {
      const double nd = du + arc.w;
      if (nd + eps < state.get_dist(arc.to)) {
        state.set_dist(arc.to, nd, u);
        heap.push_or_decrease(arc.to, nd);
      }
    }
  }
}

} // ROUTING
} // FMM
