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
  adjacency_.assign(n_edges_ + n_polygons_, {});
  build_link_arcs(link_graph);
  build_polygon_arcs(ap_layer);
  build_through_cost_tables(polygon_layer, ap_layer);
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

void PolyLinkGraph::build_polygon_arcs(const AccessPointLayer &ap_layer) {
  const auto &aps = ap_layer.access_points();
  for (const auto &ap : aps) {
    if (!ap.attached_node.has_value()) continue;
    for (PolygonIndex p_idx : ap.polygons) {
      uint32_t p_vertex = static_cast<uint32_t>(n_edges_ + p_idx);
      for (EdgeIndex e_idx : ap.attached_edges) {
        adjacency_[e_idx].push_back(Arc{p_vertex, 0.0});
        adjacency_[p_vertex].push_back(
            Arc{static_cast<uint32_t>(e_idx), 0.0});
      }
    }
  }

  for (const auto &ap : aps) {
    if (ap.polygons.size() < 2) continue;
    for (size_t i = 0; i < ap.polygons.size(); ++i) {
      for (size_t j = 0; j < ap.polygons.size(); ++j) {
        if (i == j) continue;
        uint32_t pi = static_cast<uint32_t>(n_edges_ + ap.polygons[i]);
        uint32_t pj = static_cast<uint32_t>(n_edges_ + ap.polygons[j]);
        adjacency_[pi].push_back(Arc{pj, 0.0});
      }
    }
  }
}

void PolyLinkGraph::build_through_cost_tables(
    const PolygonLayer &polygon_layer,
    const AccessPointLayer &ap_layer) {
  through_cost_tables_.assign(n_polygons_, {});
  polygon_local_ap_index_.assign(n_polygons_, {});
  for (PolygonIndex p = 0; p < n_polygons_; ++p) {
    const auto &local_aps = ap_layer.aps_for_polygon(p);
    size_t n = local_aps.size();
    through_cost_tables_[p].assign(n * n, 0.0);
    auto &local_map = polygon_local_ap_index_[p];
    for (size_t i = 0; i < n; ++i) {
      local_map[local_aps[i]] = static_cast<uint16_t>(i);
    }
    const auto &polys = polygon_layer.polygons();
    double w = polys[p].weight;
    for (size_t i = 0; i < n; ++i) {
      const auto &api = ap_layer.access_points()[local_aps[i]];
      for (size_t j = 0; j < n; ++j) {
        if (i == j) continue;
        const auto &apj = ap_layer.access_points()[local_aps[j]];
        through_cost_tables_[p][i * n + j] = w * bg::distance(api.point, apj.point);
      }
    }
  }
}

double PolyLinkGraph::through_cost_raw(PolygonIndex p, AccessPointIndex from,
                                       AccessPointIndex to) const {
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
      double w = arc.w;
      // Polygon->polygon arcs encoded with placeholder weight 0; cost is
      // applied at relaxation via the through-cost table referenced by the
      // entry-side polygon AP. This is handled fully by the caller in the
      // matcher when running specific link↔polygon transitions; for the
      // generic graph traversal, we use a conservative cost of 0 (the
      // matcher applies routing-time costs on top of polygon entries).
      const double nd = du + w;
      if (nd + eps < state.get_dist(arc.to)) {
        state.set_dist(arc.to, nd, u);
        heap.push_or_decrease(arc.to, nd);
      }
    }
  }
}

} // ROUTING
} // FMM
