#ifndef FMM_NETWORK_POLY_LINK_GRAPH_HPP_
#define FMM_NETWORK_POLY_LINK_GRAPH_HPP_

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "network/type.hpp"
#include "network/network.hpp"
#include "network/link_graph_routing.hpp"
#include "network/polygon_layer.hpp"
#include "network/access_point_layer.hpp"

namespace FMM {
namespace ROUTING {

using FMM::NETWORK::Network;
using FMM::NETWORK::PolygonLayer;
using FMM::NETWORK::AccessPointLayer;
using FMM::NETWORK::AccessPointIndex;
using FMM::NETWORK::PolygonIndex;
using FMM::NETWORK::EdgeIndex;

enum class PolyVertexKind { Edge, PolygonSubVertex };

/**
 * Polygon-aware routing graph used by polymatch's Dijkstra core.
 *
 * Vertex layout (Held-Karp-style polygon-AP expansion per R11 alternative):
 *   [0, n_edges)                                  -> road EdgeIndex
 *   [n_edges, n_edges + n_sub_vertices)           -> polygon-AP sub-vertices
 *
 * Each polygon P with n_p access points contributes n_p sub-vertices, one per
 * AP. A sub-vertex represents "currently at AP a of polygon P" — so through-
 * routing cost lives in the static arc from sub-vertex (P, a) to sub-vertex
 * (P, b). Dijkstra needs no AP-context tracking: the entry AP is implicit in
 * the sub-vertex identity.
 *
 * Arc layout (all weights baked in at construction; no runtime computation):
 *   link E1 -> link E2:                arc weight = E2.weight (copied from LinkGraph)
 *   link E    <-> sub-vertex (P, a)    arc weight = 0  (a's incident edge or the edge whose endpoint coincides with a; entry/egress cost is observation-dependent and applied by the matcher at layer level)
 *   sub-vertex (P, a) -> (P, b):       arc weight = polygon.weight * dist(a, b) * THROUGH_PENALTY_FACTOR  (the through-routing cost)
 *   sub-vertex (P, s) <-> (Q, s):      arc weight = 0  (same AP s shared between P and Q)
 *
 * The polygon sub-vertex expansion means the polygon-internal through-cost is
 * fully encoded in arc weights, so a vanilla Dijkstra over this graph picks
 * polygon shortcuts whenever they beat the link-only route.
 *
 * Construction:    O(|E| + Sum |AP.polygons| + Sum |AP.attached_links| + Sum_p n_p^2)
 * Query:           O((|E| + n_sub) log (|E| + n_sub))  reusing IndexedMinHeap/DijkstraState
 */
class PolyLinkGraph {
public:
  struct Arc { uint32_t to; double w; };

  PolyLinkGraph(const Network &network,
                const LinkGraph &link_graph,
                const PolygonLayer &polygon_layer,
                const AccessPointLayer &ap_layer,
                double through_penalty_factor = 1.5);

  size_t n_edges() const { return n_edges_; }
  size_t n_polygons() const { return n_polygons_; }
  size_t n_sub_vertices() const { return n_sub_vertices_; }
  size_t n_vertices() const { return n_edges_ + n_sub_vertices_; }
  double through_penalty_factor() const { return through_penalty_factor_; }

  PolyVertexKind vertex_kind(uint32_t v) const {
    return v < n_edges_ ? PolyVertexKind::Edge
                        : PolyVertexKind::PolygonSubVertex;
  }

  // For sub-vertices only: which polygon does this sub-vertex belong to?
  PolygonIndex polygon_of(uint32_t v) const {
    return sub_vertex_polygon_[v - n_edges_];
  }
  // For sub-vertices only: which AccessPointIndex does this sub-vertex
  // represent within its polygon?
  AccessPointIndex ap_of(uint32_t v) const {
    return sub_vertex_ap_[v - n_edges_];
  }

  // Look up the sub-vertex ID for (polygon, access-point) pair. Returns
  // n_vertices() (out-of-range sentinel) if no such pair exists.
  uint32_t sub_vertex(PolygonIndex p, AccessPointIndex ap) const;

  const std::vector<Arc> &neighbors(uint32_t v) const { return adjacency_[v]; }

  // Factor-independent quantity polygon.weight * dist(AP_from, AP_to) — kept
  // as an introspection API even though the factor is now baked into the
  // arc weight at construction. Returns infinity if either AP isn't an AP of p.
  double through_cost_raw(PolygonIndex p, AccessPointIndex from,
                          AccessPointIndex to) const;

private:
  void build_link_arcs(const LinkGraph &link_graph);
  void build_polygon_sub_vertices(const PolygonLayer &polygon_layer,
                                  const AccessPointLayer &ap_layer);
  void build_polygon_arcs(const Network &network,
                          const PolygonLayer &polygon_layer,
                          const AccessPointLayer &ap_layer);

  size_t n_edges_ = 0;
  size_t n_polygons_ = 0;
  size_t n_sub_vertices_ = 0;
  double through_penalty_factor_ = 1.5;

  std::vector<std::vector<Arc>> adjacency_;

  // Sub-vertex bookkeeping.
  std::vector<PolygonIndex> sub_vertex_polygon_;     // indexed by sub-index
  std::vector<AccessPointIndex> sub_vertex_ap_;      // indexed by sub-index
  // (polygon p, ap_global) -> sub_vertex ID (in [n_edges, n_edges+n_sub))
  std::vector<std::unordered_map<AccessPointIndex, uint32_t>>
      polygon_ap_to_subvertex_;
  // Factor-independent raw cost table (for introspection / through_cost_raw).
  std::vector<std::vector<double>> through_cost_tables_;
  std::vector<std::unordered_map<AccessPointIndex, uint16_t>>
      polygon_local_ap_index_;
};

/**
 * Dijkstra over the PolyLinkGraph from `start_v` to all `goal_vs`.
 *
 * State + heap are sized to `n_vertices() = |E| + n_sub` and reused across
 * calls. `upper_bound_factor > 0` prunes nodes whose tentative distance
 * exceeds `max_found_cost * upper_bound_factor` after at least half the
 * unique goals are settled.
 *
 * Complexity: O((|E| + n_sub) log (|E| + n_sub)). All arc weights are static,
 * so per-arc relaxation is O(1) and the algorithm is a vanilla Dijkstra.
 */
void shortest_polylink_to_polylinks(
    const PolyLinkGraph &G,
    DijkstraState &state,
    IndexedMinHeap &heap,
    uint32_t start_v,
    const std::vector<uint32_t> &goal_vs,
    std::vector<Path> &out,
    double upper_bound_factor = 0.0);

} // ROUTING
} // FMM

#endif // FMM_NETWORK_POLY_LINK_GRAPH_HPP_
