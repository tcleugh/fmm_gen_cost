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

enum class PolyVertexKind { Edge, Polygon };

/**
 * Polygon-aware routing graph used by polymatch's Dijkstra core.
 *
 * Vertex ID space packs road-edge vertices [0, n_edges) and polygon vertices
 * [n_edges, n_edges + n_polygons) into a single contiguous range so the
 * existing IndexedMinHeap / DijkstraState (sized once to |E|+|P|) can be
 * reused per Constitution Principle I.
 *
 * Construction cost (one-shot at startup):
 *   O(|E|  + Sum |AP.polygons|  + Sum |AP.attached_links|  + Sum_p n_p^2)
 *   `--link arcs           `--link<->polygon arcs    `--through-cost table (R11)
 * The Sum_p n_p^2 term precomputes per-polygon `polygon.weight x dist(AP_i, AP_j)`
 * for every ordered pair so that during Dijkstra the through-routing cost is
 * one table lookup + one multiplication by THROUGH_PENALTY_FACTOR.
 *
 * Query cost (shortest_polylink_to_polylinks):
 *   O((|E|+|P|) log (|E|+|P|))
 * No per-query allocation: state and heap are reused across calls.
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
  size_t n_vertices() const { return n_edges_ + n_polygons_; }
  double through_penalty_factor() const { return through_penalty_factor_; }
  void set_through_penalty_factor(double f) { through_penalty_factor_ = f; }

  PolyVertexKind vertex_kind(uint32_t v) const {
    return v < n_edges_ ? PolyVertexKind::Edge : PolyVertexKind::Polygon;
  }

  const std::vector<Arc> &neighbors(uint32_t v) const { return adjacency_[v]; }

  double through_cost_raw(PolygonIndex p, AccessPointIndex from,
                          AccessPointIndex to) const;

private:
  void build_link_arcs(const LinkGraph &link_graph);
  void build_polygon_arcs(const AccessPointLayer &ap_layer);
  void build_through_cost_tables(const PolygonLayer &polygon_layer,
                                 const AccessPointLayer &ap_layer);

  size_t n_edges_ = 0;
  size_t n_polygons_ = 0;
  double through_penalty_factor_ = 1.5;

  std::vector<std::vector<Arc>> adjacency_;
  std::vector<std::vector<double>> through_cost_tables_;
  std::vector<std::unordered_map<AccessPointIndex, uint16_t>>
      polygon_local_ap_index_;
};

/**
 * Dijkstra over the PolyLinkGraph from `start_v` to all `goal_vs`.
 *
 * Reuses state + heap exactly like shortest_edge_to_edges, with state sized
 * to |E|+|P| (see PolyLinkGraph's vertex-ID-space invariant). When the
 * `upper_bound_factor` is positive, the search prunes nodes whose tentative
 * distance exceeds `max_found_cost * upper_bound_factor` once at least half
 * of the unique goals are settled.
 *
 * Complexity: O((|E|+|P|) log (|E|+|P|)). Per arc relaxation: O(1).
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
