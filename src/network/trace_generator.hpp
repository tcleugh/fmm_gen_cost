#ifndef FMM_NETWORK_TRACE_GENERATOR_HPP_
#define FMM_NETWORK_TRACE_GENERATOR_HPP_

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "core/geometry.hpp"
#include "network/access_point_layer.hpp"
#include "network/link_graph_routing.hpp"
#include "network/network.hpp"
#include "network/poly_link_graph.hpp"
#include "network/polygon_layer.hpp"
#include "network/trace_category.hpp"

namespace FMM {
namespace NETWORK {

// One trace ready for serialization. Transient — produced by TraceGenerator,
// immediately written to CSV.
struct GeneratedTrace {
  int id;
  FMM::CORE::LineString geom;
  TraceCategory category;
};

// Offline generator for polymatch's real-network validation suite.
// Constructed once over a real-area fixture set; each generate_* method
// produces ~n_per_category traces for one category. The categories the
// network topology can't support emit zero traces (e.g. no two polygons
// share an AP). See specs/002-real-network-validation/research.md R4 + R7.
class TraceGenerator {
 public:
  TraceGenerator(const Network& network, const PolygonLayer& polygons,
                 const AccessPointLayer& access_points,
                 const ROUTING::LinkGraph& link_graph,
                 const ROUTING::PolyLinkGraph& poly_graph,
                 uint64_t seed);

  // Each method returns up to n_per_category traces of its category, with
  // IDs in the bucket documented in contracts/real-network-trips-csv.md.
  std::vector<GeneratedTrace> generate_link_only(int n_per_category);
  std::vector<GeneratedTrace> generate_polygon_traversal(int n_per_category);
  std::vector<GeneratedTrace> generate_polygon_shared_ap(int n_per_category);
  std::vector<GeneratedTrace> generate_mid_polygon_start(int n_per_category);
  std::vector<GeneratedTrace> generate_mid_polygon_end(int n_per_category);
  std::vector<GeneratedTrace> generate_fully_inside(int n_per_category);
  std::vector<GeneratedTrace> generate_through_routing(int n_per_category);
  std::vector<GeneratedTrace> generate_off_network_noise(int n_per_category);
  std::vector<GeneratedTrace> generate_short_trip(int n_per_category);
  std::vector<GeneratedTrace> generate_duplicate_points(int n_per_category);

  // Writes the CSV in the layout documented in
  // specs/002-real-network-validation/contracts/real-network-trips-csv.md.
  // Sorts by id ascending; semicolon-delimited; LF line endings; classic
  // locale; setprecision(9). Returns true on success.
  bool write_csv(const std::string& out_path,
                 const std::vector<GeneratedTrace>& traces) const;

 private:
  // Helper random-walk used by link-only / off-network-noise / duplicate-points.
  // Picks a start edge whose midpoint has no polygon within search_radius,
  // then walks `hops` LinkGraph arcs and samples GPS points with Gaussian noise.
  // Returns empty geom if a suitable start edge can't be found.
  FMM::CORE::LineString random_walk_trace(int hops, int n_points,
                                          double noise_sigma,
                                          double search_radius);

  const Network& network_;
  const PolygonLayer& polygons_;
  const AccessPointLayer& access_points_;
  const ROUTING::LinkGraph& link_graph_;
  const ROUTING::PolyLinkGraph& poly_graph_;
  std::mt19937_64 rng_;
};

} // namespace NETWORK
} // namespace FMM

#endif // FMM_NETWORK_TRACE_GENERATOR_HPP_
