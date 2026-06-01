#ifndef FMM_POLYMATCH_ALGORITHM_HPP
#define FMM_POLYMATCH_ALGORITHM_HPP

#include "network/network.hpp"
#include "network/network_graph.hpp"
#include "network/polygon_layer.hpp"
#include "network/access_point_layer.hpp"
#include "network/poly_link_graph.hpp"
#include "mm/transition_graph.hpp"
#include "mm/mm_type.hpp"
#include "mm/polymatch/poly_match_result.hpp"
#include "mm/polymatch/poly_candidate.hpp"
#include "mm/polymatch/poly_transition_graph.hpp"
#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "config/gps_config.hpp"
#include "config/result_config.hpp"

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cxxopts/cxxopts.hpp"

namespace FMM {
namespace MM {

struct POLYMATCHConfig {
  POLYMATCHConfig(
      int k_arg = 8,
      double r_arg = 300,
      double gps_error_arg = 50,
      int backup_k_arg = -1,
      double backup_r_arg = -1,
      double ub_factor_arg = 10.0,
      bool allow_truncation_arg = false,
      double through_penalty_factor_arg = 1.5,
      double boundary_epsilon_arg = 1e-6);

  int k;
  double radius;
  double gps_error;
  int backup_k;
  double backup_radius;
  double upper_bound_factor;
  bool allow_truncation;
  double through_penalty_factor;
  double boundary_epsilon;

  bool validate() const;
  void print() const;
  static POLYMATCHConfig load_from_arg(const cxxopts::ParseResult &arg_data);
  static POLYMATCHConfig load_from_xml(
      const boost::property_tree::ptree &xml_data);
  static void register_arg(cxxopts::Options &options);
  static void register_help(std::ostringstream &oss);
};

class POLYMATCH {
public:
  POLYMATCH(const NETWORK::Network &network,
            const NETWORK::PolygonLayer &polygon_layer,
            const NETWORK::AccessPointLayer &ap_layer,
            const ROUTING::PolyLinkGraph &poly_graph,
            const ROUTING::LinkGraph &link_graph)
      : network_(network),
        polygon_layer_(polygon_layer),
        ap_layer_(ap_layer),
        poly_graph_(poly_graph),
        link_graph_(link_graph) {};

  // Match a trajectory. In link_only_mode (or when polygon layer is empty),
  // delegates to WEIGHTMATCH for binary-identical fallback per SC-002. Else
  // runs the polygon-aware matcher.
  PolyMatchResult match_traj(
      const CORE::Trajectory &traj,
      const POLYMATCHConfig &config,
      ROUTING::DijkstraState &state,
      ROUTING::IndexedMinHeap &heap,
      bool link_only_mode = false,
      MatchTimings *timings = nullptr);

private:
  // Phase A: candidate generation.
  PolyTrajCandidates build_candidates(
      const CORE::Trajectory &traj,
      const POLYMATCHConfig &config,
      // owned by caller; we keep raw pointers into it from PolyCandidate.edge
      MM::Traj_Candidates &link_candidates_owner) const;

  // Phase B: HMM update across consecutive layers (writes cumu_prob/prev into
  // each PolyTGNode in lb).
  void update_tg(PolyTransitionGraph &tg, const CORE::Trajectory &traj,
                 const POLYMATCHConfig &config,
                 ROUTING::DijkstraState &state,
                 ROUTING::IndexedMinHeap &heap);

  void update_layer(int level, PolyTGLayer *la, PolyTGLayer *lb, double eu_dist,
                    const POLYMATCHConfig &config, ROUTING::DijkstraState &state,
                    ROUTING::IndexedMinHeap &heap);

 public:
  // Compute the transition cost from candidate a -> b for given Euclidean gap.
  // Returns infinity if unreachable. Public for direct verification of the
  // cost model in tests (the four-case rule of FR-008, R5; same-polygon
  // eu_dist override).
  double transition_cost(const PolyCandidate &a, const PolyCandidate &b,
                         double eu_dist, const POLYMATCHConfig &config,
                         ROUTING::DijkstraState &state,
                         ROUTING::IndexedMinHeap &heap) const;

 private:
  // Min cost of routing between a single link edge and a polygon's matched
  // point via any of the polygon's access points. Handles both directions:
  // link_is_source=true routes link → polygon (single-source one-to-many
  // Dijkstra to the AP's attached edges); link_is_source=false routes
  // polygon → link (many-to-one — iterated as one Dijkstra per AP-attached
  // edge, since the routing layer has no native many-to-one primitive).
  // Used by the link→polygon and polygon→link branches of transition_cost.
  double link_polygon_routing_cost(
      FMM::NETWORK::Edge *link_edge,
      double link_endpoint_adjustment,
      FMM::NETWORK::PolygonIndex polygon_index,
      const FMM::CORE::Point &polygon_endpoint,
      bool link_is_source,
      const POLYMATCHConfig &config,
      ROUTING::DijkstraState &state,
      ROUTING::IndexedMinHeap &heap) const;

  // Phase C: hybrid C_Path + PolygonSegments from the optimal opath.
  void build_hybrid_path(const PolyTGOpath &opath,
                         const CORE::Trajectory &traj,
                         const POLYMATCHConfig &config,
                         ROUTING::DijkstraState &state,
                         ROUTING::IndexedMinHeap &heap,
                         PolyMatchResult *out) const;

  // Phase D: build the hybrid matched geometry by walking cpath, emitting an
  // edge's geometry for each positive entry and the polygon traversal
  // (entry AP -> inside GPS points -> egress AP) for each negated polygon ID.
  // No clipping of the first/last edge in v1 — that's a refinement.
  void build_hybrid_geometry(PolyMatchResult *out) const;

  const NETWORK::Network &network_;
  const NETWORK::PolygonLayer &polygon_layer_;
  const NETWORK::AccessPointLayer &ap_layer_;
  const ROUTING::PolyLinkGraph &poly_graph_;
  const ROUTING::LinkGraph &link_graph_;
};

} // MM
} // FMM

#endif
