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

  PolyMatchResult match_traj(
      const CORE::Trajectory &traj,
      const POLYMATCHConfig &config,
      ROUTING::DijkstraState &state,
      ROUTING::IndexedMinHeap &heap,
      bool link_only_mode = false,
      MatchTimings *timings = nullptr);

private:
  const NETWORK::Network &network_;
  const NETWORK::PolygonLayer &polygon_layer_;
  const NETWORK::AccessPointLayer &ap_layer_;
  const ROUTING::PolyLinkGraph &poly_graph_;
  const ROUTING::LinkGraph &link_graph_;
};

} // MM
} // FMM

#endif
