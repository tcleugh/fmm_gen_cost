
#ifndef FMM_WEIGHTMATCH_ALGORITHM_HPP
#define FMM_WEIGHTMATCH_ALGORITHM_HPP

#include "network/network.hpp"
#include "network/network_graph.hpp"
#include "mm/transition_graph.hpp"
#include "mm/mm_type.hpp"
#include "python/pyfmm.hpp"
#include "config/gps_config.hpp"
#include "config/result_config.hpp"
#include "network/link_graph_routing.hpp"

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cxxopts/cxxopts.hpp"

using namespace FMM::ROUTING;
using namespace FMM::CORE;

namespace FMM {
namespace MM {

struct MatchTimings {
  double candidate_search = 0;
  double update_tg = 0;
  double backtrack = 0;
  double build_cpath = 0;
  double geometry = 0;
};

/**
 * Configuration of stmatch algorithm
 */
struct WEIGHTMATCHConfig {
  /**
   * Constructor of stmatch algorithm configuration
   * @param k_arg the number of candidates
   * @param r_arg the search radius, in map unit, which is the same as
   * GPS data and network data.
   * @param gps_error_arg the gps error, in map unit
   * @param backup_k_arg the number of candidates to use for backup search
   * @param backup_r_arg the expanded search radius if no candidates are found in the inital search, 
   * in map unit, which is the same as GPS data and network data.
   */
  WEIGHTMATCHConfig(int k_arg = 8, double r_arg = 300, double gps_error_arg = 50, int backup_k_arg = -1, double backup_r_arg = -1, double ub_factor_arg = 10.0);
  int k; /**< number of candidates */
  double radius; /**< search radius for candidates, unit is map_unit*/
  double gps_error; /**< GPS error, unit is map_unit */
  int backup_k;
  double backup_radius;
  double upper_bound_factor; /**< Dijkstra upper bound: stop after max_found_cost * factor once enough goals found. 0 = disabled. */
  /**
   * Check the validity of the configuration
   */
  bool validate() const;
  /**
   * Print configuration data
   */
  void print() const;
  /**
   * Load from argument parsed data
   */
  static WEIGHTMATCHConfig load_from_arg(const cxxopts::ParseResult &arg_data);
  /**
   * Register arguments to an option object
   */
  static void register_arg(cxxopts::Options &options);
  /**
   * Register help information to a string stream
   */
  static void register_help(std::ostringstream &oss);
};

/**
 * %WEIGHTMATCH algorithm/model
 */
class WEIGHTMATCH {
public:
  /**
   * Create a stmatch model from network and graph
   */
  WEIGHTMATCH(const NETWORK::Network &network, const ROUTING::LinkGraph &graph) :
    network_(network), graph_(graph) {
  };
  /**
   * Match a trajectory to the road network
   * @param  traj   input trajector data
   * @param  config configuration of weightmatch algorithm
   * @return map matching result
   */
  MatchResult match_traj(
    const Trajectory &traj,
    const WEIGHTMATCHConfig &config,
    DijkstraState& state,
    IndexedMinHeap& heap,
    MatchTimings *timings = nullptr
  );
protected:
  /**
   * Update probabilities in a transition graph
   * @param tg transition graph
   * @param traj raw trajectory
   * @param config map match configuration
   */
  void update_tg(
    TransitionGraph *tg, 
    const CORE::Trajectory &traj, 
    const WEIGHTMATCHConfig &config,
    DijkstraState& state, 
    IndexedMinHeap& heap
  );
  /**
   * Update probabilities between two layers a and b in the transition graph
   * @param level   the index of layer a
   * @param la_ptr  layer a
   * @param lb_ptr  layer b next to a
   * @param eu_dist Euclidean distance between two observed point
   */
  void update_layer(
    int level,
    TGLayer *la_ptr,
    TGLayer *lb_ptr,
    double eu_dist,
    DijkstraState& state,
    IndexedMinHeap& heap,
    double upper_bound_factor
  );
  /**
   * Create a topologically connected path according to each matched
   * candidate
   * @param  tg_opath A sequence of optimal candidate nodes
   * @param  indices  the indices to be updated to store the index of matched
   * edge or candidate in the returned path.
   * @return A vector of edge id representing the traversed path
   */
  C_Path build_cpath(
    const TGOpath &tg_opath, 
    std::vector<int> *indices,     
    DijkstraState& state, 
    IndexedMinHeap& heap
  );
private:
  const NETWORK::Network& network_;
  const ROUTING::LinkGraph& graph_;
}; // WEIGHTMATCH
} // MM
} // FMM

#endif
