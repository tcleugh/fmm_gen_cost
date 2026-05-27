#include "mm/polymatch/polymatch_algorithm.hpp"

#include "algorithm/geom_algorithm.hpp"
#include "util/debug.hpp"
#include "util/util.hpp"
#include "mm/weightmatch/weightmatch_algorithm.hpp"

#include <limits>

namespace FMM {
namespace MM {

POLYMATCHConfig::POLYMATCHConfig(int k_arg, double r_arg, double gps_error_arg,
                                 int backup_k_arg, double backup_r_arg,
                                 double ub_factor_arg, bool allow_truncation_arg,
                                 double through_penalty_factor_arg,
                                 double boundary_epsilon_arg)
    : k(k_arg),
      radius(r_arg),
      gps_error(gps_error_arg),
      backup_k(backup_k_arg),
      backup_radius(backup_r_arg),
      upper_bound_factor(ub_factor_arg),
      allow_truncation(allow_truncation_arg),
      through_penalty_factor(through_penalty_factor_arg),
      boundary_epsilon(boundary_epsilon_arg) {}

void POLYMATCHConfig::print() const {
  SPDLOG_INFO("POLYMATCHAlgorithmConfig");
  SPDLOG_INFO(
      "k {} radius {} gps_error {} backup_k {} backup_radius {} upper_bound_factor {} allow_truncation {} through_penalty_factor {} boundary_epsilon {}",
      k, radius, gps_error, backup_k, backup_radius, upper_bound_factor,
      allow_truncation, through_penalty_factor, boundary_epsilon);
}

POLYMATCHConfig POLYMATCHConfig::load_from_arg(
    const cxxopts::ParseResult &arg_data) {
  int k = arg_data["candidates"].as<int>();
  double radius = arg_data["radius"].as<double>();
  double gps_error = arg_data["error"].as<double>();
  int backup_k = arg_data["backup_candidates"].as<int>();
  double backup_radius = arg_data["backup_radius"].as<double>();
  double ub_factor = arg_data["upper_bound_factor"].as<double>();
  bool allow_truncation = arg_data["allow_truncation"].as<bool>();
  double tpf = arg_data["through_penalty_factor"].as<double>();
  double eps = arg_data["boundary_epsilon"].as<double>();
  return POLYMATCHConfig{k, radius, gps_error, backup_k, backup_radius,
                          ub_factor, allow_truncation, tpf, eps};
}

void POLYMATCHConfig::register_arg(cxxopts::Options &options) {
  options.add_options()
    ("k,candidates", "Number of candidates",
      cxxopts::value<int>()->default_value("8"))
    ("r,radius", "Search radius",
      cxxopts::value<double>()->default_value("300.0"))
    ("e,error", "GPS error",
      cxxopts::value<double>()->default_value("50.0"))
    ("backup_candidates", "Number of candidates in fallback radius",
      cxxopts::value<int>()->default_value("-1"))
    ("backup_radius", "Fallback search radius",
      cxxopts::value<double>()->default_value("-1"))
    ("upper_bound_factor", "Dijkstra upper bound factor (0 to disable)",
      cxxopts::value<double>()->default_value("10.0"))
    ("allow_truncation", "Allow truncation of ends of trips in search",
      cxxopts::value<bool>()->default_value("false"))
    ("through_penalty_factor", "Through-routing cost multiplier",
      cxxopts::value<double>()->default_value("1.5"))
    ("boundary_epsilon", "Polygon boundary tolerance for AP validation",
      cxxopts::value<double>()->default_value("0.000001"));
}

void POLYMATCHConfig::register_help(std::ostringstream &oss) {
  oss << "-k/--candidates (optional) <int>: number of candidates (8)\n";
  oss << "-r/--radius (optional) <double>: search radius (300)\n";
  oss << "-e/--error (optional) <double>: GPS error (50)\n";
  oss << "--backup_candidates (optional) <int>: backup candidates (-1)\n";
  oss << "--backup_radius (optional) <double>: backup search radius (-1)\n";
  oss << "--upper_bound_factor (optional) <double>: Dijkstra upper bound factor (10.0)\n";
  oss << "--allow_truncation (optional) <bool>: allow truncation (false)\n";
  oss << "--through_penalty_factor (optional) <double>: through-routing multiplier (1.5)\n";
  oss << "--boundary_epsilon (optional) <double>: AP boundary tolerance (1e-6)\n";
}

bool POLYMATCHConfig::validate() const {
  if (gps_error <= 0 || radius <= 0 || k <= 0) {
    SPDLOG_CRITICAL("Invalid polymatch parameter k {} r {} gps_error {}", k,
                    radius, gps_error);
    return false;
  }
  if (through_penalty_factor < 0) {
    SPDLOG_CRITICAL("through_penalty_factor must be >= 0, got {}",
                    through_penalty_factor);
    return false;
  }
  if (backup_radius >= 0 && backup_radius <= radius) {
    SPDLOG_CRITICAL("backup_radius must exceed radius");
    return false;
  }
  if ((backup_radius >= 0 && backup_k < 0) ||
      (backup_radius < 0 && backup_k >= 0)) {
    SPDLOG_CRITICAL("backup_radius and backup_k must both be set or both unset");
    return false;
  }
  return true;
}

PolyMatchResult POLYMATCH::match_traj(
    const CORE::Trajectory &traj,
    const POLYMATCHConfig &config,
    ROUTING::DijkstraState &state,
    ROUTING::IndexedMinHeap &heap,
    bool link_only_mode,
    MatchTimings *timings) {
  // MVP fallback: when in link-only mode (or always for this skeleton),
  // run the WEIGHTMATCH algorithm and wrap its MatchResult.
  // Full hybrid link+polygon matching is implemented in subsequent tasks.
  WEIGHTMATCHConfig wm_config{config.k, config.radius, config.gps_error,
                              config.backup_k, config.backup_radius,
                              config.upper_bound_factor,
                              config.allow_truncation};
  WEIGHTMATCH wm(network_, link_graph_);
  PolyMatchResult result;
  result.base = wm.match_traj(traj, wm_config, state, heap, timings);
  return result;
}

} // MM
} // FMM
