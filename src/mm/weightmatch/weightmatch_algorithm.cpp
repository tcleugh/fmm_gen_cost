#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "algorithm/geom_algorithm.hpp"
#include "util/debug.hpp"
#include "util/util.hpp"
#include "io/gps_reader.hpp"
#include "io/mm_writer.hpp"
#include "network/link_graph_routing.hpp"

#include <limits>

using namespace FMM;
using namespace FMM::CORE;
using namespace FMM::NETWORK;
using namespace FMM::ROUTING;
using namespace FMM::MM;
using namespace FMM::PYTHON;
using namespace FMM::ROUTING;

WEIGHTMATCHConfig::WEIGHTMATCHConfig(int k_arg, double r_arg, double gps_error_arg, int backup_k_arg, double backup_r_arg):
  k(k_arg), radius(r_arg), gps_error(gps_error_arg), backup_k(backup_k_arg), backup_radius(backup_r_arg) {
};

void WEIGHTMATCHConfig::print() const {
  SPDLOG_INFO("WEIGHTMATCHAlgorithmConfig");
  SPDLOG_INFO("k {} radius {} gps_error {} backup_k {} backup_radius {}", k, radius, gps_error, backup_k, backup_radius);
};

WEIGHTMATCHConfig WEIGHTMATCHConfig::load_from_arg(
  const cxxopts::ParseResult &arg_data) {
  int k = arg_data["candidates"].as<int>();
  double radius = arg_data["radius"].as<double>();
  double gps_error = arg_data["error"].as<double>();
  int backup_k = arg_data["backup_candidates"].as<int>();
  double backup_radius = arg_data["backup_radius"].as<double>();
  return WEIGHTMATCHConfig{k, radius, gps_error, backup_k, backup_radius};
};

void WEIGHTMATCHConfig::register_arg(cxxopts::Options &options){
  options.add_options()
    ("k,candidates","Number of candidates",
    cxxopts::value<int>()->default_value("8"))
    ("r,radius","Search radius",
    cxxopts::value<double>()->default_value("300.0"))
    ("e,error","GPS error",
    cxxopts::value<double>()->default_value("50.0"))
    ("backup_candidates","Number of candidates in fallback radius",
    cxxopts::value<int>()->default_value("-1"))
    ("backup_radius","Fallback search radius",
    cxxopts::value<double>()->default_value("-1"));
}

void WEIGHTMATCHConfig::register_help(std::ostringstream &oss){
  oss<<"-k/--candidates (optional) <int>: number of candidates (8)\n";
  oss<<"-r/--radius (optional) <double>: search "
    "radius (network data unit) (300)\n";
  oss<<"-e/--error (optional) <double>: GPS error "
    "(network data unit) (50)\n";
  oss<<"--backup_candidates (optional) <int>: number of candidates to use in backup search (-1)\n";
  oss<<"--backup_radius (optional) <double>: search "
    "radius to use if initial search fails (network data unit) (-1)\n";  
};

bool WEIGHTMATCHConfig::validate() const {
  if (gps_error <= 0 || radius <= 0 || k <= 0) {
    SPDLOG_CRITICAL("Invalid mm parameter k {} r {} gps error {} ", k, radius, gps_error);
    return false;
  }
  if (backup_radius >= 0 && backup_radius <= radius) {
    SPDLOG_CRITICAL("If backup radius is set, it must be larger than radius: r {} backup_r {} ", radius, backup_radius);
    return false;
  }
  if (
      (backup_radius >= 0 && backup_k < 0)
      || (backup_radius < 0 && backup_k >= 0)
   ) {
    SPDLOG_CRITICAL("If backup radius or candidates are set, then both must be set: backup_k {} backup_r {} ", backup_k, backup_radius);
    return false;
  }
  return true;
}

// Procedure of HMM based map matching algorithm.
MatchResult WEIGHTMATCH::match_traj(
  const Trajectory &traj, 
  const WEIGHTMATCHConfig &config, 
  DijkstraState& state, 
  IndexedMinHeap& heap
) {
  SPDLOG_DEBUG("Count of points in trajectory {}", traj.geom.get_num_points());
  SPDLOG_DEBUG("Search candidates");
  Traj_Candidates tc = network_.search_tr_cs_knn_with_fallback(
    traj.geom, 
    config.k, 
    config.radius, 
    config.backup_k, 
    config.backup_radius
  );

  SPDLOG_DEBUG("Trajectory candidate {}", tc);
  if (tc.empty()) return MatchResult{};

  SPDLOG_DEBUG("Generate transition graph");
  TransitionGraph tg(tc, config.gps_error);

  SPDLOG_DEBUG("Update cost in transition graph");
  update_tg(&tg, traj, config, state, heap);

  SPDLOG_DEBUG("Optimal path inference");
  TGOpath tg_opath = tg.backtrack();

  SPDLOG_DEBUG("Optimal path size {}", tg_opath.size());
  MatchedCandidatePath matched_candidate_path(tg_opath.size());

  std::transform(
    tg_opath.begin(), 
    tg_opath.end(),
    matched_candidate_path.begin(),
    [](const TGNode *a) {
      return MatchedCandidate{
        *(a->c), a->ep, a->tp, a->sp_dist
      };
    }
  );
  O_Path opath(tg_opath.size());
  std::transform(
    tg_opath.begin(), 
    tg_opath.end(),
    opath.begin(),
    [](const TGNode *a) {
      return a->c->edge->id;
    }
  );

  std::vector<int> indices;
  C_Path cpath = build_cpath(tg_opath, &indices, state, heap);
  SPDLOG_DEBUG("Opath is {}", opath);
  SPDLOG_DEBUG("Indices is {}", indices);
  SPDLOG_DEBUG("Complete path is {}", cpath);

  LineString mgeom = network_.complete_path_to_geometry(traj.geom, cpath);

  return MatchResult{traj.id, matched_candidate_path, opath, cpath, indices, mgeom};
}

void WEIGHTMATCH::update_tg(
  TransitionGraph *tg, 
  const Trajectory &traj, 
  const WEIGHTMATCHConfig &config, 
  DijkstraState& state, 
  IndexedMinHeap& heap
) {
  SPDLOG_DEBUG("Update transition graph");
  std::vector<TGLayer> &layers = tg->get_layers();
  std::vector<double> eu_dists = ALGORITHM::cal_eu_dist(traj.geom);
  int N = layers.size();
  for (int i = 0; i < N - 1; ++i) {
    // Routing from current_layer to next_layer
    update_layer(i, &(layers[i]), &(layers[i + 1]), eu_dists[i], state, heap);
  }
  SPDLOG_DEBUG("Update transition graph done");
}

void WEIGHTMATCH::update_layer(
  int level, 
  TGLayer *la_ptr, 
  TGLayer *lb_ptr, 
  double eu_dist, 
  DijkstraState& state, 
  IndexedMinHeap& heap
) {
  SPDLOG_DEBUG("Update layer {} starts", level);
  TGLayer &lb = *lb_ptr;
  for (auto iter_a = la_ptr->begin(); iter_a != la_ptr->end(); ++iter_a) {
    NodeIndex source = iter_a->c->edge->index;
    std::vector<NodeIndex> targets(lb.size());
    std::transform(lb.begin(), lb.end(), targets.begin(),
                   [](TGNode &a) {
      return a.c->edge->index;
    });

    std::unordered_map<EdgeIndex, Path> paths = shortest_edge_to_edges(
        graph_,
        state,
        heap,
        source,
        targets
    );

    for (auto iter_b = lb_ptr->begin(); iter_b != lb_ptr->end(); ++iter_b) {
      int i = std::distance(lb_ptr->begin(), iter_b);
      
      double path_distance;
      Edge* source_edge = iter_a->c->edge;
      Edge* target_edge = iter_b->c->edge;
      if (source_edge->id == target_edge->id) {
        path_distance = eu_dist;
      } else {
        path_distance = (
          source_edge->weight * (source_edge->length - iter_a->c->offset)// first_link_adjustment 
          + paths.at(target_edge->index).total_cost 
          + target_edge->weight * (iter_b->c->offset - target_edge->length) // last_link_adjustment
        );
      }

      double tp = TransitionGraph::calc_tp(path_distance, eu_dist);
      double temp = iter_a->cumu_prob + log(tp) + log(iter_b->ep);
      SPDLOG_TRACE("L {} f {} t {} sp {} dist {} tp {} ep {} fcp {} tcp {}",
        level, iter_a->c->edge->id,iter_b->c->edge->id,
        path_distance, eu_dist, tp, iter_b->ep, iter_a->cumu_prob,
        temp);
      if (temp >= iter_b->cumu_prob) {
        iter_b->cumu_prob = temp;
        iter_b->prev = &(*iter_a);
        iter_b->sp_dist = path_distance;
        iter_b->tp = tp;
      }
    }
  }
  SPDLOG_DEBUG("Update layer done");
}

C_Path WEIGHTMATCH::build_cpath(
  const TGOpath &opath, 
  std::vector<int> *indices, 
  DijkstraState& state, 
  IndexedMinHeap& heap
) {
  SPDLOG_DEBUG("Build cpath from optimal candidate path");
  C_Path cpath;
  if (!indices->empty()) indices->clear();
  if (opath.empty()) return cpath;
  const std::vector<Edge> &edges = network_.get_edges();
  int N = opath.size();
  cpath.push_back(opath[0]->c->edge->id);
  int current_idx = 0;
  // SPDLOG_TRACE("Insert index {}", current_idx);
  indices->push_back(current_idx);
  for (int i = 0; i < N - 1; ++i) {
    const Candidate *a = opath[i]->c;
    const Candidate *b = opath[i + 1]->c;
    
    Edge* source_edge = a->edge;
    Edge* target_edge = b->edge;
    if (source_edge->id != target_edge->id) {
      std::vector<EdgeIndex> target;
      target.push_back(b->edge->index);
      std::unordered_map<EdgeIndex, Path> paths = shortest_edge_to_edges(
        graph_,
        state,
        heap,
        source_edge->index,
        target
      );
      std::vector<EdgeIndex> segs = paths.at(target_edge->index).edges;

      // No transition found
      if (segs.empty() && source_edge->target != target_edge->source) {
        SPDLOG_TRACE(
          "Candidate {} has disconnected edge {} to {}",
          i, source_edge->id, target_edge->id
        );
        indices->clear();
        return {};
      }

      if (segs.size() > 2) {
        for (auto it = std::next(segs.begin());it != std::prev(segs.end()); ++it) {
          int e = *it;
          cpath.push_back(edges[e].id);
          ++current_idx;
        }
      }
      
      cpath.push_back(target_edge->id);
      ++current_idx;
      indices->push_back(current_idx);
    } else {
      indices->push_back(current_idx);
    }
  }
  SPDLOG_DEBUG("Build cpath from optimal candidate path done");
  return cpath;
}
