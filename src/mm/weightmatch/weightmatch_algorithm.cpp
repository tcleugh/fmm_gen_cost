#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "algorithm/geom_algorithm.hpp"
#include "util/debug.hpp"
#include "util/util.hpp"
#include "io/gps_reader.hpp"
#include "io/mm_writer.hpp"
#include "network/link_graph_routing.cpp"

#include <limits>

using namespace FMM;
using namespace FMM::CORE;
using namespace FMM::NETWORK;
using namespace FMM::ROUTING;
using namespace FMM::MM;
using namespace FMM::PYTHON;

WEIGHTMATCHConfig::WEIGHTMATCHConfig(int k_arg, double r_arg, double gps_error_arg):
  k(k_arg), radius(r_arg), gps_error(gps_error_arg) {
};

void WEIGHTMATCHConfig::print() const {
  SPDLOG_INFO("WEIGHTMATCHAlgorithmConfig");
  SPDLOG_INFO("k {} radius {} gps_error {}", k, radius, gps_error);
};

WEIGHTMATCHConfig WEIGHTMATCHConfig::load_from_arg(
  const cxxopts::ParseResult &arg_data) {
  int k = arg_data["candidates"].as<int>();
  double radius = arg_data["radius"].as<double>();
  double gps_error = arg_data["error"].as<double>();
  return WEIGHTMATCHConfig{k, radius, gps_error};
};

void WEIGHTMATCHConfig::register_arg(cxxopts::Options &options){
  options.add_options()
    ("k,candidates","Number of candidates",
    cxxopts::value<int>()->default_value("8"))
    ("r,radius","Search radius",
    cxxopts::value<double>()->default_value("300.0"))
    ("e,error","GPS error",
    cxxopts::value<double>()->default_value("50.0"));
}

void WEIGHTMATCHConfig::register_help(std::ostringstream &oss){
  oss<<"-k/--candidates (optional) <int>: number of candidates (8)\n";
  oss<<"-r/--radius (optional) <double>: search "
    "radius (network data unit) (300)\n";
  oss<<"-e/--error (optional) <double>: GPS error "
    "(network data unit) (50)\n";
};

bool WEIGHTMATCHConfig::validate() const {
  if (gps_error <= 0 || radius <= 0 || k <= 0) {
    SPDLOG_CRITICAL("Invalid mm parameter k {} r {} gps error {} ", k, radius, gps_error);
    return false;
  }
  return true;
}

PyMatchResult WEIGHTMATCH::match_wkt(
  const std::string &wkt, const WEIGHTMATCHConfig &config) {
  LineString line = wkt2linestring(wkt);
  std::vector<double> timestamps;
  Trajectory traj{0, line, timestamps};
  MatchResult result = match_traj(traj, config);
  PyMatchResult output;
  output.id = result.id;
  output.opath = result.opath;
  output.cpath = result.cpath;
  output.mgeom = result.mgeom;
  output.indices = result.indices;
  for (int i = 0; i < result.opt_candidate_path.size(); ++i) {
    const MatchedCandidate &mc = result.opt_candidate_path[i];
    output.candidates.push_back(
      {i,
       mc.c.edge->id,
       graph_.get_node_id(mc.c.edge->source),
       graph_.get_node_id(mc.c.edge->target),
       mc.c.dist,
       mc.c.offset,
       mc.c.edge->length,
       mc.ep,
       mc.tp,
       mc.sp_dist}
      );
    output.pgeom.add_point(mc.c.point);
  }
  return output;
};

// Procedure of HMM based map matching algorithm.
MatchResult WEIGHTMATCH::match_traj(const Trajectory &traj, const WEIGHTMATCHConfig &config) {
  SPDLOG_DEBUG("Count of points in trajectory {}", traj.geom.get_num_points());
  SPDLOG_DEBUG("Search candidates");
  Traj_Candidates tc = network_.search_tr_cs_knn(traj.geom, config.k, config.radius);

  SPDLOG_DEBUG("Trajectory candidate {}", tc);
  if (tc.empty()) return MatchResult{};

  SPDLOG_DEBUG("Generate transition graph");
  TransitionGraph tg(tc, config.gps_error);

  SPDLOG_DEBUG("Update cost in transition graph");
  update_tg(&tg, traj, config);

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
  C_Path cpath = build_cpath(tg_opath, &indices);
  SPDLOG_DEBUG("Opath is {}", opath);
  SPDLOG_DEBUG("Indices is {}", indices);
  SPDLOG_DEBUG("Complete path is {}", cpath);

  LineString mgeom = network_.complete_path_to_geometry(traj.geom, cpath);

  return MatchResult{traj.id, matched_candidate_path, opath, cpath, indices, mgeom};
}

std::string WEIGHTMATCH::match_gps_file(
  const FMM::CONFIG::GPSConfig &gps_config,
  const FMM::CONFIG::ResultConfig &result_config,
  const WEIGHTMATCHConfig &weightmatch_config,
  bool use_omp
  ){
  std::ostringstream oss;
  std::string status;
  bool validate = true;
  if (!gps_config.validate()) {
    oss<<"gps_config invalid\n";
    validate = false;
  }
  if (!result_config.validate()) {
    oss<<"result_config invalid\n";
    validate = false;
  }
  if (!weightmatch_config.validate()) {
    oss<<"WEIGHTMATCH_config invalid\n";
    validate = false;
  }
  if (!validate) {
    oss<<"match_gps_file canceled\n";
    return oss.str();
  }
  // Start map matching
  int progress = 0;
  int points_matched = 0;
  int total_points = 0;
  int step_size = 1000;
  auto begin_time = UTIL::get_current_time();
  FMM::IO::GPSReader reader(gps_config);
  FMM::IO::CSVMatchResultWriter writer(result_config.file, result_config.output_config);
  if (use_omp) {
    int buffer_trajectories_size = 100000;
    while (reader.has_next_trajectory()) {
      std::vector<Trajectory> trajectories =
        reader.read_next_N_trajectories(buffer_trajectories_size);
      int trajectories_fetched = trajectories.size();
      #pragma omp parallel for
      for (int i = 0; i < trajectories_fetched; ++i) {
        Trajectory &trajectory = trajectories[i];
        int points_in_tr = trajectory.geom.get_num_points();
        MM::MatchResult result = match_traj(trajectory, weightmatch_config);
        writer.write_result(trajectory,result);
        #pragma omp critical
        if (!result.cpath.empty()) {
          points_matched += points_in_tr;
        }
        total_points += points_in_tr;
        ++progress;
        if (progress % step_size == 0) {
          std::stringstream buf;
          buf << "Progress " << progress << '\n';
          std::cout << buf.rdbuf();
        }
      }
    }
  } else {
    while (reader.has_next_trajectory()) {
      if (progress % step_size == 0) {
        SPDLOG_INFO("Progress {}", progress);
      }
      Trajectory trajectory = reader.read_next_trajectory();
      int points_in_tr = trajectory.geom.get_num_points();
      MM::MatchResult result = match_traj(trajectory, weightmatch_config);
      writer.write_result(trajectory,result);
      if (!result.cpath.empty()) {
        points_matched += points_in_tr;
      }
      total_points += points_in_tr;
      ++progress;
    }
  }
  auto end_time = UTIL::get_current_time();
  double duration = UTIL::get_duration(begin_time,end_time);
  oss<<"Status: success\n";
  oss<<"Time takes " << duration << " seconds\n";
  oss<<"Total points " << total_points << " matched "<< points_matched <<"\n";
  oss<<"Map match speed " << points_matched / duration << " points/s \n";
  return oss.str();
};

void WEIGHTMATCH::update_tg(TransitionGraph *tg, const Trajectory &traj, const WEIGHTMATCHConfig &config) {
  SPDLOG_DEBUG("Update transition graph");
  std::vector<TGLayer> &layers = tg->get_layers();
  std::vector<double> eu_dists = ALGORITHM::cal_eu_dist(traj.geom);
  int N = layers.size();
  for (int i = 0; i < N - 1; ++i) {
    // Routing from current_layer to next_layer
    update_layer(i, &(layers[i]), &(layers[i + 1]), eu_dists[i]);
  }
  SPDLOG_DEBUG("Update transition graph done");
}

void WEIGHTMATCH::update_layer(int level, TGLayer *la_ptr, TGLayer *lb_ptr, double eu_dist) {
  SPDLOG_DEBUG("Update layer {} starts", level);
  TGLayer &lb = *lb_ptr;
  for (auto iter_a = la_ptr->begin(); iter_a != la_ptr->end(); ++iter_a) {
    NodeIndex source = iter_a->c->edge->index;
    // SPDLOG_TRACE("  Calculate distance from source {}", source);
    // single source upper bound routing
    std::vector<NodeIndex> targets(lb.size());
    std::transform(lb.begin(), lb.end(), targets.begin(),
                   [](TGNode &a) {
      return a.c->edge->index;
    });

    std::unordered_map<EdgeIndex, Path> paths = shortest_edge_to_edges(
        network_,
        source,
        targets);
    for (auto& it: paths) {
        // Do stuff
        SPDLOG_TRACE("{} {}", it.first, it.second.total_cost);
        SPDLOG_TRACE("PATH: {}", it.second.edges);
    }

    for (auto iter_b = lb_ptr->begin(); iter_b != lb_ptr->end(); ++iter_b) {
      int i = std::distance(lb_ptr->begin(), iter_b);
      
      double path_distance;
      if (iter_a->c->edge->id == iter_b->c->edge->id) {
        path_distance = eu_dist;
      } else {
        double first_link_adjustment = iter_a->c->edge->weight * (iter_a->c->edge->length - iter_a->c->offset);
        double last_link_adjustment = iter_b->c->edge->weight * (iter_b->c->offset - iter_b->c->edge->length);
        path_distance = first_link_adjustment + paths.at(iter_b->c->edge->index).total_cost + last_link_adjustment;
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

C_Path WEIGHTMATCH::build_cpath(const TGOpath &opath, std::vector<int> *indices) {
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
    // SPDLOG_TRACE("Check a {} b {}", a->edge->id, b->edge->id);
    if (a->edge->id != b->edge->id) {
      Path path = shortest_edge_to_edge(
        network_,
        a->edge->index,
        b->edge->index);
      std::vector<EdgeIndex> segs = path.edges;
      // No transition found
      if (segs.empty() && a->edge->target != b->edge->source) {
        SPDLOG_TRACE("Candidate {} has disconnected edge {} to {}",
          i, a->edge->id, b->edge->id);
        indices->clear();
        return {};
      }
      for (int e:segs) {
        cpath.push_back(edges[e].id);
        ++current_idx;
      }
      cpath.push_back(b->edge->id);
      ++current_idx;
      indices->push_back(current_idx);
    } else {
      indices->push_back(current_idx);
    }
  }
  SPDLOG_DEBUG("Build cpath from optimal candidate path done");
  return cpath;
}
