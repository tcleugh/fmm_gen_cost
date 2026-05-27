#include "mm/polymatch/polymatch_app.hpp"

#include "io/gps_reader.hpp"
#include "io/poly_mm_writer.hpp"
#include "util/debug.hpp"
#include "util/util.hpp"

#include <algorithm>

namespace FMM {
namespace MM {

void POLYMATCHApp::run() {
  auto start_time = UTIL::get_current_time();

  NETWORK::PolygonLayer polygon_layer;
  NETWORK::AccessPointLayer ap_layer;
  bool link_only_mode = config_.polygon_config.file.empty();

  if (!link_only_mode) {
    if (!polygon_layer.load(config_.polygon_config)) {
      SPDLOG_CRITICAL("Failed to load polygon layer");
      return;
    }
    if (polygon_layer.empty()) {
      SPDLOG_WARN("Polygon layer is empty after load; falling back to link-only mode (FR-012)");
      link_only_mode = true;
    } else {
      if (!ap_layer.load(config_.access_point_config, polygon_layer, network_,
                         config_.polymatch_config.boundary_epsilon)) {
        SPDLOG_CRITICAL("Failed to load access point layer");
        return;
      }
    }
  } else {
    SPDLOG_INFO("Running polymatch in link-only mode (no polygon layer)");
  }

  ROUTING::PolyLinkGraph poly_graph(network_, link_graph_, polygon_layer,
                                    ap_layer,
                                    config_.polymatch_config.through_penalty_factor);
  POLYMATCH mm_model(network_, polygon_layer, ap_layer, poly_graph,
                     link_graph_);

  IO::GPSReader reader(config_.gps_config);
  IO::PolyMMWriter writer(config_.result_config.file,
                          config_.result_config.output_config,
                          !link_only_mode);

  int progress = 0;
  int points_matched = 0;
  int total_points = 0;
  int step_size = config_.step > 0 ? config_.step : 100;
  SPDLOG_INFO("Progress report step {}", step_size);

  if (config_.use_omp) {
    SPDLOG_INFO("Run polymatch parallelly");
    int buffer_size = 100000;
    while (reader.has_next_trajectory()) {
      auto trajectories = reader.read_next_N_trajectories(buffer_size);
      int n = trajectories.size();
      #pragma omp parallel
      {
        ROUTING::DijkstraState state;
        ROUTING::IndexedMinHeap heap;
        MatchTimings timings;
        #pragma omp for schedule(guided, 1)
        for (int i = 0; i < n; ++i) {
          CORE::Trajectory &traj = trajectories[i];
          int npts = traj.geom.get_num_points();
          if (npts < 2) {
            SPDLOG_WARN("Trajectory {} has only {} points; skipping (FR-019)",
                        traj.id, npts);
            continue;
          }
          PolyMatchResult result = mm_model.match_traj(
              traj, config_.polymatch_config, state, heap, link_only_mode,
              &timings);
          writer.write_result(traj, result);
          #pragma omp critical
          {
            if (!result.base.cpath.empty()) points_matched += npts;
            total_points += npts;
            ++progress;
          }
          if (progress % step_size == 0) {
            std::stringstream buf;
            buf << "Progress " << progress << '\n';
            std::cout << buf.rdbuf();
          }
        }
      }
    }
  } else {
    SPDLOG_INFO("Run polymatch in single thread");
    ROUTING::DijkstraState state;
    ROUTING::IndexedMinHeap heap;
    MatchTimings timings;
    while (reader.has_next_trajectory()) {
      if (progress % step_size == 0) {
        SPDLOG_INFO("Progress {}", progress);
      }
      CORE::Trajectory traj = reader.read_next_trajectory();
      int npts = traj.geom.get_num_points();
      if (npts < 2) {
        SPDLOG_WARN("Trajectory {} has only {} points; skipping (FR-019)",
                    traj.id, npts);
        ++progress;
        continue;
      }
      PolyMatchResult result = mm_model.match_traj(
          traj, config_.polymatch_config, state, heap, link_only_mode,
          &timings);
      writer.write_result(traj, result);
      if (!result.base.cpath.empty()) points_matched += npts;
      total_points += npts;
      ++progress;
    }
  }

  auto end_time = UTIL::get_current_time();
  double time_spent = UTIL::get_duration(start_time, end_time);
  SPDLOG_INFO("Polymatch finished; total trajectories {} matched points {} / {}",
              progress, points_matched, total_points);
  SPDLOG_INFO("Time {}", time_spent);
}

} // MM
} // FMM
