
#include "mm/weightmatch/weightmatch_app.hpp"

#include <algorithm>

using namespace FMM;
using namespace FMM::CORE;
using namespace FMM::NETWORK;
using namespace FMM::MM;
using namespace FMM::ROUTING;

void WEIGHTMATCHApp::run() {
  auto start_time = UTIL::get_current_time();
  WEIGHTMATCH mm_model(network_, lg_);
  const WEIGHTMATCHConfig &weightmatch_config = config_.weightmatch_config;
  IO::GPSReader reader(config_.gps_config);
  IO::CSVMatchResultWriter writer(config_.result_config.file, config_.result_config.output_config);

  // Start map matching
  int progress = 0;
  int points_matched = 0;
  int total_points = 0;
  int step_size = 100;
  if (config_.step > 0) step_size = config_.step;
  SPDLOG_INFO("Progress report step {}", step_size);
  auto corrected_begin = UTIL::get_current_time();

  MatchTimings total_timings;
  size_t total_dijkstra_calls = 0;
  size_t total_nodes_explored = 0;
  std::vector<size_t> all_per_call_nodes;

  SPDLOG_INFO("Start to match trajectories");
  if (config_.use_omp){
    SPDLOG_INFO("Run map matching parallelly");

    int buffer_trajectories_size = 100000;
    while (reader.has_next_trajectory()) {
      std::vector<Trajectory> trajectories =
        reader.read_next_N_trajectories(buffer_trajectories_size);
      int trajectories_fetched = trajectories.size();

      #pragma omp parallel
      {
      DijkstraState state;
      IndexedMinHeap heap;
      MatchTimings thread_timings;

      #pragma omp for schedule(guided, 1)
      for (int i = 0; i < trajectories_fetched; ++i) {
        Trajectory &trajectory = trajectories[i];
        int points_in_tr = trajectory.geom.get_num_points();
        MM::MatchResult result = mm_model.match_traj(
            trajectory,
            weightmatch_config,
            state,
            heap,
            &thread_timings
          );
        writer.write_result(trajectory, result);

        #pragma omp critical
        {
        if (!result.cpath.empty()) {
          points_matched += points_in_tr;
        }
        total_points += points_in_tr;
        ++progress;
        }

        if (progress % step_size == 0) {
          std::stringstream buf;
          buf << "Progress " << progress << '\n';
          std::cout << buf.rdbuf();
        }
      }

      #pragma omp critical
      {
        total_timings.candidate_search += thread_timings.candidate_search;
        total_timings.update_tg += thread_timings.update_tg;
        total_timings.backtrack += thread_timings.backtrack;
        total_timings.build_cpath += thread_timings.build_cpath;
        total_timings.geometry += thread_timings.geometry;
        total_dijkstra_calls += state.dijkstra_calls;
        total_nodes_explored += state.nodes_explored;
        all_per_call_nodes.insert(all_per_call_nodes.end(),
          state.per_call_nodes.begin(), state.per_call_nodes.end());
      }
    }
  }
  } else {
    SPDLOG_INFO("Run map matching in single thread");
    DijkstraState state;
    IndexedMinHeap heap;
    while (reader.has_next_trajectory()) {
      if (progress % step_size == 0) {
        SPDLOG_INFO("Progress {}", progress);
      }

      Trajectory trajectory = reader.read_next_trajectory();
      int points_in_tr = trajectory.geom.get_num_points();

      MM::MatchResult result = mm_model.match_traj(
        trajectory,
        weightmatch_config,
        state,
        heap,
        &total_timings
      );
      writer.write_result(trajectory,result);

      if (!result.cpath.empty()) {
        points_matched += points_in_tr;
      }
      total_points += points_in_tr;
      ++progress;
    }
    total_dijkstra_calls = state.dijkstra_calls;
    total_nodes_explored = state.nodes_explored;
    all_per_call_nodes = std::move(state.per_call_nodes);
  }
  SPDLOG_INFO("MM process finished");
  auto end_time = UTIL::get_current_time();
  double time_spent = UTIL::get_duration(start_time, end_time);
  double time_spent_exclude_input = UTIL::get_duration(corrected_begin, end_time);
  SPDLOG_INFO("Time takes {}", time_spent);
  SPDLOG_INFO("Time takes excluding input {}", time_spent_exclude_input);
  SPDLOG_INFO("Finish map match total points {} matched {}", total_points, points_matched);
  SPDLOG_INFO("Matched percentage: {}", points_matched / (double) total_points);
  SPDLOG_INFO("Point match speed: {}", points_matched / time_spent);
  SPDLOG_INFO("Point match speed (excluding input): {}", points_matched / time_spent_exclude_input);
  SPDLOG_INFO("Time takes {}", time_spent);
  double timing_total = total_timings.candidate_search + total_timings.update_tg
    + total_timings.backtrack + total_timings.build_cpath + total_timings.geometry;
  SPDLOG_INFO("--- Phase timing breakdown (thread-seconds) ---");
  SPDLOG_INFO("  candidate_search: {:.3f}s ({:.1f}%)", total_timings.candidate_search,
    timing_total > 0 ? 100.0 * total_timings.candidate_search / timing_total : 0);
  SPDLOG_INFO("  update_tg:        {:.3f}s ({:.1f}%)", total_timings.update_tg,
    timing_total > 0 ? 100.0 * total_timings.update_tg / timing_total : 0);
  SPDLOG_INFO("  backtrack:        {:.3f}s ({:.1f}%)", total_timings.backtrack,
    timing_total > 0 ? 100.0 * total_timings.backtrack / timing_total : 0);
  SPDLOG_INFO("  build_cpath:      {:.3f}s ({:.1f}%)", total_timings.build_cpath,
    timing_total > 0 ? 100.0 * total_timings.build_cpath / timing_total : 0);
  SPDLOG_INFO("  geometry:         {:.3f}s ({:.1f}%)", total_timings.geometry,
    timing_total > 0 ? 100.0 * total_timings.geometry / timing_total : 0);
  SPDLOG_INFO("  total:            {:.3f}s", timing_total);
  SPDLOG_INFO("--- Dijkstra stats ---");
  SPDLOG_INFO("  total calls:      {}", total_dijkstra_calls);
  SPDLOG_INFO("  total nodes:      {}", total_nodes_explored);
  if (total_dijkstra_calls > 0) {
    SPDLOG_INFO("  mean nodes/call:  {:.1f}",
      static_cast<double>(total_nodes_explored) / total_dijkstra_calls);
  }
  if (!all_per_call_nodes.empty()) {
    std::sort(all_per_call_nodes.begin(), all_per_call_nodes.end());
    size_t n = all_per_call_nodes.size();
    auto pct = [&](double p) -> size_t {
      size_t idx = static_cast<size_t>(p * (n - 1));
      return all_per_call_nodes[idx];
    };
    SPDLOG_INFO("  p50 nodes/call:   {}", pct(0.50));
    SPDLOG_INFO("  p75 nodes/call:   {}", pct(0.75));
    SPDLOG_INFO("  p90 nodes/call:   {}", pct(0.90));
    SPDLOG_INFO("  p95 nodes/call:   {}", pct(0.95));
    SPDLOG_INFO("  p99 nodes/call:   {}", pct(0.99));
    SPDLOG_INFO("  max nodes/call:   {}", all_per_call_nodes.back());
  }
};
