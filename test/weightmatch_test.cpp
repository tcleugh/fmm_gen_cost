#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "util/debug.hpp"
#include "network/network.hpp"
#include "network/link_graph_routing.hpp"
#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "mm/transition_graph.hpp"
#include "core/gps.hpp"
#include "io/gps_reader.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace FMM;
using namespace FMM::IO;
using namespace FMM::CORE;
using namespace FMM::NETWORK;
using namespace FMM::ROUTING;
using namespace FMM::MM;

struct GoldenRow {
  int id;
  std::string opath;
  std::string cpath;
};

static std::vector<GoldenRow> read_golden(const std::string &path) {
  std::vector<GoldenRow> rows;
  std::ifstream f(path);
  if (!f.is_open()) return rows;
  std::string line;
  std::getline(f, line); // header
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    GoldenRow row;
    std::stringstream ss(line);
    std::string token;
    std::getline(ss, token, ';'); row.id = std::stoi(token);
    std::getline(ss, token, ';'); row.opath = token;
    std::getline(ss, token, ';'); row.cpath = token;
    rows.push_back(row);
  }
  return rows;
}

template<typename T>
static std::string vec2str(const std::vector<T> &v) {
  std::string s;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i > 0) s += ",";
    s += std::to_string(v[i]);
  }
  return s;
}

TEST_CASE("weightmatch basic matching", "[weightmatch]") {
  spdlog::set_level(spdlog::level::warn);
  Network network("../../example/data/edges.shp", "NO_TURN_BANS",
    "id", "source", "target", "NO_WEIGHT");
  LinkGraph graph(network);

  WEIGHTMATCHConfig config(4, 0.5, 0.2);
  WEIGHTMATCH model(network, graph);

  CSVTrajectoryReader reader(
    "../../test/data/weightmatch/trips_basic.csv", "id", "geom");
  std::vector<Trajectory> trajectories = reader.read_all_trajectories();

  SECTION("basic_match against golden file") {
    auto golden = read_golden("../../test/data/weightmatch/expected_basic.csv");
    REQUIRE(!golden.empty());
    REQUIRE(trajectories.size() == golden.size());

    DijkstraState state;
    IndexedMinHeap heap;
    for (size_t i = 0; i < trajectories.size(); ++i) {
      MatchResult result = model.match_traj(
        trajectories[i], config, state, heap);
      CAPTURE(i, trajectories[i].id);
      CHECK(vec2str(result.opath) == golden[i].opath);
      CHECK(vec2str(result.cpath) == golden[i].cpath);
    }
  }

  SECTION("consistency - same result on repeated run") {
    DijkstraState state;
    IndexedMinHeap heap;
    // Run first 5 trajectories twice and compare
    size_t n = std::min<size_t>(5, trajectories.size());
    for (size_t i = 0; i < n; ++i) {
      MatchResult r1 = model.match_traj(
        trajectories[i], config, state, heap);
      MatchResult r2 = model.match_traj(
        trajectories[i], config, state, heap);
      CAPTURE(i);
      CHECK(r1.opath == r2.opath);
      CHECK(r1.cpath == r2.cpath);
    }
  }
}

TEST_CASE("weightmatch edge cases", "[weightmatch]") {
  spdlog::set_level(spdlog::level::warn);
  Network network("../../example/data/edges.shp", "NO_TURN_BANS",
    "id", "source", "target", "NO_WEIGHT");
  LinkGraph graph(network);

  WEIGHTMATCHConfig config(4, 0.5, 0.2);
  WEIGHTMATCH model(network, graph);

  CSVTrajectoryReader reader(
    "../../test/data/weightmatch/trips_edge_cases.csv", "id", "geom");
  std::vector<Trajectory> trajectories = reader.read_all_trajectories();

  SECTION("edge_cases against golden file") {
    auto golden = read_golden(
      "../../test/data/weightmatch/expected_edge_cases.csv");
    REQUIRE(!golden.empty());
    REQUIRE(trajectories.size() == golden.size());

    DijkstraState state;
    IndexedMinHeap heap;
    for (size_t i = 0; i < trajectories.size(); ++i) {
      MatchResult result = model.match_traj(
        trajectories[i], config, state, heap);
      CAPTURE(i, trajectories[i].id);
      CHECK(vec2str(result.opath) == golden[i].opath);
      CHECK(vec2str(result.cpath) == golden[i].cpath);
    }
  }

  SECTION("no crash on any edge case") {
    DijkstraState state;
    IndexedMinHeap heap;
    for (size_t i = 0; i < trajectories.size(); ++i) {
      REQUIRE_NOTHROW(model.match_traj(
        trajectories[i], config, state, heap));
    }
  }
}

// ---------------------------------------------------------------------------
// Real network tests: SA4=212 (small, 3826 edges)
// ---------------------------------------------------------------------------
static void run_golden_check(
    const std::string &network_path,
    const std::string &turn_ban_path,
    const std::string &trips_path,
    const std::string &golden_path,
    const WEIGHTMATCHConfig &config,
    const std::string &weight_name = "NO_WEIGHT") {
  Network network(network_path, turn_ban_path, "id", "source", "target", weight_name);
  LinkGraph graph(network);
  WEIGHTMATCH model(network, graph);

  CSVTrajectoryReader reader(trips_path, "id", "geom");
  std::vector<Trajectory> trajectories = reader.read_all_trajectories();
  auto golden = read_golden(golden_path);
  REQUIRE(!golden.empty());
  REQUIRE(trajectories.size() == golden.size());

  DijkstraState state;
  IndexedMinHeap heap;
  for (size_t i = 0; i < trajectories.size(); ++i) {
    MatchResult result = model.match_traj(
      trajectories[i], config, state, heap);
    CAPTURE(i, trajectories[i].id);
    CHECK(vec2str(result.opath) == golden[i].opath);
    CHECK(vec2str(result.cpath) == golden[i].cpath);
  }
}

TEST_CASE("weightmatch SA4=212 small real network", "[weightmatch][sa4_212]") {
  spdlog::set_level(spdlog::level::warn);
  WEIGHTMATCHConfig config(4, 300, 50);
  std::string base = "../../test/data/weightmatch/sa4_212/";

  SECTION("basic without turn bans") {
    run_golden_check(
      base + "links.shp", "NO_TURN_BANS",
      base + "trips_basic.csv", base + "expected_basic.csv", config, "weight");
  }

  SECTION("edge cases without turn bans") {
    run_golden_check(
      base + "links.shp", "NO_TURN_BANS",
      base + "trips_edge_cases.csv", base + "expected_edge_cases.csv", config, "weight");
  }

  SECTION("basic with turn bans") {
    run_golden_check(
      base + "links.shp", base + "turn_bans.csv",
      base + "trips_basic.csv", base + "expected_basic_tb.csv", config, "weight");
  }

  SECTION("consistency") {
    Network network(base + "links.shp", "NO_TURN_BANS", "id", "source", "target", "weight");
    LinkGraph graph(network);
    WEIGHTMATCH model(network, graph);

    CSVTrajectoryReader reader(base + "trips_basic.csv", "id", "geom");
    std::vector<Trajectory> trajectories = reader.read_all_trajectories();
    DijkstraState state;
    IndexedMinHeap heap;
    size_t n = std::min<size_t>(5, trajectories.size());
    for (size_t i = 0; i < n; ++i) {
      MatchResult r1 = model.match_traj(trajectories[i], config, state, heap);
      MatchResult r2 = model.match_traj(trajectories[i], config, state, heap);
      CAPTURE(i);
      CHECK(r1.opath == r2.opath);
      CHECK(r1.cpath == r2.cpath);
    }
  }
}

// ---------------------------------------------------------------------------
// Real network tests: SA4=17 (medium, 65270 edges)
// ---------------------------------------------------------------------------
TEST_CASE("weightmatch SA4=17 medium real network", "[weightmatch][sa4_17]") {
  spdlog::set_level(spdlog::level::warn);
  WEIGHTMATCHConfig config(4, 300, 50);
  std::string base = "../../test/data/weightmatch/sa4_17/";

  SECTION("basic without turn bans") {
    run_golden_check(
      base + "links.shp", "NO_TURN_BANS",
      base + "trips_basic.csv", base + "expected_basic.csv", config, "weight");
  }

  SECTION("edge cases without turn bans") {
    run_golden_check(
      base + "links.shp", "NO_TURN_BANS",
      base + "trips_edge_cases.csv", base + "expected_edge_cases.csv", config, "weight");
  }

  SECTION("basic with turn bans") {
    run_golden_check(
      base + "links.shp", base + "turn_bans.csv",
      base + "trips_basic.csv", base + "expected_basic_tb.csv", config, "weight");
  }

  SECTION("consistency") {
    Network network(base + "links.shp", "NO_TURN_BANS", "id", "source", "target", "weight");
    LinkGraph graph(network);
    WEIGHTMATCH model(network, graph);

    CSVTrajectoryReader reader(base + "trips_basic.csv", "id", "geom");
    std::vector<Trajectory> trajectories = reader.read_all_trajectories();
    DijkstraState state;
    IndexedMinHeap heap;
    size_t n = std::min<size_t>(5, trajectories.size());
    for (size_t i = 0; i < n; ++i) {
      MatchResult r1 = model.match_traj(trajectories[i], config, state, heap);
      MatchResult r2 = model.match_traj(trajectories[i], config, state, heap);
      CAPTURE(i);
      CHECK(r1.opath == r2.opath);
      CHECK(r1.cpath == r2.cpath);
    }
  }
}

// ---------------------------------------------------------------------------
// Tier 2 stress tests: many synthetic trips on real networks
// ---------------------------------------------------------------------------
TEST_CASE("weightmatch SA4=212 stress (1000 trips)", "[weightmatch][sa4_212][stress]") {
  spdlog::set_level(spdlog::level::warn);
  WEIGHTMATCHConfig config(4, 300, 50);
  std::string base = "../../test/data/weightmatch/sa4_212/";

  run_golden_check(
    base + "links.shp", "NO_TURN_BANS",
    base + "trips_stress.csv", base + "expected_stress.csv", config, "weight");
}

TEST_CASE("weightmatch SA4=17 stress (500 trips)", "[weightmatch][sa4_17][stress]") {
  spdlog::set_level(spdlog::level::warn);
  WEIGHTMATCHConfig config(4, 300, 50);
  std::string base = "../../test/data/weightmatch/sa4_17/";

  run_golden_check(
    base + "links.shp", "NO_TURN_BANS",
    base + "trips_stress.csv", base + "expected_stress.csv", config, "weight");
}
