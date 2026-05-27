#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "util/debug.hpp"
#include "network/network.hpp"
#include "network/link_graph_routing.hpp"
#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "mm/transition_graph.hpp"
#include "core/gps.hpp"
#include "core/geometry.hpp"
#include "io/gps_reader.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifndef FMM_TEST_DATA_DIR
#error "FMM_TEST_DATA_DIR must be set via target_compile_definitions"
#endif
#ifndef FMM_EXAMPLE_DATA_DIR
#error "FMM_EXAMPLE_DATA_DIR must be set via target_compile_definitions"
#endif

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
  Network network(std::string(FMM_EXAMPLE_DATA_DIR) + "/edges.shp", "NO_TURN_BANS",
    "id", "source", "target", "NO_WEIGHT");
  LinkGraph graph(network);

  WEIGHTMATCHConfig config(4, 0.5, 0.2);
  WEIGHTMATCH model(network, graph);

  CSVTrajectoryReader reader(
    std::string(FMM_TEST_DATA_DIR) + "/weightmatch/trips_basic.csv", "id", "geom");
  std::vector<Trajectory> trajectories = reader.read_all_trajectories();

  SECTION("basic_match against golden file") {
    auto golden = read_golden(std::string(FMM_TEST_DATA_DIR) + "/weightmatch/expected_basic.csv");
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
  Network network(std::string(FMM_EXAMPLE_DATA_DIR) + "/edges.shp", "NO_TURN_BANS",
    "id", "source", "target", "NO_WEIGHT");
  LinkGraph graph(network);

  WEIGHTMATCHConfig config(4, 0.5, 0.2);
  WEIGHTMATCH model(network, graph);

  CSVTrajectoryReader reader(
    std::string(FMM_TEST_DATA_DIR) + "/weightmatch/trips_edge_cases.csv", "id", "geom");
  std::vector<Trajectory> trajectories = reader.read_all_trajectories();

  SECTION("edge_cases against golden file") {
    auto golden = read_golden(
      std::string(FMM_TEST_DATA_DIR) + "/weightmatch/expected_edge_cases.csv");
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
  std::string base = std::string(FMM_TEST_DATA_DIR) + "/weightmatch/sa4_212/";

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
  std::string base = std::string(FMM_TEST_DATA_DIR) + "/weightmatch/sa4_17/";

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
  std::string base = std::string(FMM_TEST_DATA_DIR) + "/weightmatch/sa4_212/";

  run_golden_check(
    base + "links.shp", "NO_TURN_BANS",
    base + "trips_stress.csv", base + "expected_stress.csv", config, "weight");
}

TEST_CASE("weightmatch SA4=17 stress (500 trips)", "[weightmatch][sa4_17][stress]") {
  spdlog::set_level(spdlog::level::warn);
  WEIGHTMATCHConfig config(4, 300, 50);
  std::string base = std::string(FMM_TEST_DATA_DIR) + "/weightmatch/sa4_17/";

  run_golden_check(
    base + "links.shp", "NO_TURN_BANS",
    base + "trips_stress.csv", base + "expected_stress.csv", config, "weight");
}

// ---------------------------------------------------------------------------
// allow_truncation tests
//
// The example network has coordinates roughly in [0,5] x [0,5].
// Far points (100, 100) will have no candidates at any reasonable radius.
//
// "v" = matchable point (on the network), "_" = unmatchable (far away)
//
// GOOD — allow_truncation=true returns a valid match:
//   v v v v v   (all matchable — baseline)
//   _ _ v v v   (unmatchable head, contiguous matchable tail)
//   v v v _ _   (unmatchable tail, contiguous matchable head)
//   _ _ v v _   (unmatchable head and tail, contiguous middle)
//
// BAD — allow_truncation=true still returns empty:
//   _ _ _ _ _   (no matchable points)
//   v _ _ _ _   (only 1 matchable point — too few)
//   _ _ _ _ v   (only 1 matchable point — too few)
//   _ _ v _ _   (only 1 matchable point — too few)
//   v _ _ _ v   (2 matchable but non-contiguous — gap in middle)
//   _ v v _ v   (3 matchable but non-contiguous — gap in middle)
//   _ v _ v _   (2 matchable but non-contiguous — gap in middle)
//
// BAD — allow_truncation=false returns empty for any unmatchable point:
//   All cases above with an unmatchable point produce empty results.
// ---------------------------------------------------------------------------
TEST_CASE("weightmatch allow_truncation", "[weightmatch][truncation]") {
  spdlog::set_level(spdlog::level::warn);

  Network network(std::string(FMM_EXAMPLE_DATA_DIR) + "/edges.shp", "NO_TURN_BANS",
    "id", "source", "target", "NO_WEIGHT");
  LinkGraph graph(network);
  WEIGHTMATCH model(network, graph);

  // k=4, radius=0.5, error=0.2 matches the basic test config
  WEIGHTMATCHConfig config_no_trunc(4, 0.5, 0.2);
  WEIGHTMATCHConfig config_trunc(
    4, 0.5, 0.2,
    /*backup_k=*/-1, /*backup_radius=*/-1,
    /*ub_factor=*/10.0,
    /*allow_truncation=*/true
  );

  DijkstraState state;
  IndexedMinHeap heap;

  // Points from basic test trip 1 — known to produce a valid 7-point match.
  // Each coordinate is near the example network.
  const std::string p0 = "3.9892913938695225 0.9931575358298917";
  const std::string p1 = "2.8635963260768915 1.020206027395302";
  const std::string p2 = "1.9485599541199097 1.294953743205594";
  const std::string p3 = "1.4835299900169499 1.9240346487177404";
  const std::string p4 = "0.2855598394999991 2.0264517101104342";
  const std::string p5 = "0.8155601969034265 2.0992893810300375";
  const std::string p6 = "0.0016311686882696217 1.9705537989403834";
  const std::string far = "100.0 100.0";

  // Helper: check that opath has the expected pattern of -1s at head/tail
  // and valid (>=0) edge IDs in the middle.
  auto check_truncated_opath = [](
      const O_Path &opath,
      int expected_total,
      int unmatched_head,
      int unmatched_tail
  ) {
    CHECK((int)opath.size() == expected_total);
    for (int i = 0; i < unmatched_head; ++i)
      CHECK(opath[i] == -1);
    for (int i = unmatched_head; i < expected_total - unmatched_tail; ++i)
      CHECK(opath[i] >= 0);
    for (int i = expected_total - unmatched_tail; i < expected_total; ++i)
      CHECK(opath[i] == -1);
  };

  // --- GOOD cases: allow_truncation=true produces a valid match ---

  SECTION("good: all matchable (baseline sanity)") {
    // v v v v v v v — no truncation needed; both configs give the same non-empty result
    Trajectory inner{1, wkt2linestring(
      "LINESTRING(" + p0 + "," + p1 + "," + p2 + "," + p3 + "," + p4 + "," + p5 + "," + p6 + ")")};
    MatchResult r_no  = model.match_traj(inner, config_no_trunc, state, heap);
    MatchResult r_yes = model.match_traj(inner, config_trunc,    state, heap);
    REQUIRE(!r_no.opath.empty());
    CHECK(r_no.opath  == r_yes.opath);
    CHECK(r_no.cpath  == r_yes.cpath);
  }

  SECTION("good: unmatchable head (_ _ v v v v v v v)") {
    // Two far points prepended; inner 7 points are matchable
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," + far + "," +
      p0 + "," + p1 + "," + p2 + "," + p3 + "," + p4 + "," + p5 + "," + p6 + ")")};
    MatchResult with_trunc  = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc    = model.match_traj(t, config_no_trunc, state, heap);
    // With truncation: opath size = 9, first 2 entries = -1, rest are valid edge IDs
    check_truncated_opath(with_trunc.opath, /*total=*/9, /*head=*/2, /*tail=*/0);
    CHECK(!with_trunc.cpath.empty());
    // Without truncation: fails because of the unmatchable head points
    CHECK(no_trunc.opath.empty());
    CHECK(no_trunc.cpath.empty());
  }

  SECTION("good: unmatchable tail (v v v v v v v _ _)") {
    // Inner 7 points followed by two far points
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + p0 + "," + p1 + "," + p2 + "," + p3 + "," + p4 + "," + p5 + "," + p6 + "," +
      far + "," + far + ")")};
    MatchResult with_trunc  = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc    = model.match_traj(t, config_no_trunc, state, heap);
    check_truncated_opath(with_trunc.opath, /*total=*/9, /*head=*/0, /*tail=*/2);
    CHECK(!with_trunc.cpath.empty());
    CHECK(no_trunc.opath.empty());
    CHECK(no_trunc.cpath.empty());
  }

  SECTION("good: unmatchable head and tail (_ v v v v v v v _)") {
    // One far point at each end, inner 7 matchable
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," +
      p0 + "," + p1 + "," + p2 + "," + p3 + "," + p4 + "," + p5 + "," + p6 + "," +
      far + ")")};
    MatchResult with_trunc  = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc    = model.match_traj(t, config_no_trunc, state, heap);
    check_truncated_opath(with_trunc.opath, /*total=*/9, /*head=*/1, /*tail=*/1);
    CHECK(!with_trunc.cpath.empty());
    CHECK(no_trunc.opath.empty());
    CHECK(no_trunc.cpath.empty());
  }

  // --- BAD cases: allow_truncation=true still returns empty ---

  SECTION("bad: all unmatchable (_ _ _ _ _)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," + far + "," + far + "," + far + "," + far + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }

  SECTION("bad: only 1 matchable point at head (v _ _ _ _)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + p0 + "," + far + "," + far + "," + far + "," + far + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }

  SECTION("bad: only 1 matchable point at tail (_ _ _ _ v)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," + far + "," + far + "," + far + "," + p6 + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }

  SECTION("bad: only 1 matchable point in middle (_ _ v _ _)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," + far + "," + p3 + "," + far + "," + far + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }

  SECTION("bad: gap in middle — 2 matchable non-contiguous (v _ _ _ v)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + p0 + "," + far + "," + far + "," + far + "," + p6 + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }

  SECTION("bad: gap in middle — 3 matchable non-contiguous (_ v v _ v)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," + p0 + "," + p1 + "," + far + "," + p6 + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }

  SECTION("bad: gap in middle — alternating (_ v _ v _)") {
    Trajectory t{1, wkt2linestring(
      "LINESTRING(" + far + "," + p1 + "," + far + "," + p4 + "," + far + ")")};
    MatchResult with_trunc = model.match_traj(t, config_trunc,    state, heap);
    MatchResult no_trunc   = model.match_traj(t, config_no_trunc, state, heap);
    CHECK(with_trunc.opath.empty());
    CHECK(no_trunc.opath.empty());
  }
}
