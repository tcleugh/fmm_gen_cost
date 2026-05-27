#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "util/debug.hpp"
#include "network/network.hpp"
#include "network/link_graph_routing.hpp"
#include "network/polygon_layer.hpp"
#include "network/access_point_layer.hpp"
#include "network/poly_link_graph.hpp"
#include "mm/polymatch/polymatch_algorithm.hpp"
#include "mm/polymatch/poly_match_result.hpp"
#include "config/polygon_config.hpp"
#include "config/access_point_config.hpp"

using namespace FMM;
using namespace FMM::CORE;
using namespace FMM::NETWORK;
using namespace FMM::ROUTING;
using namespace FMM::MM;

// Placeholder smoke test: verifies the polymatch headers link together.
// Full TDD test suite per tasks.md T020-T088 deferred to next implementation
// pass (requires the Python fixture generator in T003 to produce test data).
TEST_CASE("polymatch headers link and core types are accessible", "[polymatch][smoke]") {
  CONFIG::PolygonConfig poly_cfg;
  REQUIRE(poly_cfg.id_name == "id");
  REQUIRE(poly_cfg.cost_name == "cost");
  REQUIRE(poly_cfg.file.empty());
  REQUIRE(poly_cfg.validate());

  CONFIG::AccessPointConfig ap_cfg;
  REQUIRE(ap_cfg.node_id_name == "node_id");
  REQUIRE(ap_cfg.polygon_id_name == "polygon_id");
  REQUIRE(ap_cfg.validate());

  POLYMATCHConfig algo_cfg;
  REQUIRE(algo_cfg.k == 8);
  REQUIRE(algo_cfg.through_penalty_factor == Approx(1.5));
  REQUIRE(algo_cfg.boundary_epsilon == Approx(1e-6));
  REQUIRE(algo_cfg.validate());

  PolygonLayer layer;
  REQUIRE(layer.empty());
  REQUIRE(layer.load(poly_cfg));  // no file = success, empty
  REQUIRE(layer.empty());

  PolyMatchResult result;
  REQUIRE(result.polygon_segments.empty());
  REQUIRE(kNoAccessPoint == -1);
}
