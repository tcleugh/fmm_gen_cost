// Offline trace generator for polymatch's real-network validation suite.
// Loads the real-area fixtures, runs each per-category constructor in
// TraceGenerator, sorts the combined batch by trace ID, and writes the CSV
// at the configured output path.
//
// See specs/002-real-network-validation/research.md R3 (determinism strategy)
// and contracts/real-network-trips-csv.md (CSV layout).

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "cxxopts/cxxopts.hpp"

#include "network/access_point_layer.hpp"
#include "network/link_graph_routing.hpp"
#include "network/network.hpp"
#include "network/poly_link_graph.hpp"
#include "network/polygon_layer.hpp"
#include "network/trace_generator.hpp"
#include "util/debug.hpp"

int main(int argc, char **argv) {
  spdlog::set_pattern("[%^%l%$][%s:%-3#] %v");
  spdlog::set_level(spdlog::level::info);

  cxxopts::Options options(
      "polymatch_traces_gen",
      "Generate the deterministic real-network trace CSV for polymatch.");
  options.add_options()
      ("network", "Network shapefile path",
       cxxopts::value<std::string>())
      ("polygons", "Polygon shapefile path",
       cxxopts::value<std::string>())
      ("access_points", "Access-point shapefile path",
       cxxopts::value<std::string>())
      ("output", "Output trips.csv path",
       cxxopts::value<std::string>())
      ("seed", "RNG seed (determinism gate)",
       cxxopts::value<uint64_t>()->default_value("2026"))
      ("n_per_category",
       "Traces emitted per category (≥ 20 required by SC-003)",
       cxxopts::value<int>()->default_value("20"))
      ("through_penalty_factor",
       "PolyLinkGraph through-routing factor at generation time",
       cxxopts::value<double>()->default_value("0.5"))
      ("h,help", "Print help and exit");

  if (argc == 1) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  std::string net_path, poly_path, ap_path, out_path;
  uint64_t seed = 2026;
  int n_per = 20;
  double tpf = 0.5;
  try {
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      return 0;
    }
    for (const char *required :
         {"network", "polygons", "access_points", "output"}) {
      if (!result.count(required)) {
        std::cerr << "missing required flag: --" << required << "\n"
                  << options.help() << std::endl;
        return 1;
      }
    }
    net_path = result["network"].as<std::string>();
    poly_path = result["polygons"].as<std::string>();
    ap_path = result["access_points"].as<std::string>();
    out_path = result["output"].as<std::string>();
    seed = result["seed"].as<uint64_t>();
    n_per = result["n_per_category"].as<int>();
    tpf = result["through_penalty_factor"].as<double>();
  } catch (const cxxopts::OptionException &e) {
    std::cerr << e.what() << "\n" << options.help() << std::endl;
    return 1;
  }

  SPDLOG_INFO("Loading network {}", net_path);
  FMM::NETWORK::Network net(net_path, "NO_TURN_BANS", "id", "source", "target",
                            "cost");
  FMM::ROUTING::LinkGraph link_graph(net);

  SPDLOG_INFO("Loading polygons {}", poly_path);
  FMM::NETWORK::PolygonLayer polygons;
  if (!polygons.load({poly_path, "id", "cost"})) {
    SPDLOG_CRITICAL("failed to load polygons");
    return 1;
  }

  SPDLOG_INFO("Loading access points {}", ap_path);
  FMM::NETWORK::AccessPointLayer aps;
  if (!aps.load({ap_path, "node_id", "polygon_id"}, polygons, net, 1e-6)) {
    SPDLOG_CRITICAL("failed to load access points");
    return 1;
  }

  FMM::ROUTING::PolyLinkGraph poly_graph(net, link_graph, polygons, aps, tpf);

  FMM::NETWORK::TraceGenerator gen(net, polygons, aps, link_graph, poly_graph,
                                   seed);

  std::vector<FMM::NETWORK::GeneratedTrace> all;
  auto append = [&](std::vector<FMM::NETWORK::GeneratedTrace> v) {
    for (auto &t : v) all.push_back(std::move(t));
  };
  SPDLOG_INFO("Generating link-only");
  append(gen.generate_link_only(n_per));
  SPDLOG_INFO("Generating polygon-traversal");
  append(gen.generate_polygon_traversal(n_per));
  SPDLOG_INFO("Generating polygon-shared-ap");
  append(gen.generate_polygon_shared_ap(n_per));
  SPDLOG_INFO("Generating mid-polygon-start");
  append(gen.generate_mid_polygon_start(n_per));
  SPDLOG_INFO("Generating mid-polygon-end");
  append(gen.generate_mid_polygon_end(n_per));
  SPDLOG_INFO("Generating fully-inside");
  append(gen.generate_fully_inside(n_per));
  SPDLOG_INFO("Generating through-routing");
  append(gen.generate_through_routing(n_per));
  SPDLOG_INFO("Generating off-network-noise");
  append(gen.generate_off_network_noise(n_per));
  SPDLOG_INFO("Generating short-trip");
  append(gen.generate_short_trip(n_per));
  SPDLOG_INFO("Generating duplicate-points");
  append(gen.generate_duplicate_points(n_per));

  SPDLOG_INFO("Writing CSV {} ({} traces)", out_path, all.size());
  if (!gen.write_csv(out_path, all)) {
    SPDLOG_CRITICAL("failed to write CSV");
    return 1;
  }
  std::cout << "Wrote " << all.size() << " traces to " << out_path
            << std::endl;
  return 0;
}
