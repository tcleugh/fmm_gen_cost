#include "mm/polymatch/polymatch_app_config.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "util/debug.hpp"
#include "util/util.hpp"

namespace FMM {
namespace MM {

POLYMATCHAppConfig::POLYMATCHAppConfig(int argc, char **argv) {
  spdlog::set_pattern("[%^%l%$][%s:%-3#] %v");
  // Mirror stmatch's auto-detect: a single XML argument means "load config
  // from XML file" (per `polymatch <config.xml>`).
  if (argc == 2) {
    std::string maybe(argv[1]);
    if (UTIL::check_file_extension(maybe, "xml,XML")) {
      load_xml(maybe);
    } else {
      load_arg(argc, argv);
    }
  } else {
    load_arg(argc, argv);
  }
  spdlog::set_level((spdlog::level::level_enum) log_level);
  if (!help_specified) print();
}

void POLYMATCHAppConfig::load_xml(const std::string &file) {
  SPDLOG_INFO("Start reading polymatch XML configuration {}", file);
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(file, tree);
  network_config = CONFIG::NetworkConfig::load_from_xml(tree);
  gps_config = CONFIG::GPSConfig::load_from_xml(tree);
  result_config = CONFIG::ResultConfig::load_from_xml(tree);
  polygon_config = CONFIG::PolygonConfig::load_from_xml(tree);
  access_point_config = CONFIG::AccessPointConfig::load_from_xml(tree);
  polymatch_config = POLYMATCHConfig::load_from_xml(tree);
  log_level = tree.get("config.other.log_level", 2);
  step = tree.get("config.other.step", 100);
  use_omp = !(!tree.get_child_optional("config.other.use_omp"));
  SPDLOG_INFO("Finish reading polymatch XML configuration");
}

void POLYMATCHAppConfig::load_arg(int argc, char **argv) {
  SPDLOG_INFO("Start reading polymatch configuration from arguments");
  cxxopts::Options options("polymatch_config", "Configuration parser");
  CONFIG::NetworkConfig::register_arg(options);
  CONFIG::GPSConfig::register_arg(options);
  CONFIG::ResultConfig::register_arg(options);
  CONFIG::PolygonConfig::register_arg(options);
  CONFIG::AccessPointConfig::register_arg(options);
  POLYMATCHConfig::register_arg(options);
  options.add_options()
    ("l,log_level", "Log level", cxxopts::value<int>()->default_value("2"))
    ("s,step", "Step report", cxxopts::value<int>()->default_value("100"))
    ("h,help", "Help information")
    ("use_omp", "Use omp or not");
  if (argc == 1) {
    help_specified = true;
    return;
  }
  auto result = options.parse(argc, argv);
  network_config = CONFIG::NetworkConfig::load_from_arg(result);
  gps_config = CONFIG::GPSConfig::load_from_arg(result);
  result_config = CONFIG::ResultConfig::load_from_arg(result);
  polygon_config = CONFIG::PolygonConfig::load_from_arg(result);
  access_point_config = CONFIG::AccessPointConfig::load_from_arg(result);
  polymatch_config = POLYMATCHConfig::load_from_arg(result);
  log_level = result["log_level"].as<int>();
  step = result["step"].as<int>();
  use_omp = result.count("use_omp") > 0;
  if (result.count("help") > 0) help_specified = true;
  SPDLOG_INFO("Finish reading polymatch arg configuration");
}

void POLYMATCHAppConfig::print() const {
  SPDLOG_INFO("---- Print configuration ----");
  network_config.print();
  gps_config.print();
  result_config.print();
  polygon_config.print();
  access_point_config.print();
  polymatch_config.print();
  SPDLOG_INFO("Log level {}", UTIL::LOG_LEVESLS[log_level]);
  SPDLOG_INFO("Step {}", step);
  SPDLOG_INFO("Use omp {}", (use_omp ? "true" : "false"));
  SPDLOG_INFO("---- Print configuration done ----");
}

void POLYMATCHAppConfig::print_help() {
  std::ostringstream oss;
  oss << "polymatch argument lists:\n";
  CONFIG::NetworkConfig::register_help(oss);
  CONFIG::GPSConfig::register_help(oss);
  POLYMATCHConfig::register_help(oss);
  CONFIG::ResultConfig::register_help(oss);
  CONFIG::PolygonConfig::register_help(oss);
  CONFIG::AccessPointConfig::register_help(oss);
  oss << "-l/--log_level (optional) <int>: log level (2)\n";
  oss << "-s/--step (optional) <int>: progress report step (100)\n";
  oss << "--use_omp: use OpenMP for multithreaded matching\n";
  oss << "-h/--help: print help information\n";
  std::cout << oss.str();
}

bool POLYMATCHAppConfig::validate() const {
  if (!gps_config.validate()) return false;
  if (!result_config.validate()) return false;
  if (!network_config.validate()) return false;
  if (!polygon_config.validate()) return false;
  if (!access_point_config.validate()) return false;
  if (!polymatch_config.validate()) return false;

  // Require both polygon + AP layers or neither
  bool poly_set = !polygon_config.file.empty();
  bool ap_set = !access_point_config.file.empty();
  if (poly_set != ap_set) {
    SPDLOG_CRITICAL("Polygon and access point layers must both be provided or both omitted");
    return false;
  }
  return true;
}

} // MM
} // FMM
