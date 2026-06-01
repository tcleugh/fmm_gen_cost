#include "config/access_point_config.hpp"
#include "util/util.hpp"
#include "util/debug.hpp"

namespace FMM {
namespace CONFIG {

bool AccessPointConfig::is_shapefile_format() const {
  return UTIL::check_file_extension(file, "shp");
}

void AccessPointConfig::print() const {
  SPDLOG_INFO("AccessPointConfig");
  SPDLOG_INFO("File name: {} ", file);
  SPDLOG_INFO("Node ID name: {} ", node_id_name);
  SPDLOG_INFO("Polygon ID name: {} ", polygon_id_name);
}

AccessPointConfig AccessPointConfig::load_from_xml(
  const boost::property_tree::ptree &xml_data) {
  std::string file = xml_data.get("config.input.access_point.file", "");
  std::string node_id_name =
    xml_data.get("config.input.access_point.node_id_name", "node_id");
  std::string polygon_id_name =
    xml_data.get("config.input.access_point.polygon_id_name", "polygon_id");
  return AccessPointConfig{file, node_id_name, polygon_id_name};
}

AccessPointConfig AccessPointConfig::load_from_arg(
  const cxxopts::ParseResult &arg_data) {
  std::string file = arg_data["access_points"].as<std::string>();
  std::string node_id_name = arg_data["ap_node_id_name"].as<std::string>();
  std::string polygon_id_name =
    arg_data["ap_polygon_id_name"].as<std::string>();
  return AccessPointConfig{file, node_id_name, polygon_id_name};
}

void AccessPointConfig::register_arg(cxxopts::Options &options) {
  options.add_options()
    ("access_points", "Access point shapefile",
      cxxopts::value<std::string>()->default_value(""))
    ("ap_node_id_name", "Access point node ID field name",
      cxxopts::value<std::string>()->default_value("node_id"))
    ("ap_polygon_id_name", "Access point polygon ID field name",
      cxxopts::value<std::string>()->default_value("polygon_id"));
}

void AccessPointConfig::register_help(std::ostringstream &oss) {
  oss << "--access_points (optional) <string>: access point shapefile\n";
  oss << "--ap_node_id_name (optional) <string>: AP node ID field name (node_id)\n";
  oss << "--ap_polygon_id_name (optional) <string>: AP polygon ID field name (polygon_id)\n";
}

bool AccessPointConfig::validate() const {
  if (file.empty()) {
    return true;
  }
  if (!UTIL::file_exists(file)) {
    SPDLOG_CRITICAL("Access point file not found {}", file);
    return false;
  }
  if (!is_shapefile_format()) {
    SPDLOG_CRITICAL("Access point file format not recognized {}", file);
    return false;
  }
  return true;
}

} // CONFIG
} // FMM
