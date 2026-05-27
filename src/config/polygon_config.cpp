#include "config/polygon_config.hpp"
#include "util/util.hpp"
#include "util/debug.hpp"

namespace FMM {
namespace CONFIG {

bool PolygonConfig::is_shapefile_format() const {
  return UTIL::check_file_extension(file, "shp");
}

void PolygonConfig::print() const {
  SPDLOG_INFO("PolygonConfig");
  SPDLOG_INFO("File name: {} ", file);
  SPDLOG_INFO("ID name: {} ", id_name);
  SPDLOG_INFO("Cost name: {} ", cost_name);
}

PolygonConfig PolygonConfig::load_from_xml(
  const boost::property_tree::ptree &xml_data) {
  std::string file = xml_data.get("config.input.polygon.file", "");
  std::string id_name = xml_data.get("config.input.polygon.id_name", "id");
  std::string cost_name = xml_data.get("config.input.polygon.cost_name", "cost");
  return PolygonConfig{file, id_name, cost_name};
}

PolygonConfig PolygonConfig::load_from_arg(
  const cxxopts::ParseResult &arg_data) {
  std::string file = arg_data["polygons"].as<std::string>();
  std::string id_name = arg_data["polygon_id_name"].as<std::string>();
  std::string cost_name = arg_data["polygon_cost_name"].as<std::string>();
  return PolygonConfig{file, id_name, cost_name};
}

void PolygonConfig::register_arg(cxxopts::Options &options) {
  options.add_options()
    ("polygons", "Polygon shapefile",
      cxxopts::value<std::string>()->default_value(""))
    ("polygon_id_name", "Polygon ID field name",
      cxxopts::value<std::string>()->default_value("id"))
    ("polygon_cost_name", "Polygon cost field name",
      cxxopts::value<std::string>()->default_value("cost"));
}

void PolygonConfig::register_help(std::ostringstream &oss) {
  oss << "--polygons (optional) <string>: polygon shapefile (empty: link-only fallback)\n";
  oss << "--polygon_id_name (optional) <string>: polygon ID field name (id)\n";
  oss << "--polygon_cost_name (optional) <string>: polygon cost field name (cost)\n";
}

bool PolygonConfig::validate() const {
  if (file.empty()) {
    return true;
  }
  if (!UTIL::file_exists(file)) {
    SPDLOG_CRITICAL("Polygon file not found {}", file);
    return false;
  }
  if (!is_shapefile_format()) {
    SPDLOG_CRITICAL("Polygon file format not recognized {}", file);
    return false;
  }
  return true;
}

} // CONFIG
} // FMM
