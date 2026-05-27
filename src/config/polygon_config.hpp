#ifndef FMM_SRC_CONFIG_POLYGON_CONFIG_HPP_
#define FMM_SRC_CONFIG_POLYGON_CONFIG_HPP_

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cxxopts/cxxopts.hpp"

namespace FMM {
namespace CONFIG {

struct PolygonConfig {
  PolygonConfig(const std::string &file_arg = "",
                const std::string &id_name_arg = "id",
                const std::string &cost_name_arg = "cost") :
    file(file_arg), id_name(id_name_arg), cost_name(cost_name_arg) {};

  std::string file;
  std::string id_name;
  std::string cost_name;

  bool is_shapefile_format() const;
  bool validate() const;
  void print() const;

  static PolygonConfig load_from_xml(
    const boost::property_tree::ptree &xml_data);
  static PolygonConfig load_from_arg(
    const cxxopts::ParseResult &arg_data);
  static void register_arg(cxxopts::Options &options);
  static void register_help(std::ostringstream &oss);
};

} // CONFIG
} // FMM

#endif // FMM_SRC_CONFIG_POLYGON_CONFIG_HPP_
