#ifndef FMM_SRC_CONFIG_ACCESS_POINT_CONFIG_HPP_
#define FMM_SRC_CONFIG_ACCESS_POINT_CONFIG_HPP_

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cxxopts/cxxopts.hpp"

namespace FMM {
namespace CONFIG {

struct AccessPointConfig {
  AccessPointConfig(const std::string &file_arg = "",
                    const std::string &node_id_name_arg = "node_id",
                    const std::string &polygon_id_name_arg = "polygon_id") :
    file(file_arg),
    node_id_name(node_id_name_arg),
    polygon_id_name(polygon_id_name_arg) {};

  std::string file;
  std::string node_id_name;
  std::string polygon_id_name;

  bool is_shapefile_format() const;
  bool validate() const;
  void print() const;

  static AccessPointConfig load_from_xml(
    const boost::property_tree::ptree &xml_data);
  static AccessPointConfig load_from_arg(
    const cxxopts::ParseResult &arg_data);
  static void register_arg(cxxopts::Options &options);
  static void register_help(std::ostringstream &oss);
};

} // CONFIG
} // FMM

#endif // FMM_SRC_CONFIG_ACCESS_POINT_CONFIG_HPP_
