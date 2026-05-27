#ifndef FMM_POLYMATCH_APP_CONFIG_HPP_
#define FMM_POLYMATCH_APP_CONFIG_HPP_

#include "config/gps_config.hpp"
#include "config/network_config.hpp"
#include "config/result_config.hpp"
#include "config/polygon_config.hpp"
#include "config/access_point_config.hpp"
#include "mm/polymatch/polymatch_algorithm.hpp"

namespace FMM {
namespace MM {

class POLYMATCHAppConfig {
public:
  POLYMATCHAppConfig(int argc, char **argv);

  void load_arg(int argc, char **argv);
  void load_xml(const std::string &file);
  static void print_help();
  void print() const;
  bool validate() const;

  CONFIG::NetworkConfig network_config;
  CONFIG::GPSConfig gps_config;
  CONFIG::ResultConfig result_config;
  CONFIG::PolygonConfig polygon_config;
  CONFIG::AccessPointConfig access_point_config;
  POLYMATCHConfig polymatch_config;
  bool use_omp = false;
  bool help_specified = false;
  int log_level = 2;
  int step = 100;
};

} // MM
} // FMM

#endif // FMM_POLYMATCH_APP_CONFIG_HPP_
