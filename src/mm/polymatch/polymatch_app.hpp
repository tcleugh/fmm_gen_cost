#ifndef FMM_POLYMATCH_APP_H_
#define FMM_POLYMATCH_APP_H_

#include "network/link_graph_routing.hpp"
#include "network/polygon_layer.hpp"
#include "network/access_point_layer.hpp"
#include "network/poly_link_graph.hpp"
#include "mm/polymatch/polymatch_app_config.hpp"
#include "mm/polymatch/polymatch_algorithm.hpp"
#include "io/gps_reader.hpp"
#include "io/poly_mm_writer.hpp"

namespace FMM {
namespace MM {

class POLYMATCHApp {
public:
  POLYMATCHApp(const POLYMATCHAppConfig &config)
      : config_(config),
        network_(config_.network_config),
        link_graph_(network_) {};

  void run();

private:
  const POLYMATCHAppConfig &config_;
  NETWORK::Network network_;
  ROUTING::LinkGraph link_graph_;
};

} // MM
} // FMM

#endif
