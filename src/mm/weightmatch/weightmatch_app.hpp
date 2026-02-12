
#ifndef FMM_WEIGHTMATCH_APP_H_
#define FMM_WEIGHTMATCH_APP_H_

#include "mm/weightmatch/weightmatch_app_config.hpp"
#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "io/gps_reader.hpp"
#include "io/mm_writer.hpp"

namespace FMM {
namespace MM{
/**
 * Class of weightmatch command line program
 */
class WEIGHTMATCHApp {
public:
  /**
   * Create weightmatch command application from configuration
   * @param config configuration of the command line app
   */
  WEIGHTMATCHApp(const WEIGHTMATCHAppConfig &config):
      config_(config),
      network_(config_.network_config),
      ng_(network_){};
  /**
   * Run the weightmatch program
   */
  void run();
 private:
  const WEIGHTMATCHAppConfig &config_;
  NETWORK::Network network_;
  NETWORK::NetworkGraph ng_;
};
}
}

#endif
