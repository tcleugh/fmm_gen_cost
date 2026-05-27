#include "mm/polymatch/polymatch_app.hpp"

using namespace FMM;
using namespace FMM::MM;

int main(int argc, char **argv) {
  POLYMATCHAppConfig config(argc, argv);
  if (config.help_specified) {
    POLYMATCHAppConfig::print_help();
    return 0;
  }
  if (!config.validate()) {
    return 0;
  }
  POLYMATCHApp app(config);
  app.run();
  return 0;
}
