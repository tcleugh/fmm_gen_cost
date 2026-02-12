
#include "mm/weightmatch/weightmatch_app.hpp"

using namespace FMM;
using namespace FMM::MM;

int main(int argc, char **argv){
  WEIGHTMATCHAppConfig config(argc,argv);
  if (config.help_specified) {
    WEIGHTMATCHAppConfig::print_help();
    return 0;
  }
  if (!config.validate()){
    return 0;
  }
  WEIGHTMATCHApp app(config);
  app.run();
  return 0;
};
