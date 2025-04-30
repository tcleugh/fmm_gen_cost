#include "config/network_config.hpp"
#include "util/util.hpp"
#include "util/debug.hpp"

void FMM::CONFIG::NetworkConfig::print() const{
  SPDLOG_INFO("NetworkConfig");
  SPDLOG_INFO("File name: {} ",file);
  SPDLOG_INFO("ID name: {} ",id);
  SPDLOG_INFO("Source name: {} ",source);
  SPDLOG_INFO("Target name: {} ",target);
  SPDLOG_INFO("Weight name: {} ",weight);
};

FMM::CONFIG::NetworkConfig FMM::CONFIG::NetworkConfig::load_from_xml(
  const boost::property_tree::ptree &xml_data){
  std::string file = xml_data.get<std::string>("config.input.network.file");
  std::string id = xml_data.get("config.input.network.id", "id");
  std::string source = xml_data.get("config.input.network.source","source");
  std::string target = xml_data.get("config.input.network.target","target");
  std::string weight = xml_data.get("config.input.network.weight","weight");
  return FMM::CONFIG::NetworkConfig{file, id, source, target, weight};
};

FMM::CONFIG::NetworkConfig FMM::CONFIG::NetworkConfig::load_from_arg(
  const cxxopts::ParseResult &arg_data){
  std::string file = arg_data["network"].as<std::string>();
  std::string id = arg_data["network_id"].as<std::string>();
  std::string source = arg_data["source"].as<std::string>();
  std::string target = arg_data["target"].as<std::string>();
  std::string weight = arg_data["weight"].as<std::string>();
  return FMM::CONFIG::NetworkConfig{file, id, source, target, weight};
};

void FMM::CONFIG::NetworkConfig::register_arg(cxxopts::Options &options){
  options.add_options()
  ("network","Network file name",
  cxxopts::value<std::string>()->default_value(""))
  ("network_id","Network id name",
  cxxopts::value<std::string>()->default_value("id"))
  ("source","Network source name",
  cxxopts::value<std::string>()->default_value("source"))
  ("target","Network target name",
  cxxopts::value<std::string>()->default_value("target"))
  ("weight","Network weight name",
    cxxopts::value<std::string>()->default_value("weight"));
};

void FMM::CONFIG::NetworkConfig::register_help(std::ostringstream &oss){
  oss<<"--network (required) <string>: Network file name\n";
  oss<<"--network_id (optional) <string>: Network id name (id)\n";
  oss<<"--source (optional) <string>: Network source name (source)\n";
  oss<<"--target (optional) <string>: Network target name (target)\n";
  oss<<"--weight (optional) <string>: Network weight name (weight)\n";
};

bool FMM::CONFIG::NetworkConfig::is_shapefile_format() const {
  if (FMM::UTIL::check_file_extension(file,"shp"))
    return true;
  return false;
};

bool FMM::CONFIG::NetworkConfig::validate() const {
  if (!UTIL::file_exists(file)){
    SPDLOG_CRITICAL("Network file not found {}",file);
    return false;
  }
  bool shapefile_format = is_shapefile_format();
  if (shapefile_format){
    return true;
  }
  SPDLOG_CRITICAL("Network format not recognized {}",file);
  return false;
}
