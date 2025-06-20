/**
 * Fast map matching.
 *
 * Network configuration class
 *
 * @author Can Yang
 */

#ifndef FMM_PRIORITY_NETWORK_CONFIG_HPP_
#define FMM_PRIORITY_NETWORK_CONFIG_HPP_

#include <string>
#include "cxxopts/cxxopts.hpp"
// #include <boost/predef.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace FMM{

namespace CONFIG {
/**
 *  Network configuration class for reading network from a file.
 */
struct PriorityNetworkConfig{
  std::string file; /**< filename */
  std::string id; /**< id field/column name */
  std::string source; /**< source field/column name */
  std::string target; /**< target field/column name */
  std::string weight; /**< weight field/column name */

  /**
   * Check if the input is shapefile format
   */
  bool is_shapefile_format() const;
  /**
   * Validate the GPS configuration for file existence.
   * @return if file exists returns true, otherwise return false
   */
  bool validate() const;
    /**
   * Check if a priority network was provided.
   * @return if a priority network was provided retruns true, otherwise return false
   */
  bool use_priority_links() const;
  /**
   * Print informaiton
   */
  void print() const;
  /**
   * Load PriorityNetworkConfig from xml data
   * @param  xml_data XML data read from a file
   * @return PriorityNetworkConfig object containing information stored in the
   * xml file.
   */
  static PriorityNetworkConfig load_from_xml(
    const boost::property_tree::ptree &xml_data);
  /**
   * Load PriorityNetworkConfig from argument parsed data
   * @param  arg_data Argument data parsed from command line
   * @return PriorityNetworkConfig object containing information parsed from
   * command line argument.
   */
  static PriorityNetworkConfig load_from_arg(
    const cxxopts::ParseResult &arg_data);
  /**
   * Register arguments to an option object
   */
  static void register_arg(cxxopts::Options &options);
  /**
   * Register help information to a string stream
   */
  static void register_help(std::ostringstream &oss);
};

} // CONFIG

} // FMM

#endif //FMM_PRIORITY_NETWORK_CONFIG_HPP_
