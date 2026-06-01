#ifndef FMM_POLY_MM_WRITER_HPP
#define FMM_POLY_MM_WRITER_HPP

#include "mm/polymatch/poly_match_result.hpp"
#include "config/result_config.hpp"
#include "core/gps.hpp"

#include <fstream>
#include <mutex>
#include <string>

namespace FMM {
namespace IO {

class PolyMMWriter {
public:
  PolyMMWriter(const std::string &filename,
               const CONFIG::OutputConfig &output_config,
               bool include_polygon_columns);

  void write_result(const FMM::CORE::Trajectory &traj,
                    const FMM::MM::PolyMatchResult &result);

private:
  void write_header();

  std::ofstream m_fstream;
  const CONFIG::OutputConfig config_;
  bool include_polygon_columns_;
  std::mutex write_mutex_;
};

} // IO
} // FMM

#endif
