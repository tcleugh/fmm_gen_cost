#include "io/poly_mm_writer.hpp"

#include "util/util.hpp"
#include "util/debug.hpp"

#include <sstream>

namespace FMM {
namespace IO {

PolyMMWriter::PolyMMWriter(const std::string &filename,
                           const CONFIG::OutputConfig &output_config,
                           bool include_polygon_columns)
    : m_fstream(filename),
      config_(output_config),
      include_polygon_columns_(include_polygon_columns) {
  write_header();
}

void PolyMMWriter::write_header() {
  std::string header = "id";
  if (config_.write_opath) header += ";opath";
  if (config_.write_error) header += ";error";
  if (config_.write_offset) header += ";offset";
  if (config_.write_spdist) header += ";spdist";
  if (config_.write_pgeom) header += ";pgeom";
  if (config_.write_cpath) header += ";cpath";
  if (config_.write_tpath) header += ";tpath";
  if (config_.write_mgeom) header += ";mgeom";
  if (include_polygon_columns_) {
    header += ";polygon_ids;entry_aps;egress_aps;is_through;polygon_distances";
  }
  m_fstream << header << '\n';
}

void PolyMMWriter::write_result(const FMM::CORE::Trajectory &traj,
                                const FMM::MM::PolyMatchResult &result) {
  std::lock_guard<std::mutex> guard(write_mutex_);
  std::stringstream buf;
  buf << result.base.id;
  if (config_.write_opath) buf << ";" << result.base.opath;
  if (config_.write_error) {
    buf << ";";
    if (!result.base.opt_candidate_path.empty()) {
      int N = result.base.opt_candidate_path.size();
      for (int i = 0; i < N - 1; ++i) {
        buf << result.base.opt_candidate_path[i].c.dist << ",";
      }
      buf << result.base.opt_candidate_path[N - 1].c.dist;
    }
  }
  if (config_.write_offset) {
    buf << ";";
    if (!result.base.opt_candidate_path.empty()) {
      int N = result.base.opt_candidate_path.size();
      for (int i = 0; i < N - 1; ++i) {
        buf << result.base.opt_candidate_path[i].c.offset << ",";
      }
      buf << result.base.opt_candidate_path[N - 1].c.offset;
    }
  }
  if (config_.write_spdist) {
    buf << ";";
    if (!result.base.opt_candidate_path.empty()) {
      int N = result.base.opt_candidate_path.size();
      for (int i = 1; i < N; ++i) {
        buf << result.base.opt_candidate_path[i].sp_dist
            << (i == N - 1 ? "" : ",");
      }
    }
  }
  if (config_.write_pgeom) {
    buf << ";";
    if (!result.base.opt_candidate_path.empty()) {
      FMM::CORE::LineString pline;
      for (const auto &m : result.base.opt_candidate_path) {
        pline.add_point(m.c.point);
      }
      buf << pline;
    }
  }
  if (config_.write_cpath) buf << ";" << result.base.cpath;
  if (config_.write_tpath) {
    buf << ";";
    // tpath: reuse cpath-style writing via indices (per existing convention)
  }
  if (config_.write_mgeom) buf << ";" << result.base.mgeom;

  if (include_polygon_columns_) {
    // polygon_ids
    buf << ";";
    for (size_t i = 0; i < result.polygon_segments.size(); ++i) {
      if (i > 0) buf << ",";
      buf << result.polygon_segments[i].polygon_id;
    }
    // entry_aps
    buf << ";";
    for (size_t i = 0; i < result.polygon_segments.size(); ++i) {
      if (i > 0) buf << ",";
      if (result.polygon_segments[i].entry_ap == FMM::MM::kNoAccessPoint) {
        buf << "-";
      } else {
        buf << result.polygon_segments[i].entry_ap;
      }
    }
    // egress_aps
    buf << ";";
    for (size_t i = 0; i < result.polygon_segments.size(); ++i) {
      if (i > 0) buf << ",";
      if (result.polygon_segments[i].egress_ap == FMM::MM::kNoAccessPoint) {
        buf << "-";
      } else {
        buf << result.polygon_segments[i].egress_ap;
      }
    }
    // is_through
    buf << ";";
    for (size_t i = 0; i < result.polygon_segments.size(); ++i) {
      if (i > 0) buf << ",";
      buf << (result.polygon_segments[i].is_through ? 1 : 0);
    }
    // polygon_distances
    buf << ";";
    for (size_t i = 0; i < result.polygon_segments.size(); ++i) {
      if (i > 0) buf << ",";
      buf << result.polygon_segments[i].distance_inside;
    }
  }
  buf << '\n';
  m_fstream << buf.str();
}

} // IO
} // FMM
