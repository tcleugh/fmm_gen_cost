#ifndef FMM_MM_POLYMATCH_POLY_MATCH_RESULT_HPP_
#define FMM_MM_POLYMATCH_POLY_MATCH_RESULT_HPP_

#include <cstdint>
#include <vector>

#include "mm/mm_type.hpp"
#include "network/type.hpp"
#include "network/polygon_layer.hpp"

namespace FMM {
namespace MM {

constexpr FMM::NETWORK::NodeID kNoAccessPoint = -1;

struct PolygonSegment {
  FMM::NETWORK::PolygonID polygon_id;
  FMM::NETWORK::NodeID entry_ap;
  FMM::NETWORK::NodeID egress_ap;
  bool is_through;
  double distance_inside;
  size_t position_in_cpath;
};

struct PolyMatchResult {
  MatchResult base;
  std::vector<PolygonSegment> polygon_segments;
};

} // MM
} // FMM

#endif // FMM_MM_POLYMATCH_POLY_MATCH_RESULT_HPP_
