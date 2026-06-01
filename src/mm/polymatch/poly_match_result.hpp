#ifndef FMM_MM_POLYMATCH_POLY_MATCH_RESULT_HPP_
#define FMM_MM_POLYMATCH_POLY_MATCH_RESULT_HPP_

#include <cstdint>
#include <vector>

#include "core/geometry.hpp"
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

  // Geometry of the polygon traversal — populated by build_hybrid_path so
  // the matcher's mgeom builder can stitch together a continuous LineString
  // across edge and polygon segments. entry_point / egress_point are valid
  // when the corresponding *_ap field is not kNoAccessPoint. inside_points
  // are the matched-point coordinates of GPS observations that fell inside
  // this segment's polygon, in trajectory order.
  FMM::CORE::Point entry_point;
  FMM::CORE::Point egress_point;
  std::vector<FMM::CORE::Point> inside_points;
};

struct PolyMatchResult {
  MatchResult base;
  std::vector<PolygonSegment> polygon_segments;
};

} // MM
} // FMM

#endif // FMM_MM_POLYMATCH_POLY_MATCH_RESULT_HPP_
