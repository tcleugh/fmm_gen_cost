#ifndef FMM_MM_POLYMATCH_POLY_CANDIDATE_HPP_
#define FMM_MM_POLYMATCH_POLY_CANDIDATE_HPP_

#include <vector>

#include "core/geometry.hpp"
#include "network/polygon_layer.hpp"
#include "network/type.hpp"

namespace FMM {
namespace MM {

enum class PolyCandidateKind { Link, Polygon };

// A unified GPS-point candidate for polymatch. Link candidates carry an edge
// pointer and offset (same as existing FMM::MM::Candidate); polygon candidates
// carry a polygon index and the matched location inside or on the polygon. Both
// share `ep_distance` (distance from observed GPS to the matched location) for
// emission probability computation.
struct PolyCandidate {
  PolyCandidateKind kind = PolyCandidateKind::Link;

  // For both kinds:
  double ep_distance = 0.0;        // distance from observed point to matched point
  FMM::CORE::Point matched_point;  // matched location (on edge, inside polygon, or on boundary)

  // Link-only fields (valid when kind == Link)
  FMM::NETWORK::Edge *edge = nullptr;
  double offset = 0.0;  // arc-length offset along edge

  // Polygon-only fields (valid when kind == Polygon)
  FMM::NETWORK::PolygonIndex polygon_index = 0;
  bool inside = false;  // true if GPS observation is inside / on boundary of polygon

  bool is_link() const { return kind == PolyCandidateKind::Link; }
  bool is_polygon() const { return kind == PolyCandidateKind::Polygon; }
};

using PolyPointCandidates = std::vector<PolyCandidate>;
using PolyTrajCandidates = std::vector<PolyPointCandidates>;

} // MM
} // FMM

#endif // FMM_MM_POLYMATCH_POLY_CANDIDATE_HPP_
