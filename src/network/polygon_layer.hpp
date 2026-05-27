#ifndef FMM_NETWORK_POLYGON_LAYER_HPP_
#define FMM_NETWORK_POLYGON_LAYER_HPP_

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "core/geometry.hpp"
#include "config/polygon_config.hpp"

namespace FMM {
namespace NETWORK {

typedef int64_t PolygonID;
typedef uint32_t PolygonIndex;

struct Polygon {
  PolygonIndex index;
  PolygonID id;
  boost::geometry::model::polygon<FMM::CORE::Point> geom;
  double weight = 1.0;
  boost::geometry::model::box<FMM::CORE::Point> bbox;
};

class PolygonLayer {
public:
  typedef boost::geometry::model::box<FMM::CORE::Point> boost_box;
  typedef std::pair<boost_box, PolygonIndex> RTreeItem;
  typedef boost::geometry::index::rtree<
      RTreeItem, boost::geometry::index::quadratic<16> > Rtree;

  PolygonLayer() = default;

  bool load(const CONFIG::PolygonConfig &config);

  const std::vector<Polygon> &polygons() const { return polygons_; }
  size_t size() const { return polygons_.size(); }
  bool empty() const { return polygons_.empty(); }

  PolygonIndex id_to_index_or_throw(PolygonID id) const {
    return id_to_index_.at(id);
  }
  bool has_id(PolygonID id) const {
    return id_to_index_.find(id) != id_to_index_.end();
  }

  std::vector<PolygonIndex> polygons_containing(const FMM::CORE::Point &p) const;
  std::vector<PolygonIndex> polygons_within_radius(
      const FMM::CORE::Point &p, double radius) const;
  double min_boundary_distance(PolygonIndex idx,
                               const FMM::CORE::Point &p) const;

private:
  std::vector<Polygon> polygons_;
  std::unordered_map<PolygonID, PolygonIndex> id_to_index_;
  Rtree rtree_;
};

} // NETWORK
} // FMM

#endif // FMM_NETWORK_POLYGON_LAYER_HPP_
