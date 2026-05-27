#ifndef FMM_NETWORK_ACCESS_POINT_LAYER_HPP_
#define FMM_NETWORK_ACCESS_POINT_LAYER_HPP_

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "core/geometry.hpp"
#include "network/type.hpp"
#include "network/polygon_layer.hpp"
#include "network/network.hpp"
#include "config/access_point_config.hpp"

namespace FMM {
namespace NETWORK {

typedef uint32_t AccessPointIndex;

struct AccessPointFeature {
  NodeID node_id;
  PolygonID polygon_id;
  FMM::CORE::Point point;
};

struct AccessPoint {
  AccessPointIndex index;
  NodeID node_id;
  FMM::CORE::Point point;
  std::vector<PolygonIndex> polygons;
  boost::optional<NodeIndex> attached_node;
  std::vector<EdgeIndex> attached_edges;
};

class AccessPointLayer {
public:
  AccessPointLayer() = default;

  bool load(const CONFIG::AccessPointConfig &config,
            const PolygonLayer &polygon_layer,
            const Network &network,
            double boundary_epsilon = 1e-6);

  const std::vector<AccessPoint> &access_points() const {
    return access_points_;
  }
  size_t size() const { return access_points_.size(); }

  bool has_node_id(NodeID node_id) const {
    return node_id_to_index_.find(node_id) != node_id_to_index_.end();
  }
  AccessPointIndex node_id_to_index(NodeID node_id) const {
    return node_id_to_index_.at(node_id);
  }

  const std::vector<AccessPointIndex> &aps_for_polygon(PolygonIndex p) const {
    static const std::vector<AccessPointIndex> kEmpty;
    auto it = polygon_to_aps_.find(p);
    return it == polygon_to_aps_.end() ? kEmpty : it->second;
  }

private:
  std::vector<AccessPoint> access_points_;
  std::unordered_map<NodeID, AccessPointIndex> node_id_to_index_;
  std::unordered_map<PolygonIndex, std::vector<AccessPointIndex>>
      polygon_to_aps_;
};

} // NETWORK
} // FMM

#endif // FMM_NETWORK_ACCESS_POINT_LAYER_HPP_
