#include "network/access_point_layer.hpp"

#include <ogrsf_frmts.h>
#include <boost/geometry.hpp>
#include <unordered_map>
#include <sstream>
#include <cmath>

#include "util/debug.hpp"

namespace FMM {
namespace NETWORK {

namespace bg = boost::geometry;

bool AccessPointLayer::load(const CONFIG::AccessPointConfig &config,
                            const PolygonLayer &polygon_layer,
                            const Network &network,
                            double boundary_epsilon) {
  access_points_.clear();
  node_id_to_index_.clear();
  polygon_to_aps_.clear();

  if (config.file.empty()) {
    SPDLOG_INFO("AccessPointLayer: no file configured, skipping load");
    return true;
  }

  GDALAllRegister();
  GDALDataset *ds = static_cast<GDALDataset *>(
      GDALOpenEx(config.file.c_str(), GDAL_OF_VECTOR, nullptr, nullptr,
                 nullptr));
  if (!ds) {
    SPDLOG_CRITICAL("Failed to open access point shapefile {}", config.file);
    return false;
  }

  OGRLayer *layer = ds->GetLayer(0);
  if (!layer) {
    SPDLOG_CRITICAL("Access point shapefile has no layer");
    GDALClose(ds);
    return false;
  }

  OGRFeatureDefn *defn = layer->GetLayerDefn();
  int node_field = defn->GetFieldIndex(config.node_id_name.c_str());
  int poly_field = defn->GetFieldIndex(config.polygon_id_name.c_str());
  if (node_field < 0 || poly_field < 0) {
    SPDLOG_CRITICAL("Access point fields missing: node_id={} polygon_id={}",
                    config.node_id_name, config.polygon_id_name);
    GDALClose(ds);
    return false;
  }

  std::vector<AccessPointFeature> features;
  layer->ResetReading();
  OGRFeature *feature = nullptr;
  while ((feature = layer->GetNextFeature()) != nullptr) {
    NodeID node_id = feature->GetFieldAsInteger64(node_field);
    PolygonID polygon_id = feature->GetFieldAsInteger64(poly_field);

    OGRGeometry *ogr_geom = feature->GetGeometryRef();
    if (!ogr_geom || wkbFlatten(ogr_geom->getGeometryType()) != wkbPoint) {
      SPDLOG_CRITICAL("Access point feature has invalid geometry type");
      OGRFeature::DestroyFeature(feature);
      GDALClose(ds);
      return false;
    }
    OGRPoint *ogr_pt = ogr_geom->toPoint();

    AccessPointFeature f;
    f.node_id = node_id;
    f.polygon_id = polygon_id;
    f.point = FMM::CORE::Point(ogr_pt->getX(), ogr_pt->getY());

    if (!polygon_layer.has_id(polygon_id)) {
      SPDLOG_CRITICAL(
          "Access point references unknown polygon ID {} (FR-005)", polygon_id);
      OGRFeature::DestroyFeature(feature);
      GDALClose(ds);
      return false;
    }
    PolygonIndex p_idx = polygon_layer.id_to_index_or_throw(polygon_id);
    const Polygon &poly = polygon_layer.polygons()[p_idx];
    double d = bg::distance(f.point, poly.geom);
    if (d > boundary_epsilon) {
      SPDLOG_CRITICAL(
          "Access point node_id {} is {} from polygon {} boundary, exceeds epsilon {} (FR-005)",
          node_id, d, polygon_id, boundary_epsilon);
      OGRFeature::DestroyFeature(feature);
      GDALClose(ds);
      return false;
    }

    features.push_back(f);
    OGRFeature::DestroyFeature(feature);
  }
  GDALClose(ds);

  std::unordered_map<NodeID, std::vector<size_t>> by_node;
  for (size_t i = 0; i < features.size(); ++i) {
    by_node[features[i].node_id].push_back(i);
  }

  for (const auto &kv : by_node) {
    NodeID node_id = kv.first;
    const auto &idxs = kv.second;
    const auto &first = features[idxs[0]];
    for (size_t k = 1; k < idxs.size(); ++k) {
      const auto &other = features[idxs[k]];
      double dx = bg::get<0>(first.point) - bg::get<0>(other.point);
      double dy = bg::get<1>(first.point) - bg::get<1>(other.point);
      if (std::sqrt(dx * dx + dy * dy) > boundary_epsilon) {
        SPDLOG_CRITICAL(
            "Access point node_id {} has contradictory geometries across features (FR-005)",
            node_id);
        return false;
      }
    }

    AccessPoint ap;
    ap.index = static_cast<AccessPointIndex>(access_points_.size());
    ap.node_id = node_id;
    ap.point = first.point;
    for (size_t i : idxs) {
      PolygonIndex p_idx =
          polygon_layer.id_to_index_or_throw(features[i].polygon_id);
      ap.polygons.push_back(p_idx);
    }

    try {
      NodeIndex n_idx = network.get_node_index(node_id);
      ap.attached_node = n_idx;
      const auto &edges = network.get_edges();
      for (const auto &e : edges) {
        if (e.source == n_idx || e.target == n_idx) {
          ap.attached_edges.push_back(e.index);
        }
      }
    } catch (const std::out_of_range &) {
      // node_id not in network — must be polygon-shared (>= 2 polygons)
      ap.attached_node = boost::none;
    }

    if (!ap.attached_node.has_value() && ap.polygons.size() < 2) {
      SPDLOG_CRITICAL(
          "Access point node_id {} is neither link-attached nor shared between polygons (FR-004)",
          node_id);
      return false;
    }

    AccessPointIndex idx = ap.index;
    node_id_to_index_[node_id] = idx;
    for (PolygonIndex p : ap.polygons) {
      polygon_to_aps_[p].push_back(idx);
    }
    access_points_.push_back(std::move(ap));
  }

  // Warn for polygons with no access points (FR-014)
  for (PolygonIndex p = 0; p < polygon_layer.size(); ++p) {
    if (polygon_to_aps_.find(p) == polygon_to_aps_.end()) {
      SPDLOG_WARN("Polygon ID {} has no access points; excluded from candidates (FR-014)",
                  polygon_layer.polygons()[p].id);
    }
  }

  SPDLOG_INFO("Loaded {} access points from {}", access_points_.size(),
              config.file);
  return true;
}

} // NETWORK
} // FMM
