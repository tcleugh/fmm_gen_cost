#include "network/polygon_layer.hpp"

#include <ogrsf_frmts.h>
#include <boost/geometry.hpp>
#include <stdexcept>
#include <sstream>

#include "util/debug.hpp"

namespace FMM {
namespace NETWORK {

namespace bg = boost::geometry;

bool PolygonLayer::load(const CONFIG::PolygonConfig &config) {
  polygons_.clear();
  id_to_index_.clear();
  rtree_.clear();

  if (config.file.empty()) {
    SPDLOG_INFO("PolygonLayer: no file configured, skipping load");
    return true;
  }

  GDALAllRegister();
  GDALDataset *ds = static_cast<GDALDataset *>(
      GDALOpenEx(config.file.c_str(), GDAL_OF_VECTOR, nullptr, nullptr,
                 nullptr));
  if (!ds) {
    SPDLOG_CRITICAL("Failed to open polygon shapefile {}", config.file);
    return false;
  }

  OGRLayer *layer = ds->GetLayer(0);
  if (!layer) {
    SPDLOG_CRITICAL("Polygon shapefile contains no layer");
    GDALClose(ds);
    return false;
  }

  OGRFeatureDefn *defn = layer->GetLayerDefn();
  int id_field = defn->GetFieldIndex(config.id_name.c_str());
  int cost_field = defn->GetFieldIndex(config.cost_name.c_str());
  if (id_field < 0) {
    SPDLOG_CRITICAL("Polygon ID field {} not found", config.id_name);
    GDALClose(ds);
    return false;
  }

  layer->ResetReading();
  OGRFeature *feature = nullptr;
  while ((feature = layer->GetNextFeature()) != nullptr) {
    PolygonID id = feature->GetFieldAsInteger64(id_field);
    if (id == 0) {
      SPDLOG_CRITICAL("Polygon feature has invalid ID 0 (FR-018)");
      OGRFeature::DestroyFeature(feature);
      GDALClose(ds);
      return false;
    }
    if (id_to_index_.find(id) != id_to_index_.end()) {
      SPDLOG_CRITICAL("Duplicate polygon ID {} (FR-018)", id);
      OGRFeature::DestroyFeature(feature);
      GDALClose(ds);
      return false;
    }

    OGRGeometry *ogr_geom = feature->GetGeometryRef();
    if (!ogr_geom ||
        wkbFlatten(ogr_geom->getGeometryType()) != wkbPolygon) {
      SPDLOG_WARN("Polygon feature ID {} has invalid geometry type, skipping",
                  id);
      OGRFeature::DestroyFeature(feature);
      continue;
    }
    OGRPolygon *ogr_poly = ogr_geom->toPolygon();

    Polygon poly;
    poly.id = id;
    poly.weight = (cost_field >= 0) ? feature->GetFieldAsDouble(cost_field)
                                    : 1.0;

    OGRLinearRing *ring = ogr_poly->getExteriorRing();
    if (!ring) {
      SPDLOG_WARN("Polygon feature ID {} has no exterior ring, skipping", id);
      OGRFeature::DestroyFeature(feature);
      continue;
    }
    int n_pts = ring->getNumPoints();
    auto &outer = poly.geom.outer();
    outer.reserve(n_pts);
    for (int i = 0; i < n_pts; ++i) {
      outer.emplace_back(ring->getX(i), ring->getY(i));
    }
    bg::correct(poly.geom);
    if (!bg::is_valid(poly.geom)) {
      SPDLOG_WARN("Polygon feature ID {} is invalid (self-intersect or zero area), skipping (FR-013)",
                  id);
      OGRFeature::DestroyFeature(feature);
      continue;
    }
    if (poly.weight < 0) {
      SPDLOG_WARN("Polygon ID {} has negative weight, skipping", id);
      OGRFeature::DestroyFeature(feature);
      continue;
    }

    poly.index = static_cast<PolygonIndex>(polygons_.size());
    poly.bbox = bg::return_envelope<boost_box>(poly.geom);
    id_to_index_[id] = poly.index;
    rtree_.insert(std::make_pair(poly.bbox, poly.index));
    polygons_.push_back(std::move(poly));

    OGRFeature::DestroyFeature(feature);
  }

  GDALClose(ds);
  SPDLOG_INFO("Loaded {} polygons from {}", polygons_.size(), config.file);
  return true;
}

std::vector<PolygonIndex> PolygonLayer::polygons_containing(
    const FMM::CORE::Point &p) const {
  std::vector<PolygonIndex> out;
  std::vector<RTreeItem> hits;
  rtree_.query(boost::geometry::index::intersects(p),
               std::back_inserter(hits));
  for (const auto &hit : hits) {
    const Polygon &poly = polygons_[hit.second];
    if (bg::covered_by(p, poly.geom)) {
      out.push_back(hit.second);
    }
  }
  return out;
}

std::vector<PolygonIndex> PolygonLayer::polygons_within_radius(
    const FMM::CORE::Point &p, double radius) const {
  std::vector<PolygonIndex> out;
  boost_box query_box(
      FMM::CORE::Point(bg::get<0>(p) - radius, bg::get<1>(p) - radius),
      FMM::CORE::Point(bg::get<0>(p) + radius, bg::get<1>(p) + radius));
  std::vector<RTreeItem> hits;
  rtree_.query(boost::geometry::index::intersects(query_box),
               std::back_inserter(hits));
  for (const auto &hit : hits) {
    const Polygon &poly = polygons_[hit.second];
    if (bg::distance(p, poly.geom) <= radius) {
      out.push_back(hit.second);
    }
  }
  return out;
}

double PolygonLayer::min_boundary_distance(PolygonIndex idx,
                                           const FMM::CORE::Point &p) const {
  const Polygon &poly = polygons_[idx];
  if (bg::covered_by(p, poly.geom)) {
    return 0.0;
  }
  return bg::distance(p, poly.geom);
}

} // NETWORK
} // FMM
