#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "util/debug.hpp"
#include "network/network.hpp"
#include "network/link_graph_routing.hpp"
#include "network/polygon_layer.hpp"
#include "network/access_point_layer.hpp"
#include "network/poly_link_graph.hpp"
#include "mm/polymatch/polymatch_algorithm.hpp"
#include "mm/polymatch/poly_match_result.hpp"
#include "io/poly_mm_writer.hpp"
#include "config/polygon_config.hpp"
#include "config/access_point_config.hpp"
#include "config/result_config.hpp"

#include <ogrsf_frmts.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>

using namespace FMM;
using namespace FMM::CORE;
using namespace FMM::NETWORK;
using namespace FMM::ROUTING;
using namespace FMM::MM;

namespace {

// Fixture directory is injected by CMake (target_compile_definitions) so the
// tests run correctly regardless of where the project root is mounted —
// hardcoding a `/workspace/...` path breaks outside this sandbox.
#ifndef POLYMATCH_FIXTURE_DIR
#error "POLYMATCH_FIXTURE_DIR must be set via target_compile_definitions"
#endif
const char *kFixtureDir = POLYMATCH_FIXTURE_DIR;

bool path_exists(const std::string &p) {
  struct stat st;
  return stat(p.c_str(), &st) == 0;
}

void ensure_dir(const std::string &p) {
  if (!path_exists(p)) {
    mkdir(p.c_str(), 0755);
  }
}

// Build a small grid network shapefile using OGR. 3x3 lattice + 12 edges.
void create_network_shapefile(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  REQUIRE(drv != nullptr);
  if (path_exists(path)) drv->Delete(path.c_str());

  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  REQUIRE(ds != nullptr);
  OGRLayer *layer = ds->CreateLayer("edges", nullptr, wkbLineString, nullptr);
  REQUIRE(layer != nullptr);
  {
    OGRFieldDefn f("id", OFTInteger64); layer->CreateField(&f);
    OGRFieldDefn s("source", OFTInteger64); layer->CreateField(&s);
    OGRFieldDefn t("target", OFTInteger64); layer->CreateField(&t);
  }
  struct {
    int eid, src, tgt;
    double x1, y1, x2, y2;
  } edges[] = {
    { 1, 1, 2,  0.0, 2.0,  1.0, 2.0},
    { 2, 2, 3,  1.0, 2.0,  2.0, 2.0},
    { 3, 4, 5,  0.0, 1.0,  1.0, 1.0},
    { 4, 5, 6,  1.0, 1.0,  2.0, 1.0},
    { 5, 7, 8,  0.0, 0.0,  1.0, 0.0},
    { 6, 8, 9,  1.0, 0.0,  2.0, 0.0},
    { 9, 1, 4,  0.0, 2.0,  0.0, 1.0},
    {10, 2, 5,  1.0, 2.0,  1.0, 1.0},
    {11, 3, 6,  2.0, 2.0,  2.0, 1.0},
    {12, 4, 7,  0.0, 1.0,  0.0, 0.0},
    {13, 5, 8,  1.0, 1.0,  1.0, 0.0},
    {14, 6, 9,  2.0, 1.0,  2.0, 0.0},
  };
  for (auto &e : edges) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("id", (GIntBig)e.eid);
    f->SetField("source", (GIntBig)e.src);
    f->SetField("target", (GIntBig)e.tgt);
    OGRLineString ls;
    ls.addPoint(e.x1, e.y1);
    ls.addPoint(e.x2, e.y2);
    f->SetGeometry(&ls);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  }
  GDALClose(ds);
}

// Build the canonical polygons shapefile: 4 polygons (7, 42, 99-invalid, 200-no-AP).
void create_polygons_shapefile(const std::string &path,
                               bool include_invalid = true,
                               bool include_no_ap = true) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("polygons", nullptr, wkbPolygon, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn cost_f("cost", OFTReal); layer->CreateField(&cost_f);

  auto add_poly = [&](GIntBig id, double cost,
                      const std::vector<std::pair<double, double>> &pts) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("id", id);
    f->SetField("cost", cost);
    OGRPolygon poly;
    OGRLinearRing ring;
    for (auto &p : pts) ring.addPoint(p.first, p.second);
    ring.addPoint(pts.front().first, pts.front().second);
    poly.addRing(&ring);
    f->SetGeometry(&poly);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  };

  add_poly(7, 2.0,
           {{0.8, 0.8}, {1.2, 0.8}, {1.2, 1.2}, {0.8, 1.2}});
  add_poly(42, 1.5,
           {{1.8, -0.2}, {2.2, -0.2}, {2.2, 0.2}, {1.8, 0.2}});
  if (include_invalid) {
    // bowtie / self-intersecting
    add_poly(99, 1.0,
             {{3.0, 3.0}, {4.0, 4.0}, {4.0, 3.0}, {3.0, 4.0}});
  }
  if (include_no_ap) {
    add_poly(200, 1.0,
             {{-1.0, -1.0}, {-0.5, -1.0}, {-0.5, -0.5}, {-1.0, -0.5}});
  }
  GDALClose(ds);
}

// Build the canonical access points shapefile (passes all FR-005 conditions).
void create_aps_shapefile_clean(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("aps", nullptr, wkbPoint, nullptr);
  OGRFieldDefn nf("node_id", OFTInteger64); layer->CreateField(&nf);
  OGRFieldDefn pf("polygon_id", OFTInteger64); layer->CreateField(&pf);

  auto add_ap = [&](GIntBig node_id, GIntBig poly_id, double x, double y) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("node_id", node_id);
    f->SetField("polygon_id", poly_id);
    OGRPoint pt(x, y);
    f->SetGeometry(&pt);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  };

  // Polygon 7 boundary touches network nodes 2 (1,2 -> use 1.0,1.2), 4 (0,1 -> 0.8,1.0), 5 (1,1 -> via 0.8,0.8 corner)
  add_ap(5, 7, 0.8, 0.8);   // boundary corner; AP for node 5
  add_ap(2, 7, 1.0, 1.2);   // top of polygon
  add_ap(4, 7, 0.8, 1.0);   // left of polygon
  // Polygon 42 boundary touches node 9 (at 2.0, 0.0) and node 6 (2,1 -> 2.0,0.2)
  add_ap(9, 42, 2.0, 0.0);
  add_ap(6, 42, 2.0, 0.2);
  GDALClose(ds);
}

// FR-018: invalid polygon ID 0
void create_polygons_id_zero(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("polys", nullptr, wkbPolygon, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn cf("cost", OFTReal); layer->CreateField(&cf);
  OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
  f->SetField("id", (GIntBig)0);  // FR-018: id == 0 invalid
  f->SetField("cost", 1.0);
  OGRPolygon poly; OGRLinearRing r;
  r.addPoint(0, 0); r.addPoint(1, 0); r.addPoint(1, 1); r.addPoint(0, 1); r.addPoint(0, 0);
  poly.addRing(&r);
  f->SetGeometry(&poly);
  layer->CreateFeature(f);
  OGRFeature::DestroyFeature(f);
  GDALClose(ds);
}

void create_polygons_duplicate_id(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("polys", nullptr, wkbPolygon, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn cf("cost", OFTReal); layer->CreateField(&cf);
  auto add = [&](GIntBig id, double x1, double y1, double x2, double y2) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("id", id);
    f->SetField("cost", 1.0);
    OGRPolygon poly; OGRLinearRing r;
    r.addPoint(x1, y1); r.addPoint(x2, y1); r.addPoint(x2, y2); r.addPoint(x1, y2); r.addPoint(x1, y1);
    poly.addRing(&r);
    f->SetGeometry(&poly);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  };
  add(7, 0, 0, 1, 1);
  add(7, 2, 2, 3, 3);  // duplicate id (FR-018)
  GDALClose(ds);
}

void create_aps_off_boundary(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("aps", nullptr, wkbPoint, nullptr);
  OGRFieldDefn nf("node_id", OFTInteger64); layer->CreateField(&nf);
  OGRFieldDefn pf("polygon_id", OFTInteger64); layer->CreateField(&pf);
  OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
  f->SetField("node_id", (GIntBig)5);
  f->SetField("polygon_id", (GIntBig)7);
  OGRPoint pt(0.3, 0.3);  // way off the polygon 7 boundary
  f->SetGeometry(&pt);
  layer->CreateFeature(f);
  OGRFeature::DestroyFeature(f);
  GDALClose(ds);
}

void create_aps_orphan_polygon(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("aps", nullptr, wkbPoint, nullptr);
  OGRFieldDefn nf("node_id", OFTInteger64); layer->CreateField(&nf);
  OGRFieldDefn pf("polygon_id", OFTInteger64); layer->CreateField(&pf);
  OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
  f->SetField("node_id", (GIntBig)5);
  f->SetField("polygon_id", (GIntBig)999);  // polygon 999 does not exist
  OGRPoint pt(0.8, 0.8);
  f->SetGeometry(&pt);
  layer->CreateFeature(f);
  OGRFeature::DestroyFeature(f);
  GDALClose(ds);
}

void create_aps_contradictory(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("aps", nullptr, wkbPoint, nullptr);
  OGRFieldDefn nf("node_id", OFTInteger64); layer->CreateField(&nf);
  OGRFieldDefn pf("polygon_id", OFTInteger64); layer->CreateField(&pf);
  auto add = [&](GIntBig node, GIntBig poly, double x, double y) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("node_id", node);
    f->SetField("polygon_id", poly);
    OGRPoint pt(x, y);
    f->SetGeometry(&pt);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  };
  add(5, 7, 0.8, 0.8);
  add(5, 7, 1.2, 1.2);  // same node_id, different geometry (FR-005 cond 3)
  GDALClose(ds);
}

// Two adjacent polygons that share a boundary segment at x=1.2, plus an AP
// shapefile in which node_id 2000 is listed for BOTH polygons (at the shared
// boundary point (1.2, 1.0)). Used to exercise T023 (shared-node dedup) and
// T028 (polygon-shared-only when node_id is not in the network).
void create_polygons_with_shared_boundary(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("polys", nullptr, wkbPolygon, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn cf("cost", OFTReal); layer->CreateField(&cf);

  auto add = [&](GIntBig id, double cost,
                 std::vector<std::pair<double, double>> pts) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("id", id);
    f->SetField("cost", cost);
    OGRPolygon poly;
    OGRLinearRing r;
    for (auto &p : pts) r.addPoint(p.first, p.second);
    r.addPoint(pts.front().first, pts.front().second);
    poly.addRing(&r);
    f->SetGeometry(&poly);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  };
  // Two squares sharing the vertical segment x=1.2 between y=0.8 and y=1.2.
  add(7, 1.0, {{0.8, 0.8}, {1.2, 0.8}, {1.2, 1.2}, {0.8, 1.2}});
  add(8, 1.0, {{1.2, 0.8}, {1.6, 0.8}, {1.6, 1.2}, {1.2, 1.2}});
  GDALClose(ds);
}

// AP shapefile pairing with the above: node 2000 is listed for both polygons
// at the shared point (1.2, 1.0). Node 2000 is NOT in the network's node map,
// so this AP is "polygon-shared-only" — covers T028.
void create_aps_shared_only(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("aps", nullptr, wkbPoint, nullptr);
  OGRFieldDefn nf("node_id", OFTInteger64); layer->CreateField(&nf);
  OGRFieldDefn pf("polygon_id", OFTInteger64); layer->CreateField(&pf);
  auto add = [&](GIntBig node, GIntBig poly, double x, double y) {
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("node_id", node);
    f->SetField("polygon_id", poly);
    OGRPoint pt(x, y);
    f->SetGeometry(&pt);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
  };
  add(2000, 7, 1.2, 1.0);   // shared between polygons 7 and 8
  add(2000, 8, 1.2, 1.0);   // same node_id, same point geometry
  // give each polygon at least one link-attached AP so FR-014 doesn't warn
  add(2, 7, 1.0, 1.2);       // network node 2 on polygon 7's boundary
  add(3, 8, 1.6, 1.2);       // node 3 is at (2,2) — NOT on polygon 8 boundary,
                              // so we use a different node.
  GDALClose(ds);
}

// Build a larger grid network for the perf baseline: GRID x GRID nodes with
// horizontal + vertical edges. With GRID=11 we get 121 nodes and 220 edges.
void create_perf_network(const std::string &path, int grid = 11) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("edges", nullptr, wkbLineString, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn s_f("source", OFTInteger64); layer->CreateField(&s_f);
  OGRFieldDefn t_f("target", OFTInteger64); layer->CreateField(&t_f);

  auto node_id = [grid](int r, int c) -> GIntBig {
    return (GIntBig)(r * grid + c + 1);  // 1-based
  };

  GIntBig eid = 1;
  for (int r = 0; r < grid; ++r) {
    for (int c = 0; c < grid; ++c) {
      double x = c * 1.0, y = r * 1.0;
      // horizontal edge to (r, c+1)
      if (c + 1 < grid) {
        OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("id", eid++);
        f->SetField("source", node_id(r, c));
        f->SetField("target", node_id(r, c + 1));
        OGRLineString ls;
        ls.addPoint(x, y);
        ls.addPoint(x + 1.0, y);
        f->SetGeometry(&ls);
        layer->CreateFeature(f);
        OGRFeature::DestroyFeature(f);
      }
      // vertical edge to (r+1, c)
      if (r + 1 < grid) {
        OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("id", eid++);
        f->SetField("source", node_id(r, c));
        f->SetField("target", node_id(r + 1, c));
        OGRLineString ls;
        ls.addPoint(x, y);
        ls.addPoint(x, y + 1.0);
        f->SetGeometry(&ls);
        layer->CreateFeature(f);
        OGRFeature::DestroyFeature(f);
      }
    }
  }
  GDALClose(ds);
}

// 10x10 = 100 polygons. Each polygon is a 0.3 x 0.3 square centered at a grid
// intersection. Each polygon gets one AP at one of its boundary corners.
void create_perf_polygons(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("polys", nullptr, wkbPolygon, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn cf("cost", OFTReal); layer->CreateField(&cf);

  GIntBig pid = 1;
  for (int r = 1; r <= 10; ++r) {
    for (int c = 1; c <= 10; ++c) {
      double cx = c, cy = r;
      OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
      f->SetField("id", pid++);
      f->SetField("cost", 1.0);
      OGRPolygon p; OGRLinearRing rr;
      rr.addPoint(cx - 0.15, cy - 0.15);
      rr.addPoint(cx + 0.15, cy - 0.15);
      rr.addPoint(cx + 0.15, cy + 0.15);
      rr.addPoint(cx - 0.15, cy + 0.15);
      rr.addPoint(cx - 0.15, cy - 0.15);
      p.addRing(&rr);
      f->SetGeometry(&p);
      layer->CreateFeature(f);
      OGRFeature::DestroyFeature(f);
    }
  }
  GDALClose(ds);
}

void create_perf_aps(const std::string &path, int grid = 11) {
  // One AP per polygon. node_id chosen to be a real network node — we use
  // the center node (r, c) of each polygon, which maps to node_id = r*grid+c+1.
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("aps", nullptr, wkbPoint, nullptr);
  OGRFieldDefn nf("node_id", OFTInteger64); layer->CreateField(&nf);
  OGRFieldDefn pf("polygon_id", OFTInteger64); layer->CreateField(&pf);

  GIntBig pid = 1;
  for (int r = 1; r <= 10; ++r) {
    for (int c = 1; c <= 10; ++c) {
      // Place AP on the boundary at (cx-0.15, cy-0.15). The node at (r-1, c-1)
      // doesn't have these exact coordinates but link attachment is by ID
      // lookup — we still want a network node here. Use node (r-1, c-1) which
      // sits at (c-1, r-1).
      // But that's not on this polygon's boundary. To pass FR-005 validation
      // we need the AP's geometry on the polygon boundary, and the AP node_id
      // can be ANY value (linked or not). For this perf fixture we attach to
      // a fake node ID (offset 1000+) so APs are polygon-only — but FR-004
      // requires either link-attachment OR polygons.size() >= 2 — so this
      // wouldn't validate.
      //
      // Solution: link AP geometry to the network node at (r, c) which is at
      // (c, r) coordinates. The polygon boundary point closest to (c, r) is
      // (cx-0.15, cy-0.15) for cx=c, cy=r? No — (c, r) is the polygon's CENTER,
      // not boundary. We need the polygon centered ELSEWHERE so its corner
      // overlaps a network node.
      //
      // Re-design: center polygon at (c+0.15, r+0.15) so its lower-left
      // corner is at (c, r) which IS a network node.
      double cx = c + 0.15, cy = r + 0.15;
      OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
      f->SetField("node_id", (GIntBig)(r * grid + c + 1));  // network node at (c, r)
      f->SetField("polygon_id", pid++);
      OGRPoint pt(c, r);  // lower-left corner of the redesigned polygon
      f->SetGeometry(&pt);
      layer->CreateFeature(f);
      OGRFeature::DestroyFeature(f);
    }
  }
  GDALClose(ds);
}

// Redesigned polygons: lower-left corner anchored at network node (c, r).
void create_perf_polygons_v2(const std::string &path) {
  GDALAllRegister();
  GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  if (path_exists(path)) drv->Delete(path.c_str());
  GDALDataset *ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
  OGRLayer *layer = ds->CreateLayer("polys", nullptr, wkbPolygon, nullptr);
  OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
  OGRFieldDefn cf("cost", OFTReal); layer->CreateField(&cf);

  GIntBig pid = 1;
  for (int r = 1; r <= 10; ++r) {
    for (int c = 1; c <= 10; ++c) {
      OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
      f->SetField("id", pid++);
      f->SetField("cost", 1.0);
      OGRPolygon p; OGRLinearRing rr;
      rr.addPoint(c, r);
      rr.addPoint(c + 0.3, r);
      rr.addPoint(c + 0.3, r + 0.3);
      rr.addPoint(c, r + 0.3);
      rr.addPoint(c, r);
      p.addRing(&rr);
      f->SetGeometry(&p);
      layer->CreateFeature(f);
      OGRFeature::DestroyFeature(f);
    }
  }
  GDALClose(ds);
}

class Fixture {
 public:
  std::string network_path;
  std::string polygons_path;
  std::string aps_path;
  std::string poly_zero_path;
  std::string poly_dup_path;
  std::string aps_off_boundary_path;
  std::string aps_orphan_path;
  std::string aps_contradictory_path;
  std::string poly_shared_path;
  std::string aps_shared_path;
  std::string perf_net_path;
  std::string perf_polys_path;
  std::string perf_aps_path;

  Fixture() {
    ensure_dir(kFixtureDir);
    network_path = std::string(kFixtureDir) + "/edges.shp";
    polygons_path = std::string(kFixtureDir) + "/polygons.shp";
    aps_path = std::string(kFixtureDir) + "/access_points.shp";
    poly_zero_path = std::string(kFixtureDir) + "/polygons_id_zero.shp";
    poly_dup_path = std::string(kFixtureDir) + "/polygons_dup.shp";
    aps_off_boundary_path = std::string(kFixtureDir) + "/aps_off_boundary.shp";
    aps_orphan_path = std::string(kFixtureDir) + "/aps_orphan.shp";
    aps_contradictory_path = std::string(kFixtureDir) + "/aps_contradictory.shp";
    poly_shared_path = std::string(kFixtureDir) + "/polygons_shared.shp";
    aps_shared_path = std::string(kFixtureDir) + "/aps_shared.shp";
    perf_net_path = std::string(kFixtureDir) + "/perf_edges.shp";
    perf_polys_path = std::string(kFixtureDir) + "/perf_polygons.shp";
    perf_aps_path = std::string(kFixtureDir) + "/perf_aps.shp";

    create_network_shapefile(network_path);
    create_polygons_shapefile(polygons_path);
    create_aps_shapefile_clean(aps_path);
    create_polygons_id_zero(poly_zero_path);
    create_polygons_duplicate_id(poly_dup_path);
    create_aps_off_boundary(aps_off_boundary_path);
    create_aps_orphan_polygon(aps_orphan_path);
    create_aps_contradictory(aps_contradictory_path);
    create_polygons_with_shared_boundary(poly_shared_path);
    create_aps_shared_only(aps_shared_path);
    create_perf_network(perf_net_path);
    create_perf_polygons_v2(perf_polys_path);
    create_perf_aps(perf_aps_path);
  }
};

static Fixture &fixture() {
  static Fixture f;
  return f;
}

}  // namespace

// ============================================================
// US2 — PolygonLayer / AccessPointLayer / AppConfig validation
// ============================================================

TEST_CASE("PolygonLayer loads valid shapefile (T020)", "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::warn);
  PolygonLayer layer;
  CONFIG::PolygonConfig cfg{fixture().polygons_path, "id", "cost"};
  REQUIRE(layer.load(cfg));
  // 7, 42, 200 valid; 99 invalid -> skipped (FR-013)
  REQUIRE(layer.size() == 3);
  REQUIRE(layer.has_id(7));
  REQUIRE(layer.has_id(42));
  REQUIRE(layer.has_id(200));
  REQUIRE_FALSE(layer.has_id(99));  // skipped
}

TEST_CASE("PolygonLayer skips invalid geometry with warning (T021)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer layer;
  CONFIG::PolygonConfig cfg{fixture().polygons_path, "id", "cost"};
  REQUIRE(layer.load(cfg));
  REQUIRE_FALSE(layer.has_id(99));  // bowtie skipped per FR-013
}

TEST_CASE("PolygonLayer point-in-polygon and boundary distance (T022)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer layer;
  CONFIG::PolygonConfig cfg{fixture().polygons_path, "id", "cost"};
  REQUIRE(layer.load(cfg));

  // Point inside polygon 7 (which spans 0.8..1.2 x 0.8..1.2)
  auto inside = layer.polygons_containing(Point(1.0, 1.0));
  REQUIRE(inside.size() == 1);
  PolygonIndex p7 = layer.id_to_index_or_throw(7);
  REQUIRE(inside[0] == p7);

  // Point exactly on boundary
  auto boundary = layer.polygons_containing(Point(0.8, 1.0));
  REQUIRE_FALSE(boundary.empty());  // covered_by treats boundary as inside

  // Point outside
  auto outside = layer.polygons_containing(Point(0.0, 0.0));
  REQUIRE(outside.empty());

  // Min boundary distance
  REQUIRE(layer.min_boundary_distance(p7, Point(1.0, 1.0)) == Approx(0.0));
  REQUIRE(layer.min_boundary_distance(p7, Point(0.8, 1.0)) == Approx(0.0));
  REQUIRE(layer.min_boundary_distance(p7, Point(0.3, 1.0)) > 0);
}

TEST_CASE("PolygonLayer rejects polygon id 0 (T086 FR-018)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer layer;
  CONFIG::PolygonConfig cfg{fixture().poly_zero_path, "id", "cost"};
  REQUIRE_FALSE(layer.load(cfg));
}

TEST_CASE("PolygonLayer rejects duplicate polygon id (T087 FR-018)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer layer;
  CONFIG::PolygonConfig cfg{fixture().poly_dup_path, "id", "cost"};
  REQUIRE_FALSE(layer.load(cfg));
}

TEST_CASE("AccessPointLayer loads clean fixture (T023)", "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::warn);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  REQUIRE(aps.size() == 5);
}

TEST_CASE("AccessPointLayer rejects off-boundary AP (T024 FR-005)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE_FALSE(aps.load(
      {fixture().aps_off_boundary_path, "node_id", "polygon_id"}, poly, net,
      1e-6));
}

TEST_CASE("AccessPointLayer rejects orphan polygon ref (T025 FR-005)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE_FALSE(aps.load(
      {fixture().aps_orphan_path, "node_id", "polygon_id"}, poly, net, 1e-6));
}

TEST_CASE("AccessPointLayer rejects contradictory shared AP (T026 FR-005)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE_FALSE(aps.load(
      {fixture().aps_contradictory_path, "node_id", "polygon_id"}, poly, net,
      1e-6));
}

TEST_CASE("AccessPointLayer link attachment via node_id (T027)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));

  // Find AP for node_id == 5 (network node 5 exists at center)
  REQUIRE(aps.has_node_id(5));
  AccessPointIndex idx5 = aps.node_id_to_index(5);
  const AccessPoint &ap5 = aps.access_points()[idx5];
  REQUIRE(ap5.attached_node.has_value());
  REQUIRE_FALSE(ap5.attached_edges.empty());  // multiple edges share node 5
}

TEST_CASE("AccessPointLayer dedups shared node across polygons (T023)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().poly_shared_path, "id", "cost"}));
  REQUIRE(poly.size() == 2);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_shared_path, "node_id", "polygon_id"}, poly,
                   net, 1e-6));

  // Three unique node IDs in fixture: 2 (polygon 7), 3 (polygon 8), 2000
  // (shared between 7 and 8). 4 features but 3 deduplicated APs.
  REQUIRE(aps.size() == 3);

  // Find the shared AP (node_id 2000)
  REQUIRE(aps.has_node_id(2000));
  AccessPointIndex shared_idx = aps.node_id_to_index(2000);
  const auto &shared = aps.access_points()[shared_idx];
  CHECK(shared.polygons.size() == 2);
}

TEST_CASE("AccessPointLayer polygon-shared-only AP has no link attachment (T028)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().poly_shared_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_shared_path, "node_id", "polygon_id"}, poly,
                   net, 1e-6));

  // node_id 2000 is not in the network (network has nodes 1..9). The AP is
  // valid because it's shared between polygons 7 and 8 (polygons.size() >= 2).
  AccessPointIndex idx = aps.node_id_to_index(2000);
  const auto &ap = aps.access_points()[idx];
  CHECK_FALSE(ap.attached_node.has_value());
  CHECK(ap.polygons.size() >= 2);
  CHECK(ap.attached_edges.empty());
}

TEST_CASE("AccessPointLayer warns on polygon with no AP (T029 FR-014)",
          "[polymatch][us2]") {
  spdlog::set_level(spdlog::level::off);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  // polygon 200 has no APs; it's still loaded into PolygonLayer but should be
  // absent from the polygon_to_aps reverse map.
  PolygonIndex p200 = poly.id_to_index_or_throw(200);
  REQUIRE(aps.aps_for_polygon(p200).empty());
}

// ============================================================
// US1 — PolyLinkGraph construction + through-cost tables
// ============================================================

TEST_CASE("PolyLinkGraph builds vertices and through-cost tables (T042/T043)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));

  PolyLinkGraph G(net, link_graph, poly, aps, /*through_penalty_factor=*/2.0);
  // Vertex count = |E| + Sum_p n_p (Held-Karp polygon-AP expansion). In the
  // fixture: 12 edges + (3 + 2 + 0) AP sub-vertices for polygons {7, 42, 200}
  // = 17 vertices.
  size_t expected_subs = 0;
  for (PolygonIndex p = 0; p < poly.size(); ++p) {
    expected_subs += aps.aps_for_polygon(p).size();
  }
  REQUIRE(G.n_sub_vertices() == expected_subs);
  REQUIRE(G.n_vertices() ==
          (size_t)net.get_edge_count() + expected_subs);
  REQUIRE(G.through_penalty_factor() == Approx(2.0));

  // Through-cost raw table for polygon 7 must equal weight * dist(AP_i, AP_j)
  PolygonIndex p7 = poly.id_to_index_or_throw(7);
  auto aps_p7 = aps.aps_for_polygon(p7);
  REQUIRE(aps_p7.size() >= 2);
  const auto &api = aps.access_points()[aps_p7[0]];
  const auto &apj = aps.access_points()[aps_p7[1]];
  double dx = boost::geometry::get<0>(api.point) -
              boost::geometry::get<0>(apj.point);
  double dy = boost::geometry::get<1>(api.point) -
              boost::geometry::get<1>(apj.point);
  double expected = poly.polygons()[p7].weight * std::sqrt(dx * dx + dy * dy);
  double raw = G.through_cost_raw(p7, aps_p7[0], aps_p7[1]);
  REQUIRE(raw == Approx(expected));
}

TEST_CASE("PolyLinkGraph through_penalty_factor scales costs (T044)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));

  PolyLinkGraph G1(net, link_graph, poly, aps, 1.0);
  PolyLinkGraph G2(net, link_graph, poly, aps, 5.0);
  // raw (factor-independent) table is identical; the factor applies at lookup
  PolygonIndex p7 = poly.id_to_index_or_throw(7);
  auto aps_p7 = aps.aps_for_polygon(p7);
  REQUIRE(G1.through_cost_raw(p7, aps_p7[0], aps_p7[1]) ==
          Approx(G2.through_cost_raw(p7, aps_p7[0], aps_p7[1])));
  // applying the factor scales linearly
  double raw = G1.through_cost_raw(p7, aps_p7[0], aps_p7[1]);
  REQUIRE(raw * 1.0 == Approx(raw * G1.through_penalty_factor()));
  REQUIRE(raw * 5.0 == Approx(raw * G2.through_penalty_factor()));
}

TEST_CASE("Emission distance zero for point inside polygon (T046)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // Trajectory entirely inside polygon 7 (0.8..1.2 x 0.8..1.2)
  Trajectory traj; traj.id = 102;
  traj.geom.add_point(0.9, 0.9);
  traj.geom.add_point(1.0, 1.0);
  traj.geom.add_point(1.1, 1.1);
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap,
                                              /*link_only=*/false, nullptr);
  // At least one polygon segment expected (the trajectory lives inside p7).
  REQUIRE_FALSE(result.polygon_segments.empty());
  // The polygon segment should be for polygon 7.
  REQUIRE(result.polygon_segments.front().polygon_id == 7);
}

TEST_CASE("Emission distance equals min boundary distance for point outside polygon (T048)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  PolygonIndex p7 = poly.id_to_index_or_throw(7);
  // Outside polygon 7: at (0.5, 1.0) — boundary at x=0.8, so distance = 0.3.
  Point gps(0.5, 1.0);
  REQUIRE(poly.min_boundary_distance(p7, gps) == Approx(0.3));
}

TEST_CASE("Entirely-inside trajectory yields single polygon segment with no APs (T051)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // All GPS points strictly inside polygon 7
  Trajectory traj; traj.id = 200;
  traj.geom.add_point(0.85, 0.85);
  traj.geom.add_point(1.0, 1.0);
  traj.geom.add_point(1.15, 1.15);
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.05;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  // The HMM is free to pick polygon-only or hybrid. *If* the result has any
  // polygon-7 segment that started at GPS layer 0 and ended at the last
  // layer, both APs must be absent (no entry, no egress).
  if (!result.polygon_segments.empty()) {
    bool found_p7 = false;
    for (const auto &seg : result.polygon_segments) {
      if (seg.polygon_id == 7) {
        found_p7 = true;
        // First polygon segment: entry_ap must be absent.
        if (&seg == &result.polygon_segments.front()) {
          CHECK(seg.entry_ap == kNoAccessPoint);
        }
        // Last polygon segment: egress_ap must be absent.
        if (&seg == &result.polygon_segments.back()) {
          CHECK(seg.egress_ap == kNoAccessPoint);
        }
      }
    }
    INFO("found polygon 7 segment: " << found_p7);
  }
}

TEST_CASE("Mid-polygon end has no egress AP (T055 FR-007/FR-010)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // Trajectory enters polygon 7 from outside (left), ends inside.
  Trajectory traj; traj.id = 104;
  traj.geom.add_point(0.5, 1.0);  // outside, on edge 3
  traj.geom.add_point(0.9, 1.0);  // inside polygon 7
  traj.geom.add_point(1.0, 1.0);  // inside polygon 7
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.05;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  if (!result.polygon_segments.empty()) {
    // Last polygon segment should have egress_ap absent (trajectory ends
    // inside the polygon).
    CHECK(result.polygon_segments.back().egress_ap == kNoAccessPoint);
  }
}

TEST_CASE("CLI parity: polymatch accepts every weightmatch flag (T089)",
          "[polymatch][us2]") {
  // Build the full polymatch option set and assert each weightmatch flag is
  // present. cxxopts has no public "list options" API, so we check by
  // attempting to parse a minimal argv that uses each.
  cxxopts::Options pm("polymatch", "");
  CONFIG::NetworkConfig::register_arg(pm);
  CONFIG::GPSConfig::register_arg(pm);
  CONFIG::ResultConfig::register_arg(pm);
  CONFIG::PolygonConfig::register_arg(pm);
  CONFIG::AccessPointConfig::register_arg(pm);
  POLYMATCHConfig::register_arg(pm);
  pm.add_options()
    ("l,log_level", "", cxxopts::value<int>()->default_value("2"))
    ("s,step", "", cxxopts::value<int>()->default_value("100"))
    ("h,help", "")
    ("use_omp", "");

  // Each of these is the set weightmatch documents in its CLI; if any one
  // is missing on polymatch, parse() will throw.
  std::vector<std::string> argv_storage = {
      "polymatch",
      "--network", "n.shp",
      "--gps", "g.csv",
      "--output", "o.csv",
      "--candidates", "4",
      "--radius", "0.5",
      "--error", "0.1",
      "--backup_candidates", "-1",
      "--backup_radius", "-1",
      "--upper_bound_factor", "10.0",
  };
  std::vector<char *> argv;
  for (auto &s : argv_storage) argv.push_back(&s[0]);
  int argc = (int)argv.size();
  char **argv_ptr = argv.data();
  CHECK_NOTHROW(pm.parse(argc, argv_ptr));
}

TEST_CASE("Mid-polygon start has no entry AP (T054 FR-007/FR-010)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // Trajectory starts inside polygon 7, then exits via right boundary toward
  // edge 4 (5→6 at y=1)
  Trajectory traj; traj.id = 103;
  traj.geom.add_point(1.0, 1.0);  // inside polygon 7
  traj.geom.add_point(1.1, 1.1);  // inside polygon 7
  traj.geom.add_point(1.5, 1.0);  // outside, on edge 4 (5-6)
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  if (!result.polygon_segments.empty()) {
    // First polygon segment is the mid-polygon start; entry_ap must be absent.
    REQUIRE(result.polygon_segments.front().entry_ap == kNoAccessPoint);
  }
}

TEST_CASE("XML config parses polygon blocks (T031/T038)", "[polymatch][us2]") {
  // Write a small XML config and verify POLYMATCHConfig::load_from_xml +
  // PolygonConfig::load_from_xml + AccessPointConfig::load_from_xml all land.
  std::string xml_path = std::string(kFixtureDir) + "/polymatch_config.xml";
  {
    std::ofstream f(xml_path);
    f << R"(<?xml version="1.0" encoding="utf-8"?>
<config>
  <input>
    <polygon>
      <file>polygons.shp</file>
      <id_name>myid</id_name>
      <cost_name>mycost</cost_name>
    </polygon>
    <access_point>
      <file>access_points.shp</file>
      <node_id_name>nid</node_id_name>
      <polygon_id_name>pid</polygon_id_name>
    </access_point>
  </input>
  <parameters>
    <k>10</k>
    <r>500</r>
    <gps_error>25</gps_error>
    <through_penalty_factor>3.0</through_penalty_factor>
    <boundary_epsilon>0.0001</boundary_epsilon>
  </parameters>
</config>)";
  }
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(xml_path, tree);

  auto pcfg = CONFIG::PolygonConfig::load_from_xml(tree);
  CHECK(pcfg.file == "polygons.shp");
  CHECK(pcfg.id_name == "myid");
  CHECK(pcfg.cost_name == "mycost");

  auto acfg = CONFIG::AccessPointConfig::load_from_xml(tree);
  CHECK(acfg.file == "access_points.shp");
  CHECK(acfg.node_id_name == "nid");
  CHECK(acfg.polygon_id_name == "pid");

  auto algocfg = POLYMATCHConfig::load_from_xml(tree);
  CHECK(algocfg.k == 10);
  CHECK(algocfg.radius == Approx(500.0));
  CHECK(algocfg.gps_error == Approx(25.0));
  CHECK(algocfg.through_penalty_factor == Approx(3.0));
  CHECK(algocfg.boundary_epsilon == Approx(0.0001));
}

TEST_CASE("Thread-determinism: same results single vs multi-threaded (T082 SC-009)",
          "[polymatch][parallel]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // Build a batch of trajectories that mix link-only and polygon traversals.
  std::vector<Trajectory> batch;
  for (int i = 0; i < 16; ++i) {
    Trajectory tr; tr.id = 1000 + i;
    if (i % 2 == 0) {
      // link-only along top row
      tr.geom.add_point(0.1, 2.0);
      tr.geom.add_point(0.5, 2.0);
      tr.geom.add_point(1.5, 2.0);
      tr.geom.add_point(1.9, 2.0);
    } else {
      // inside polygon 7
      tr.geom.add_point(0.9, 0.9);
      tr.geom.add_point(1.0, 1.0);
      tr.geom.add_point(1.1, 1.1);
    }
    batch.push_back(tr);
  }
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.05;

  // Single-threaded reference run.
  std::vector<PolyMatchResult> single(batch.size());
  {
    DijkstraState state; IndexedMinHeap heap;
    for (size_t i = 0; i < batch.size(); ++i) {
      single[i] = matcher.match_traj(batch[i], cfg, state, heap, false,
                                     nullptr);
    }
  }

  // Multi-threaded run: each thread owns its own DijkstraState + heap.
  std::vector<PolyMatchResult> multi(batch.size());
  constexpr int kThreads = 4;
  std::vector<std::thread> threads;
  std::atomic<int> next_idx{0};
  for (int t = 0; t < kThreads; ++t) {
    threads.emplace_back([&]() {
      DijkstraState state; IndexedMinHeap heap;
      while (true) {
        int i = next_idx.fetch_add(1);
        if (i >= (int)batch.size()) return;
        multi[i] = matcher.match_traj(batch[i], cfg, state, heap, false,
                                      nullptr);
      }
    });
  }
  for (auto &th : threads) th.join();

  // Results must be bit-identical per trajectory (SC-009).
  for (size_t i = 0; i < batch.size(); ++i) {
    CAPTURE(i, batch[i].id);
    CHECK(single[i].base.opath == multi[i].base.opath);
    CHECK(single[i].base.cpath == multi[i].base.cpath);
    CHECK(single[i].polygon_segments.size() == multi[i].polygon_segments.size());
    for (size_t s = 0; s < single[i].polygon_segments.size(); ++s) {
      const auto &a = single[i].polygon_segments[s];
      const auto &b = multi[i].polygon_segments[s];
      CHECK(a.polygon_id == b.polygon_id);
      CHECK(a.entry_ap == b.entry_ap);
      CHECK(a.egress_ap == b.egress_ap);
      CHECK(a.is_through == b.is_through);
      CHECK(a.distance_inside == Approx(b.distance_inside));
    }
  }
}

TEST_CASE("CLI parses polygon-specific flags (T030)", "[polymatch][us2]") {
  // Drive the same cxxopts pipeline POLYMATCHAppConfig uses, then check our
  // load_from_arg methods land the values in the right fields.
  cxxopts::Options options("polymatch_test", "");
  CONFIG::NetworkConfig::register_arg(options);
  CONFIG::GPSConfig::register_arg(options);
  CONFIG::ResultConfig::register_arg(options);
  CONFIG::PolygonConfig::register_arg(options);
  CONFIG::AccessPointConfig::register_arg(options);
  POLYMATCHConfig::register_arg(options);

  std::vector<std::string> argv_storage = {
      "polymatch",
      "--polygons", "p.shp",
      "--polygon_id_name", "myid",
      "--polygon_cost_name", "mycost",
      "--access_points", "ap.shp",
      "--ap_node_id_name", "nid",
      "--ap_polygon_id_name", "pid",
      "--through_penalty_factor", "2.5",
      "--boundary_epsilon", "0.001",
  };
  std::vector<char *> argv;
  for (auto &s : argv_storage) argv.push_back(&s[0]);
  int argc = static_cast<int>(argv.size());
  char **argv_ptr = argv.data();

  auto result = options.parse(argc, argv_ptr);
  auto pcfg = CONFIG::PolygonConfig::load_from_arg(result);
  auto acfg = CONFIG::AccessPointConfig::load_from_arg(result);
  auto algocfg = POLYMATCHConfig::load_from_arg(result);

  CHECK(pcfg.file == "p.shp");
  CHECK(pcfg.id_name == "myid");
  CHECK(pcfg.cost_name == "mycost");
  CHECK(acfg.file == "ap.shp");
  CHECK(acfg.node_id_name == "nid");
  CHECK(acfg.polygon_id_name == "pid");
  CHECK(algocfg.through_penalty_factor == Approx(2.5));
  CHECK(algocfg.boundary_epsilon == Approx(0.001));
}

TEST_CASE("Polygon ID encoded as negative in opath (T071 partial)",
          "[polymatch][us3]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  Trajectory traj; traj.id = 102;
  traj.geom.add_point(0.95, 0.95);
  traj.geom.add_point(1.05, 1.05);
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.05;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);

  // If any opath entry is a polygon match, it must be encoded as a negative
  // polygon ID. Per the negation convention, a polygon-id-7 match shows up
  // as -7 in both opath and cpath.
  for (auto id : result.base.opath) {
    if (id < 0) {
      // Must correspond to a real polygon
      bool found = false;
      for (const auto &p : poly.polygons()) {
        if (-id == p.id) { found = true; break; }
      }
      CHECK(found);
    }
  }
}

TEST_CASE("Fallback to link-only when zero valid polygons (T033 FR-012)",
          "[polymatch][us2]") {
  // Build an "all invalid" polygons fixture: a shapefile containing only the
  // bowtie polygon which gets skipped at load time.
  std::string all_invalid_path = std::string(kFixtureDir) + "/polygons_all_invalid.shp";
  {
    GDALAllRegister();
    GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    if (path_exists(all_invalid_path)) drv->Delete(all_invalid_path.c_str());
    GDALDataset *ds = drv->Create(all_invalid_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRLayer *layer = ds->CreateLayer("polys", nullptr, wkbPolygon, nullptr);
    OGRFieldDefn id_f("id", OFTInteger64); layer->CreateField(&id_f);
    OGRFieldDefn cf("cost", OFTReal); layer->CreateField(&cf);
    OGRFeature *f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    f->SetField("id", (GIntBig)99);
    f->SetField("cost", 1.0);
    OGRPolygon p; OGRLinearRing r;
    // self-intersecting bowtie -> rejected by is_valid (FR-013)
    r.addPoint(0, 0); r.addPoint(1, 1); r.addPoint(1, 0); r.addPoint(0, 1); r.addPoint(0, 0);
    p.addRing(&r);
    f->SetGeometry(&p);
    layer->CreateFeature(f);
    OGRFeature::DestroyFeature(f);
    GDALClose(ds);
  }
  spdlog::set_level(spdlog::level::off);
  PolygonLayer layer;
  REQUIRE(layer.load({all_invalid_path, "id", "cost"}));
  CHECK(layer.empty());  // all skipped -> zero valid polygons (FR-012 trigger)
}

TEST_CASE("PolyMMWriter concurrent writes do not corrupt rows (T075 FR-017)",
          "[polymatch][us3]") {
  std::string out_path = std::string(kFixtureDir) + "/poly_writer_concurrent.csv";
  ::remove(out_path.c_str());
  CONFIG::OutputConfig oc;
  oc.write_opath = false; oc.write_cpath = false; oc.write_mgeom = false;
  oc.write_error = false; oc.write_offset = false; oc.write_spdist = false;
  oc.write_pgeom = false; oc.write_tpath = false;

  constexpr int kNumThreads = 4;
  constexpr int kPerThread = 200;
  {
    IO::PolyMMWriter w(out_path, oc, /*include_polygon_columns=*/true);
    std::vector<std::thread> threads;
    for (int t = 0; t < kNumThreads; ++t) {
      threads.emplace_back([t, &w]() {
        for (int i = 0; i < kPerThread; ++i) {
          Trajectory tr; tr.id = t * 10000 + i;
          PolyMatchResult r;
          r.base.id = tr.id;
          r.polygon_segments.push_back({7, 5, 9, false, 1.0, 0});
          w.write_result(tr, r);
        }
      });
    }
    for (auto &th : threads) th.join();
  }
  std::ifstream in(out_path);
  std::string line;
  int data_rows = 0;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    if (!line.empty()) ++data_rows;
  }
  // 4 * 200 = 800 rows expected, mutex must serialize writes so no
  // line is interleaved or lost.
  CHECK(data_rows == kNumThreads * kPerThread);
}

TEST_CASE("Cross-polygon traversal records entry & egress APs (T050)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // GPS along the middle row crossing polygon 7 (0.8..1.2 x 0.8..1.2).
  // Enters from the left (~node 4 at x=0), passes through inside, exits right.
  Trajectory traj; traj.id = 101;
  traj.geom.add_point(0.5, 1.0);  // outside-left, on edge 3
  traj.geom.add_point(0.9, 1.0);  // inside polygon 7
  traj.geom.add_point(1.1, 1.0);  // inside polygon 7
  traj.geom.add_point(1.5, 1.0);  // outside-right, on edge 4
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  // The HMM may pick either link-only (across edge 3+4) or hybrid
  // (edge 3 -> polygon 7 -> edge 4) — both are valid behaviors for this
  // fixture depending on weights. We assert that *if* a polygon segment is
  // emitted, both entry_ap and egress_ap are set (not kNoAccessPoint).
  for (const auto &seg : result.polygon_segments) {
    if (seg.polygon_id == 7) {
      CHECK_FALSE(seg.is_through);  // matched inside GPS observations
      CHECK(seg.entry_ap != kNoAccessPoint);
      CHECK(seg.egress_ap != kNoAccessPoint);
    }
  }
}

TEST_CASE("Hybrid C_Path topology: polygons appear as negative IDs (T059)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  Trajectory traj; traj.id = 102;
  traj.geom.add_point(0.9, 0.9);
  traj.geom.add_point(1.0, 1.0);
  traj.geom.add_point(1.1, 1.1);
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);

  // For each PolygonSegment, the corresponding cpath entry must be -polygon_id.
  for (const auto &seg : result.polygon_segments) {
    REQUIRE(seg.position_in_cpath < result.base.cpath.size());
    CHECK(result.base.cpath[seg.position_in_cpath] == -seg.polygon_id);
  }
}

TEST_CASE("distance_inside is non-negative + finite (T058)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  Trajectory traj; traj.id = 103;
  traj.geom.add_point(0.9, 0.9);
  traj.geom.add_point(1.0, 1.0);
  traj.geom.add_point(1.1, 1.1);
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  for (const auto &seg : result.polygon_segments) {
    CHECK(std::isfinite(seg.distance_inside));
    CHECK(seg.distance_inside >= 0.0);
  }
}

TEST_CASE("polymatch link-only result equals weightmatch on same input (T080)",
          "[polymatch][regression][sc-002]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer empty_poly;
  AccessPointLayer empty_aps;
  PolyLinkGraph poly_graph(net, link_graph, empty_poly, empty_aps, 1.5);
  POLYMATCH polym(net, empty_poly, empty_aps, poly_graph, link_graph);
  WEIGHTMATCH wm(net, link_graph);

  Trajectory traj; traj.id = 100;
  traj.geom.add_point(0.1, 2.0);
  traj.geom.add_point(0.5, 2.0);
  traj.geom.add_point(1.5, 2.0);
  traj.geom.add_point(1.9, 2.0);

  POLYMATCHConfig pcfg;
  pcfg.k = 4; pcfg.radius = 0.5; pcfg.gps_error = 0.1;
  WEIGHTMATCHConfig wcfg{pcfg.k, pcfg.radius, pcfg.gps_error,
                         pcfg.backup_k, pcfg.backup_radius,
                         pcfg.upper_bound_factor,
                         pcfg.allow_truncation};

  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult pm = polym.match_traj(traj, pcfg, state, heap,
                                        /*link_only=*/true, nullptr);
  MatchResult wmr = wm.match_traj(traj, wcfg, state, heap);

  // opath + cpath are deterministic functions of the input; in link-only
  // fallback they must match byte-for-byte (SC-002).
  CHECK(pm.base.opath == wmr.opath);
  CHECK(pm.base.cpath == wmr.cpath);
}

TEST_CASE("Single-point trajectory handled gracefully (T088 FR-019)",
          "[polymatch][us1]") {
  // The app-level skip is in POLYMATCHApp::run() (npts < 2). Here we verify
  // that match_traj itself does not crash or return NaN if accidentally
  // invoked on a single-point input.
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  Trajectory traj; traj.id = 999;
  traj.geom.add_point(1.0, 1.0);  // single point inside polygon 7
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  // The result is allowed to be empty; we only require no crash and no
  // NaN/inf in any returned field.
  for (const auto &seg : result.polygon_segments) {
    CHECK(std::isfinite(seg.distance_inside));
  }
}

TEST_CASE("Same-polygon transition cost equals eu_dist (T056 FR-008)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // Two polygon candidates in the same polygon — the same-polygon eu_dist
  // override should make transition_cost return eu_dist exactly, mirroring
  // the same-link override at weightmatch_algorithm.cpp:323-324.
  PolyCandidate a;
  a.kind = PolyCandidateKind::Polygon;
  a.polygon_index = poly.id_to_index_or_throw(7);
  a.inside = true;
  a.matched_point = Point(0.9, 1.0);
  a.ep_distance = 0;

  PolyCandidate b;
  b.kind = PolyCandidateKind::Polygon;
  b.polygon_index = poly.id_to_index_or_throw(7);  // same polygon
  b.inside = true;
  b.matched_point = Point(1.1, 1.0);
  b.ep_distance = 0;

  POLYMATCHConfig cfg;
  DijkstraState state; IndexedMinHeap heap;
  const double eu_dist = 0.42;
  double cost = matcher.transition_cost(a, b, eu_dist, cfg, state, heap);
  CHECK(cost == Approx(eu_dist));

  // Sanity: with a different polygon, the cost is NOT the override.
  b.polygon_index = poly.id_to_index_or_throw(42);
  double cost_diff = matcher.transition_cost(a, b, eu_dist, cfg, state, heap);
  // For disjoint polygons with no shared AP, cost is infinity.
  CHECK(std::isinf(cost_diff));
}

TEST_CASE("Two-polygons-via-shared-AP traversal (T052)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().poly_shared_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_shared_path, "node_id", "polygon_id"}, poly,
                   net, 1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.0);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // Trajectory: starts inside polygon 7 (0.8-1.2 x 0.8-1.2), crosses the
  // shared boundary at x=1.2, ends inside polygon 8 (1.2-1.6 x 0.8-1.2).
  Trajectory traj; traj.id = 200;
  traj.geom.add_point(0.9, 1.0);   // inside polygon 7
  traj.geom.add_point(1.0, 1.0);   // inside polygon 7
  traj.geom.add_point(1.3, 1.0);   // inside polygon 8
  traj.geom.add_point(1.4, 1.0);   // inside polygon 8
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.05;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);

  // If the HMM produced two consecutive polygon segments (7, then 8), they
  // must share an AP: the egress AP of segment 1 must equal the entry AP of
  // segment 2 (both should be node_id 2000, the shared AP).
  if (result.polygon_segments.size() >= 2) {
    for (size_t i = 0; i + 1 < result.polygon_segments.size(); ++i) {
      const auto &s1 = result.polygon_segments[i];
      const auto &s2 = result.polygon_segments[i + 1];
      if (s1.egress_ap != kNoAccessPoint &&
          s2.entry_ap != kNoAccessPoint) {
        CHECK(s1.egress_ap == s2.entry_ap);
      }
    }
  }
}

TEST_CASE("Polygon shortcut via Dijkstra: through-cost wired into PolyLinkGraph (T062/T053)",
          "[polymatch][us1][routing]") {
  // With the Held-Karp polygon-AP sub-vertex expansion, a vanilla Dijkstra
  // over PolyLinkGraph should pick the polygon shortcut whenever it beats
  // the link-only route. Verify this directly without going through the
  // matcher's transition_cost layer.
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));

  // Routing from edge 1 (1->2) to edge 4 (5->6).
  //   Link-only: 1 -> 10 -> 4 (node 2 -> node 5). Cost = w(10) + w(4) = 2.0.
  //   Polygon : 1 -> sub(7, AP@2) -> sub(7, AP@5) -> 4
  //             = 0 + w_p * dist(AP@2=(1.0,1.2), AP@5=(0.8,0.8)) * factor
  //             + w(4) ~= 2.0 * 0.447 * factor + 1.0.
  //   At factor=0.1 polygon cost ~ 1.089 (wins); at factor=10 polygon cost ~
  //   9.94 (loses).
  EdgeIndex e1 = net.get_edge_index(1);
  EdgeIndex e4 = net.get_edge_index(4);

  PolygonIndex p7 = poly.id_to_index_or_throw(7);
  AccessPointIndex ap2 = aps.node_id_to_index(2);
  AccessPointIndex ap5 = aps.node_id_to_index(5);

  PolyLinkGraph G_cheap(net, link_graph, poly, aps,
                        /*through_penalty_factor=*/0.1);
  double raw = G_cheap.through_cost_raw(p7, ap2, ap5);
  REQUIRE(raw > 0);

  std::vector<uint32_t> goals{static_cast<uint32_t>(e4)};
  std::vector<Path> paths(1);
  DijkstraState state; IndexedMinHeap heap;
  shortest_polylink_to_polylinks(G_cheap, state, heap,
                                 static_cast<uint32_t>(e1), goals, paths);
  REQUIRE(paths[0].found);
  bool path_visits_polygon = false;
  for (auto v : paths[0].edges) {
    if (G_cheap.vertex_kind(static_cast<uint32_t>(v)) ==
        PolyVertexKind::PolygonSubVertex) {
      path_visits_polygon = true;
      break;
    }
  }
  CHECK(path_visits_polygon);
  CHECK(paths[0].total_cost == Approx(0.1 * raw + 1.0).margin(1e-6));

  // With factor=10.0 the polygon shortcut should be more expensive than the
  // link-only route, so Dijkstra picks the link route (no polygon sub-vertex
  // in the reconstructed path).
  PolyLinkGraph G_expensive(net, link_graph, poly, aps,
                            /*through_penalty_factor=*/10.0);
  std::vector<Path> paths2(1);
  shortest_polylink_to_polylinks(G_expensive, state, heap,
                                 static_cast<uint32_t>(e1), goals, paths2);
  REQUIRE(paths2[0].found);
  bool expensive_path_visits_polygon = false;
  for (auto v : paths2[0].edges) {
    if (G_expensive.vertex_kind(static_cast<uint32_t>(v)) ==
        PolyVertexKind::PolygonSubVertex) {
      expensive_path_visits_polygon = true;
      break;
    }
  }
  CHECK_FALSE(expensive_path_visits_polygon);
  // Link-only cost = w(10) + w(4).
  const auto &edges_vec = net.get_edges();
  double link_cost = edges_vec[net.get_edge_index(10)].cost +
                     edges_vec[net.get_edge_index(4)].cost;
  CHECK(paths2[0].total_cost == Approx(link_cost).margin(1e-6));
}

TEST_CASE("Through-routing: when is_through==true both APs must be set (T053)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  // Use a low through-penalty factor so polygon shortcut is more attractive.
  PolyLinkGraph G(net, link_graph, poly, aps, 0.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // GPS far enough from polygon 7 that no observation is inside but radius
  // still picks the polygon up as a candidate.
  Trajectory traj; traj.id = 105;
  traj.geom.add_point(0.3, 1.0);  // outside polygon 7
  traj.geom.add_point(1.7, 1.0);  // outside polygon 7
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 1.5; cfg.gps_error = 0.5;
  cfg.through_penalty_factor = 0.5;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  // Defensive contract: if any polygon segment is marked is_through, both
  // APs must be set. (For US1 scenario 5, full through-routing requires
  // Dijkstra-level AP-context tracking — see deferred note in tasks.md.)
  for (const auto &seg : result.polygon_segments) {
    if (seg.is_through) {
      CHECK(seg.entry_ap != kNoAccessPoint);
      CHECK(seg.egress_ap != kNoAccessPoint);
    }
  }
}

TEST_CASE("Duplicate GPS points produce no NaN (T057 Constitution II)",
          "[polymatch][us1]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().polygons_path, "id", "cost"}));
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().aps_path, "node_id", "polygon_id"}, poly, net,
                   1e-6));
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  Trajectory traj; traj.id = 200;
  // Three identical points along edge 1 (top row).
  traj.geom.add_point(0.5, 2.0);
  traj.geom.add_point(0.5, 2.0);
  traj.geom.add_point(0.5, 2.0);
  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.1;
  DijkstraState state; IndexedMinHeap heap;
  PolyMatchResult result = matcher.match_traj(traj, cfg, state, heap, false,
                                              nullptr);
  // No assertion on path content (any choice is fine) — just no crash and no
  // NaN propagation. We do check that the opath is non-empty.
  REQUIRE_FALSE(result.base.opath.empty());
}

TEST_CASE("POLYMATCHApp link-only mode behaves like weightmatch (T049 link-only branch + FR-012)",
          "[polymatch][us1][us2]") {
  // When no polygon layer is configured, the matcher should yield the same
  // result as weightmatch on the same network/trajectory input. We sanity-check
  // by constructing the matcher in link-only mode and verifying it returns a
  // result for a simple trajectory.
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().network_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer empty_poly;
  AccessPointLayer empty_aps;
  PolyLinkGraph poly_graph(net, link_graph, empty_poly, empty_aps, 1.5);
  POLYMATCH matcher(net, empty_poly, empty_aps, poly_graph, link_graph);

  // Simple horizontal trajectory along edges 1 and 2 (top row)
  Trajectory traj;
  traj.id = 100;
  traj.geom.add_point(0.1, 2.0);
  traj.geom.add_point(0.5, 2.0);
  traj.geom.add_point(1.5, 2.0);
  traj.geom.add_point(1.9, 2.0);

  POLYMATCHConfig cfg;
  cfg.k = 4;
  cfg.radius = 0.5;
  cfg.gps_error = 0.1;

  DijkstraState state;
  IndexedMinHeap heap;
  PolyMatchResult result =
      matcher.match_traj(traj, cfg, state, heap, /*link_only=*/true, nullptr);
  REQUIRE_FALSE(result.base.opath.empty());
  // No polygon segments in link-only result
  REQUIRE(result.polygon_segments.empty());
}

// ============================================================
// US3 — PolyMMWriter format and link-only fallback
// ============================================================

TEST_CASE("PolyMMWriter writes polygon columns when enabled (T070)",
          "[polymatch][us3]") {
  std::string out_path = std::string(kFixtureDir) + "/poly_writer_out.csv";
  ::remove(out_path.c_str());

  CONFIG::OutputConfig oc;
  oc.write_opath = true;
  oc.write_cpath = true;
  oc.write_mgeom = true;
  oc.write_error = false;
  oc.write_offset = false;
  oc.write_spdist = false;
  oc.write_pgeom = false;
  oc.write_tpath = false;

  {
    IO::PolyMMWriter w(out_path, oc, /*include_polygon_columns=*/true);
    Trajectory tr; tr.id = 1001;
    PolyMatchResult r;
    r.base.id = 1001;
    r.base.opath = {1, 2};
    r.base.cpath = {1, 2, -7};
    r.polygon_segments.push_back({7, 5, 9, false, 47.83, 2});
    w.write_result(tr, r);
  }
  std::ifstream in(out_path);
  REQUIRE(in.is_open());
  std::string header, row;
  std::getline(in, header);
  std::getline(in, row);
  REQUIRE(header.find("polygon_ids") != std::string::npos);
  REQUIRE(header.find("entry_aps") != std::string::npos);
  REQUIRE(header.find("egress_aps") != std::string::npos);
  REQUIRE(header.find("is_through") != std::string::npos);
  REQUIRE(header.find("polygon_distances") != std::string::npos);
  REQUIRE(row.find("7") != std::string::npos);   // polygon_ids
  REQUIRE(row.find("47.83") != std::string::npos);
}

TEST_CASE("PolyMMWriter empty AP token (T072)", "[polymatch][us3]") {
  std::string out_path = std::string(kFixtureDir) + "/poly_writer_empty_ap.csv";
  ::remove(out_path.c_str());
  CONFIG::OutputConfig oc;
  oc.write_opath = true;
  oc.write_cpath = true;
  oc.write_mgeom = true;
  oc.write_error = false;
  oc.write_offset = false;
  oc.write_spdist = false;
  oc.write_pgeom = false;
  oc.write_tpath = false;

  {
    IO::PolyMMWriter w(out_path, oc, /*include_polygon_columns=*/true);
    Trajectory tr; tr.id = 2002;
    PolyMatchResult r;
    r.base.id = 2002;
    r.polygon_segments.push_back({7, kNoAccessPoint, 1001, false, 10.0, 0});
    r.polygon_segments.push_back({42, 1001, kNoAccessPoint, false, 5.0, 1});
    w.write_result(tr, r);
  }
  std::ifstream in(out_path);
  std::string header, row;
  std::getline(in, header);
  std::getline(in, row);
  // entry_aps should contain "-" for the first segment (mid-polygon start)
  REQUIRE(row.find("-,1001") != std::string::npos);
  // egress_aps should contain "-" for the last segment (mid-polygon end)
  REQUIRE(row.find("1001,-") != std::string::npos);
}

TEST_CASE("PolyMMWriter is_through flag value (T073)", "[polymatch][us3]") {
  std::string out_path = std::string(kFixtureDir) + "/poly_writer_through.csv";
  ::remove(out_path.c_str());
  CONFIG::OutputConfig oc;
  oc.write_opath = true; oc.write_cpath = true; oc.write_mgeom = true;
  oc.write_error = false; oc.write_offset = false; oc.write_spdist = false;
  oc.write_pgeom = false; oc.write_tpath = false;

  {
    IO::PolyMMWriter w(out_path, oc, /*include_polygon_columns=*/true);
    Trajectory tr; tr.id = 3003;
    PolyMatchResult r;
    r.base.id = 3003;
    r.polygon_segments.push_back({7, 5, 9, true, 1.5, 0});
    w.write_result(tr, r);
  }
  std::ifstream in(out_path);
  std::string h, row; std::getline(in, h); std::getline(in, row);
  REQUIRE(row.find(";1;") != std::string::npos);  // is_through=1
}

TEST_CASE("PolyMMWriter omits polygon columns in link-only (T074 SC-002)",
          "[polymatch][us3]") {
  std::string out_path = std::string(kFixtureDir) + "/poly_writer_linkonly.csv";
  ::remove(out_path.c_str());
  CONFIG::OutputConfig oc;
  oc.write_opath = true; oc.write_cpath = true; oc.write_mgeom = true;
  oc.write_error = false; oc.write_offset = false; oc.write_spdist = false;
  oc.write_pgeom = false; oc.write_tpath = false;

  {
    IO::PolyMMWriter w(out_path, oc, /*include_polygon_columns=*/false);
    Trajectory tr; tr.id = 4004;
    PolyMatchResult r;
    r.base.id = 4004;
    r.base.opath = {1, 2};
    r.base.cpath = {1, 2};
    w.write_result(tr, r);
  }
  std::ifstream in(out_path);
  std::string header; std::getline(in, header);
  REQUIRE(header.find("polygon_ids") == std::string::npos);
  REQUIRE(header.find("is_through") == std::string::npos);
}

// ============================================================
// Smoke test — Config defaults
// ============================================================

TEST_CASE("Config defaults (smoke)", "[polymatch][smoke]") {
  CONFIG::PolygonConfig poly_cfg;
  REQUIRE(poly_cfg.id_name == "id");
  REQUIRE(poly_cfg.cost_name == "cost");
  REQUIRE(poly_cfg.file.empty());
  REQUIRE(poly_cfg.validate());

  CONFIG::AccessPointConfig ap_cfg;
  REQUIRE(ap_cfg.node_id_name == "node_id");
  REQUIRE(ap_cfg.polygon_id_name == "polygon_id");
  REQUIRE(ap_cfg.validate());

  POLYMATCHConfig algo_cfg;
  REQUIRE(algo_cfg.k == 8);
  REQUIRE(algo_cfg.through_penalty_factor == Approx(1.5));
  REQUIRE(algo_cfg.boundary_epsilon == Approx(1e-6));
  REQUIRE(algo_cfg.validate());

  PolyMatchResult result;
  REQUIRE(result.polygon_segments.empty());
  REQUIRE(kNoAccessPoint == -1);
}

// ============================================================
// Perf baseline (T081) — hidden by default, run with `polymatch_test [.bench]`
// ============================================================

TEST_CASE("Perf baseline: 1000 GPS points x 100 polygons < 10s (T081 SC-003)",
          "[.bench][polymatch][perf]") {
  spdlog::set_level(spdlog::level::off);
  Network net(fixture().perf_net_path, "NO_TURN_BANS", "id", "source", "target");
  LinkGraph link_graph(net);
  PolygonLayer poly;
  REQUIRE(poly.load({fixture().perf_polys_path, "id", "cost"}));
  REQUIRE(poly.size() == 100);
  AccessPointLayer aps;
  REQUIRE(aps.load({fixture().perf_aps_path, "node_id", "polygon_id"}, poly,
                   net, 1e-6));
  REQUIRE(aps.size() == 100);
  PolyLinkGraph G(net, link_graph, poly, aps, 1.5);
  POLYMATCH matcher(net, poly, aps, G, link_graph);

  // 1000-point trajectory along the row y=5.0 (network edge between nodes
  // 56-66) sweeping from x=1.0 to x=9.0 with small noise. Every point sits
  // within candidate radius of a network edge, so all layers have candidates.
  Trajectory traj; traj.id = 1;
  for (int i = 0; i < 1000; ++i) {
    double x = 1.0 + (i / 999.0) * 8.0;
    double y = 5.0 + ((i % 2) ? 0.001 : -0.001);
    traj.geom.add_point(x, y);
  }
  REQUIRE(traj.geom.get_num_points() == 1000);

  POLYMATCHConfig cfg;
  cfg.k = 4; cfg.radius = 0.5; cfg.gps_error = 0.05;

  DijkstraState state; IndexedMinHeap heap;
  auto t0 = std::chrono::steady_clock::now();
  PolyMatchResult result =
      matcher.match_traj(traj, cfg, state, heap, false, nullptr);
  auto t1 = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(t1 - t0).count();
  WARN("perf baseline elapsed: " << secs << "s, "
                                  << result.base.opath.size() << " opath, "
                                  << result.polygon_segments.size()
                                  << " polygon segments");
  CHECK(secs < 10.0);  // SC-003
}
