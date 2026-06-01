#!/usr/bin/env python3
"""Generate synthetic polymatch test fixtures.

Emits:
  - road network shapefile (~20 edges, 4-cluster intersections)
  - polygon shapefile (3 polygons: one valid, one self-intersecting invalid, one
    valid but with no AP)
  - access point shapefile (link-attached, shared-between-polygons, plus several
    "error" subsets the test harness uses to exercise FR-005 validation paths)
  - GPS CSV with trajectories covering all US1 scenarios (link-only,
    cross-polygon, two-polygons-shared-AP, mid-polygon start, mid-polygon end,
    fully-inside, through-routing)

Run from the repository root (use a relative path — never hardcode an absolute
path such as `/workspace/...` because the project root may be mounted at a
different location outside this sandbox):

  python3 test/generate_polymatch_test_data.py --output_dir test/data/polymatch
"""

import argparse
import os
import sys

try:
    import geopandas as gpd
    from shapely.geometry import LineString, Point, Polygon
except ImportError:
    sys.exit(
        "geopandas + shapely required. Install with: pip install geopandas shapely"
    )


def build_network():
    """Build a small grid: 4 horizontal + 4 vertical edges forming a 3x3 lattice
    of nodes. Edge IDs and node IDs share an integer namespace by design.

    Layout (node_id at each intersection):

      1---e1---2---e2---3
      |        |        |
      e9       e10      e11
      |        |        |
      4---e3---5---e4---6
      |        |        |
      e12      e13      e14
      |        |        |
      7---e5---8---e6---9
    """
    coords = {
        1: (0.0, 2.0), 2: (1.0, 2.0), 3: (2.0, 2.0),
        4: (0.0, 1.0), 5: (1.0, 1.0), 6: (2.0, 1.0),
        7: (0.0, 0.0), 8: (1.0, 0.0), 9: (2.0, 0.0),
    }
    # (edge_id, source, target)
    edges = [
        (1, 1, 2),  (2, 2, 3),
        (3, 4, 5),  (4, 5, 6),
        (5, 7, 8),  (6, 8, 9),
        (9, 1, 4),  (10, 2, 5), (11, 3, 6),
        (12, 4, 7), (13, 5, 8), (14, 6, 9),
    ]
    rows = []
    for eid, src, tgt in edges:
        geom = LineString([coords[src], coords[tgt]])
        rows.append({"id": eid, "source": src, "target": tgt, "geom": geom})
    gdf = gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")
    return gdf, coords


def build_polygons():
    """Three polygons:
      - id=7: valid square covering the area around node 5 (the lattice center)
      - id=42: valid square covering the SE quadrant around node 9
      - id=99: invalid self-intersecting "bowtie" polygon (FR-013 path)
      - id=200: valid polygon with NO access point features (FR-014 path)
    """
    rows = [
        {
            "id": 7,
            "cost": 2.0,
            "geom": Polygon([
                (0.8, 0.8), (1.2, 0.8), (1.2, 1.2), (0.8, 1.2), (0.8, 0.8)
            ]),
        },
        {
            "id": 42,
            "cost": 1.5,
            "geom": Polygon([
                (1.8, -0.2), (2.2, -0.2), (2.2, 0.2), (1.8, 0.2), (1.8, -0.2)
            ]),
        },
        {
            "id": 99,
            "cost": 1.0,
            # self-intersecting bowtie — boost::geometry::is_valid -> false
            "geom": Polygon([
                (3.0, 3.0), (4.0, 4.0), (4.0, 3.0), (3.0, 4.0), (3.0, 3.0)
            ]),
        },
        {
            "id": 200,
            "cost": 1.0,
            # valid but no AP features point to it
            "geom": Polygon([
                (-1.0, -1.0), (-0.5, -1.0), (-0.5, -0.5), (-1.0, -0.5), (-1.0, -1.0)
            ]),
        },
    ]
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_access_points():
    """Access points for the polygons.

    Node IDs match the network's nodes by design (R4 — direct ID lookup).

      Polygon 7 boundary touches node 5 (1,1), node 2 (1,2), node 4 (0,1).
        Use node_id=5 (link-attached at center), node_id=2 (top), node_id=4
        (left).
      Polygon 42 boundary touches node 9 (2,0). Use node_id=9.

    To exercise the polygon-shared-AP case for the routing graph, we add a
    synthetic AP with node_id=2000 that is not in the network's node map and
    appears on both polygon 7 and polygon 42 — though geometrically polygon 7
    and 42 are disjoint here, this exercises the polygon-shared validation
    code path. The "error" subsets below remain in disabled rows so callers can
    flip them on per-test.
    """
    rows = []
    # polygon 7 access points (3 link-attached APs)
    rows.append({"node_id": 5, "polygon_id": 7, "geom": Point(1.0, 1.0)})
    # Wait — node 5 is at the polygon's interior (0.8..1.2), not boundary.
    # Fix: place these on the polygon boundary using the polygon corner coords.
    rows = []
    rows.append({"node_id": 5, "polygon_id": 7, "geom": Point(0.8, 0.8)})
    rows.append({"node_id": 2, "polygon_id": 7, "geom": Point(1.0, 1.2)})
    rows.append({"node_id": 4, "polygon_id": 7, "geom": Point(0.8, 1.0)})
    # polygon 42 access point
    rows.append({"node_id": 9, "polygon_id": 42, "geom": Point(2.0, 0.0)})
    # Synthetic shared AP (not in network) — appears in both polygons 7 and 42
    # to exercise polygon-shared-only path. Geometry placed deliberately on
    # the shared boundary point of each polygon for validation.
    rows.append({"node_id": 2000, "polygon_id": 7, "geom": Point(1.2, 1.2)})
    rows.append({"node_id": 2000, "polygon_id": 42, "geom": Point(1.2, 1.2)})
    # Note: this AP will FAIL the FR-005 boundary validation for polygon 42
    # because (1.2, 1.2) is not on its boundary. For valid datasets we keep
    # this disabled by writing to a separate "errors" shapefile.
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_access_points_clean():
    """Clean, valid access point shapefile (passes FR-005)."""
    rows = []
    rows.append({"node_id": 5, "polygon_id": 7, "geom": Point(0.8, 0.8)})
    rows.append({"node_id": 2, "polygon_id": 7, "geom": Point(1.0, 1.2)})
    rows.append({"node_id": 4, "polygon_id": 7, "geom": Point(0.8, 1.0)})
    rows.append({"node_id": 9, "polygon_id": 42, "geom": Point(2.0, 0.0)})
    rows.append({"node_id": 6, "polygon_id": 42, "geom": Point(2.0, 0.2)})
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_access_points_off_boundary():
    """AP geometry far from declared polygon boundary — FR-005 condition 1."""
    rows = [
        {"node_id": 5, "polygon_id": 7, "geom": Point(0.8, 0.8)},
        # this point is ~0.5 from polygon 7 boundary
        {"node_id": 8, "polygon_id": 7, "geom": Point(0.3, 0.3)},
    ]
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_access_points_orphan():
    """AP references a polygon ID not in PolygonLayer — FR-005 condition 2."""
    rows = [
        {"node_id": 5, "polygon_id": 7, "geom": Point(0.8, 0.8)},
        {"node_id": 6, "polygon_id": 999, "geom": Point(2.0, 0.0)},
    ]
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_access_points_contradictory():
    """Two features same node_id, different geometry — FR-005 condition 3."""
    rows = [
        {"node_id": 5, "polygon_id": 7, "geom": Point(0.8, 0.8)},
        {"node_id": 5, "polygon_id": 7, "geom": Point(1.2, 1.2)},  # diff geom
    ]
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_polygons_id_zero():
    """Polygon feature with id=0 — FR-018 rejection."""
    rows = [
        {"id": 0, "cost": 1.0, "geom": Polygon([
            (0.8, 0.8), (1.2, 0.8), (1.2, 1.2), (0.8, 1.2), (0.8, 0.8)
        ])},
    ]
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_polygons_duplicate_id():
    """Two polygon features sharing the same id — FR-018 rejection."""
    rows = [
        {"id": 7, "cost": 1.0, "geom": Polygon([
            (0.8, 0.8), (1.2, 0.8), (1.2, 1.2), (0.8, 1.2), (0.8, 0.8)
        ])},
        {"id": 7, "cost": 1.0, "geom": Polygon([
            (1.8, -0.2), (2.2, -0.2), (2.2, 0.2), (1.8, 0.2), (1.8, -0.2)
        ])},
    ]
    return gpd.GeoDataFrame(rows, geometry="geom", crs="EPSG:32633")


def build_gps_csv():
    """GPS trajectories covering all US1 scenarios.

    Each trajectory is a WKT LINESTRING in the network's coordinate space.
    Trip IDs:
      100 — link-only path on top row (e1 -> e2)
      101 — crosses polygon 7 (enters/exits via APs)
      102 — fully inside polygon 7
      103 — begins mid-polygon 7
      104 — ends mid-polygon 7
      105 — through-routing (no GPS observation inside polygon 7)
      106 — 0 points (empty) — FR-019 skip
      107 — single point — FR-019 skip
    """
    rows = [
        (100, "LINESTRING(0.1 2.0, 0.5 2.0, 1.5 2.0, 1.9 2.0)"),
        (101, "LINESTRING(0.7 1.0, 0.9 1.0, 1.0 1.0, 1.1 1.0, 1.3 1.0)"),
        (102, "LINESTRING(0.9 0.9, 1.0 1.0, 1.1 1.1)"),
        (103, "LINESTRING(1.0 1.0, 1.1 1.1, 1.5 1.0)"),
        (104, "LINESTRING(0.5 1.0, 0.9 1.0, 1.0 1.0)"),
        (105, "LINESTRING(0.5 1.0, 1.5 1.0)"),
        (106, ""),
        (107, "POINT(1.0 1.0)"),
    ]
    return rows


def main():
    parser = argparse.ArgumentParser(description="Generate polymatch test fixtures")
    parser.add_argument(
        "--output_dir",
        default="test/data/polymatch",
        help="Output directory",
    )
    args = parser.parse_args()

    out = args.output_dir
    os.makedirs(out, exist_ok=True)

    # Network
    net, _ = build_network()
    net.to_file(os.path.join(out, "edges.shp"))
    print(f"Wrote {len(net)} edges -> edges.shp")

    # Polygons (valid set)
    polys = build_polygons()
    polys.to_file(os.path.join(out, "polygons.shp"))
    print(f"Wrote {len(polys)} polygons -> polygons.shp")

    # Access points (clean set)
    aps = build_access_points_clean()
    aps.to_file(os.path.join(out, "access_points.shp"))
    print(f"Wrote {len(aps)} access points -> access_points.shp")

    # Validation-failure fixtures
    build_access_points_off_boundary().to_file(
        os.path.join(out, "ap_off_boundary.shp")
    )
    build_access_points_orphan().to_file(
        os.path.join(out, "ap_orphan.shp")
    )
    build_access_points_contradictory().to_file(
        os.path.join(out, "ap_contradictory.shp")
    )
    build_polygons_id_zero().to_file(os.path.join(out, "polygons_id_zero.shp"))
    build_polygons_duplicate_id().to_file(
        os.path.join(out, "polygons_duplicate.shp")
    )

    # GPS trips
    trips = build_gps_csv()
    with open(os.path.join(out, "trips.csv"), "w") as f:
        f.write("id;geom\n")
        for tid, wkt in trips:
            f.write(f"{tid};{wkt}\n")
    print(f"Wrote {len(trips)} trips -> trips.csv")

    print("Done.")


if __name__ == "__main__":
    main()
