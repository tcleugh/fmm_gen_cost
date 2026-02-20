#!/usr/bin/env python3
"""Generate synthetic GPS trips on an ESRI shapefile road network.

Produces CSV files in FMM's CSVTrajectoryReader format (id;geom with WKT LINESTRING).
"""

import argparse
import math
import random
import sys
from collections import defaultdict

try:
    import geopandas as gpd
except ImportError:
    sys.exit("geopandas required. Install with: pip install geopandas")


def load_network(shp_path, id_field="id", source_field="source", target_field="target"):
    """Load edges from shapefile, return edges list and adjacency map."""
    gdf = gpd.read_file(shp_path)

    edges = []  # list of (edge_id, source, target, [(x,y), ...])
    adjacency = defaultdict(list)  # node -> [edge_index, ...]

    for _, row in gdf.iterrows():
        eid = int(row[id_field])
        src = int(row[source_field])
        tgt = int(row[target_field])
        coords = list(row.geometry.coords)
        edge_idx = len(edges)
        edges.append((eid, src, tgt, coords))
        adjacency[src].append(edge_idx)

    return edges, adjacency


def edge_length(coords):
    """Euclidean length of a polyline."""
    total = 0.0
    for i in range(len(coords) - 1):
        dx = coords[i + 1][0] - coords[i][0]
        dy = coords[i + 1][1] - coords[i][1]
        total += math.sqrt(dx * dx + dy * dy)
    return total


def interpolate_along_edge(coords, fraction):
    """Get (x, y) at a given fraction [0, 1] along the polyline."""
    total = edge_length(coords)
    if total < 1e-12:
        return coords[0]
    target_dist = fraction * total
    accum = 0.0
    for i in range(len(coords) - 1):
        dx = coords[i + 1][0] - coords[i][0]
        dy = coords[i + 1][1] - coords[i][1]
        seg_len = math.sqrt(dx * dx + dy * dy)
        if accum + seg_len >= target_dist - 1e-12:
            if seg_len < 1e-12:
                return coords[i]
            t = (target_dist - accum) / seg_len
            t = max(0.0, min(1.0, t))
            return (coords[i][0] + t * dx, coords[i][1] + t * dy)
        accum += seg_len
    return coords[-1]


def random_walk(edges, adjacency, num_edges, rng):
    """Perform a random walk through the network, returning a list of edge indices."""
    if not edges:
        return []
    start_idx = rng.randint(0, len(edges) - 1)
    path = [start_idx]
    for _ in range(num_edges - 1):
        _, _, tgt, _ = edges[path[-1]]
        neighbors = adjacency.get(tgt, [])
        if not neighbors:
            break
        path.append(rng.choice(neighbors))
    return path


def sample_points_along_path(edges, edge_indices, num_points, noise_sigma, rng):
    """Sample GPS points along a walked path with Gaussian noise."""
    # Build a continuous polyline from the path
    segments = []  # (edge_idx, coords, cumulative_start_dist)
    total_length = 0.0
    for ei in edge_indices:
        coords = edges[ei][3]
        seg_len = edge_length(coords)
        segments.append((ei, coords, total_length))
        total_length += seg_len

    if total_length < 1e-12 or num_points < 1:
        # Degenerate: just return midpoint of first edge
        coords = edges[edge_indices[0]][3]
        pt = interpolate_along_edge(coords, 0.5)
        return [pt]

    points = []
    for i in range(num_points):
        # Evenly space points along the total path length
        if num_points == 1:
            target_dist = total_length * 0.5
        else:
            target_dist = (i / (num_points - 1)) * total_length

        # Find which segment this falls in
        seg_idx = len(segments) - 1
        for j in range(len(segments) - 1):
            if segments[j + 1][2] > target_dist:
                seg_idx = j
                break

        _, coords, seg_start = segments[seg_idx]
        seg_len = edge_length(coords)
        if seg_len < 1e-12:
            frac = 0.5
        else:
            frac = (target_dist - seg_start) / seg_len
            frac = max(0.0, min(1.0, frac))

        px, py = interpolate_along_edge(coords, frac)

        # Add noise
        if noise_sigma > 0:
            px += rng.gauss(0, noise_sigma)
            py += rng.gauss(0, noise_sigma)

        points.append((px, py))

    return points


def format_linestring(points):
    """Format points as WKT LINESTRING."""
    if len(points) == 1:
        # Single point: duplicate it to form a valid linestring
        points = [points[0], points[0]]
    coords_str = ",".join(f"{x} {y}" for x, y in points)
    return f"LINESTRING({coords_str})"


def generate_basic_trips(edges, adjacency, args, rng):
    """Generate basic synthetic trips."""
    trips = []
    for trip_id in range(1, args.num_trips + 1):
        num_walk_edges = rng.randint(2, 8)
        path = random_walk(edges, adjacency, num_walk_edges, rng)
        if not path:
            continue
        num_pts = rng.randint(args.min_points, args.max_points)
        points = sample_points_along_path(edges, path, num_pts, args.noise, rng)
        trips.append((trip_id, points))
    return trips


def generate_edge_case_trips(edges, adjacency, args, rng):
    """Generate edge-case trips for corner-case testing."""
    trips = []
    tid = 1

    # 1. Single-point trajectory (3 of these)
    for _ in range(3):
        ei = rng.randint(0, len(edges) - 1)
        pt = interpolate_along_edge(edges[ei][3], rng.random())
        px = pt[0] + rng.gauss(0, args.noise)
        py = pt[1] + rng.gauss(0, args.noise)
        trips.append((tid, [(px, py)]))
        tid += 1

    # 2. Two-point trajectory (3 of these)
    for _ in range(3):
        path = random_walk(edges, adjacency, 2, rng)
        if path:
            points = sample_points_along_path(edges, path, 2, args.noise, rng)
            trips.append((tid, points))
            tid += 1

    # 3. Same-edge trips: all points near one edge (3 of these)
    for _ in range(3):
        ei = rng.randint(0, len(edges) - 1)
        num_pts = rng.randint(3, 6)
        points = sample_points_along_path(edges, [ei], num_pts, args.noise * 0.5, rng)
        trips.append((tid, points))
        tid += 1

    # 4. Long trips traversing many edges (3 of these)
    for _ in range(3):
        num_walk = rng.randint(10, 15)
        path = random_walk(edges, adjacency, num_walk, rng)
        if path:
            num_pts = rng.randint(15, 30)
            points = sample_points_along_path(edges, path, num_pts, args.noise, rng)
            trips.append((tid, points))
            tid += 1

    # 5. Trips with large gap (skip intermediate points) (3 of these)
    for _ in range(3):
        path = random_walk(edges, adjacency, 6, rng)
        if path and len(path) >= 4:
            all_pts = sample_points_along_path(edges, path, 10, args.noise, rng)
            # Keep first 2 and last 2 points, creating a gap
            gapped = all_pts[:2] + all_pts[-2:]
            trips.append((tid, gapped))
            tid += 1

    # 6. Points far from network (should fail to match) (3 of these)
    for _ in range(3):
        far_points = [
            (100 + rng.random(), 100 + rng.random()),
            (101 + rng.random(), 101 + rng.random()),
            (102 + rng.random(), 102 + rng.random()),
        ]
        trips.append((tid, far_points))
        tid += 1

    # 7. Normal trips with varying noise levels (2 of these)
    for noise_mult in [0.01, 0.3]:
        path = random_walk(edges, adjacency, 4, rng)
        if path:
            num_pts = rng.randint(5, 10)
            points = sample_points_along_path(
                edges, path, num_pts, args.noise * noise_mult, rng
            )
            trips.append((tid, points))
            tid += 1

    return trips


def write_trips(trips, output_path):
    """Write trips to CSV in FMM trajectory format."""
    with open(output_path, "w") as f:
        f.write("id;geom\n")
        for trip_id, points in trips:
            wkt = format_linestring(points)
            f.write(f"{trip_id};{wkt}\n")
    print(f"Wrote {len(trips)} trips to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic GPS trips")
    parser.add_argument("--network", required=True, help="Path to network shapefile")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument(
        "--mode",
        choices=["basic", "edge_cases"],
        default="basic",
        help="Generation mode",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--noise", type=float, default=0.05, help="GPS noise sigma (map units)"
    )
    parser.add_argument("--min_points", type=int, default=5, help="Min points per trip")
    parser.add_argument(
        "--max_points", type=int, default=20, help="Max points per trip"
    )
    parser.add_argument("--num_trips", type=int, default=50, help="Number of trips")
    parser.add_argument("--id_field", default="id", help="Edge ID field name")
    parser.add_argument("--source_field", default="source", help="Source node field")
    parser.add_argument("--target_field", default="target", help="Target node field")

    args = parser.parse_args()
    rng = random.Random(args.seed)

    edges, adjacency = load_network(
        args.network, args.id_field, args.source_field, args.target_field
    )
    print(f"Loaded network: {len(edges)} edges, {len(adjacency)} nodes")

    if args.mode == "basic":
        trips = generate_basic_trips(edges, adjacency, args, rng)
    else:
        trips = generate_edge_case_trips(edges, adjacency, args, rng)

    write_trips(trips, args.output)


if __name__ == "__main__":
    main()
