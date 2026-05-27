# Quickstart: PolyMatch

This guide walks through building polymatch and running it on a small example.

## Prerequisites

Same as the rest of the project: C++11 compiler, CMake ≥ 3.5, GDAL ≥ 2.2, Boost (Graph, Geometry,
Serialization) ≥ 1.56, OpenMP. See the project root `README.md` for installation instructions.

## Build

```bash
cd /workspace
mkdir -p build && cd build
cmake ..
make -j$(nproc) polymatch
```

The `polymatch` executable will be placed at `build/polymatch` alongside `weightmatch`, `stmatch`,
`fmm`, `h3mm`, and `ubodt_gen`.

## Build tests

```bash
cd build
make tests   # builds all test executables, including polymatch_test
```

All five test executables must pass:

```bash
./test/algorithm_test
./test/fmm_test
./test/network_test
./test/network_graph_test
./test/weightmatch_test
./test/polymatch_test          # NEW
```

## Minimal example: link-only fallback

To confirm polymatch is functionally equivalent to weightmatch when no polygon layer is provided
(spec SC-002):

```bash
./polymatch \
  --network ../test/data/edges.shp \
  --gps ../test/data/gps.csv \
  --output /tmp/polymatch_out.csv
```

The output should be identical (modulo column ordering) to:

```bash
./weightmatch \
  --network ../test/data/edges.shp \
  --gps ../test/data/gps.csv \
  --output /tmp/weightmatch_out.csv
```

A diff check is run as part of `polymatch_test`.

## Full example: polygon-aware matching

```bash
./polymatch \
  --network ../test/data/polymatch/network.shp \
  --polygons ../test/data/polymatch/polygons.shp \
  --access_points ../test/data/polymatch/access_points.shp \
  --gps ../test/data/polymatch/gps.csv \
  --output /tmp/polymatch_out.csv \
  -k 8 -r 300 -e 50 \
  --through_penalty_factor 1.5 \
  --use_omp
```

This will:
1. Load the road network, polygon layer, and access point layer (validating all three).
2. Build the `PolyLinkGraph` polygon-aware routing graph (including the precomputed through-routing
   cost table).
3. Match each GPS trajectory in parallel across OpenMP threads (FR-017) and write a CSV with both
   edge and polygon segments.

**Note**: `--use_omp` is recommended for production input batches; results are bit-identical to a
single-threaded run (SC-009).

## XML config alternative

```bash
./polymatch --config ../test/data/polymatch/config.xml
```

See `contracts/polymatch-cli.md` for the XML schema.

## Output inspection

```bash
head -2 /tmp/polymatch_out.csv
# id;opath;cpath;polygon_ids;entry_aps;egress_aps;is_through;mgeom
# 1001;1,2,-7,-7,3,4;1,2,-7,3,4;7;1001;2002;0;LINESTRING(...)
```

See `contracts/polymatch-output.md` for column definitions.

## Verifying spec acceptance criteria

| Criterion | Verification |
|-----------|--------------|
| **SC-001** (cross-polygon trip matched correctly) | `polymatch_test` golden fixture `gps_cross.csv` |
| **SC-002** (zero regression vs weightmatch) | `polymatch_test` link-only diff test |
| **SC-003** (1000 pts × 100 polys in < 10s) | `polymatch_test` perf benchmark (skipped by default; run with `--bench`) |
| **SC-004** (existing tests pass) | `make tests && for t in algorithm fmm network network_graph weightmatch; do ./test/${t}_test; done` |
| **SC-005** (polygon segments identifiable in output) | Manual inspection of CSV; covered by `polymatch_test` |
| **SC-006** (4 cost types verifiable) | `polymatch_test` cost-model unit tests |
| **SC-007** (emission distance rule) | `polymatch_test` emission-distance unit tests |
| **SC-008** (polygon distance output) | `polymatch_test` against hand-computed fixtures for each segment type |
| **SC-009** (parallel determinism) | `polymatch_test` runs same fixture with N ∈ {1,2,4,8} threads; sorts by trajectory ID; diff = 0 |

## TDD-first development checklist

Per Constitution Principle III, all algorithmic work follows TDD:

1. [ ] Write `polymatch_test.cpp` skeleton + fixture data **first** before any algorithm code.
2. [ ] Implement `PolygonLayer` + `AccessPointLayer` with their validation tests.
3. [ ] Implement `PolyLinkGraph` with construction + Dijkstra unit tests.
4. [ ] Implement `POLYMATCH::match_traj` with HMM tests (entry, egress, within, through, mid-polygon).
5. [ ] Implement `PolyMMWriter` with output-format tests **and a thread-safety test** (FR-017): two
       threads write 1000 result rows concurrently with no corruption or interleaving.
6. [ ] Implement `POLYMATCHApp` OpenMP outer loop with the parallel-determinism test (SC-009).
7. [ ] Run all 6 test executables — all must pass before any PR is merged.

## Common errors

| Error | Likely cause | Fix |
|-------|--------------|-----|
| `Access point feature X not on polygon Y boundary` | Geometry mismatch beyond `boundary_epsilon`. | Increase `--boundary_epsilon` or fix the shapefile. |
| `Polygon ID Z referenced by access point but not in polygon shapefile` | Orphaned access point. | Add the polygon or remove the AP feature. |
| `Polygon Q has no access points; excluded from matching` | Warning only. | Add ≥ 1 access point feature, or accept that Q is unreachable. |
