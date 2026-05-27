# Implementation Plan: PolyMatch — Link-Polygon Map Matching Algorithm

**Branch**: `001-polymatch-algorithm` | **Date**: 2026-05-26 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `/workspace/specs/001-polymatch-algorithm/spec.md`

## Summary

PolyMatch is a new map matching algorithm modeled after WeightMatch. It extends the link-graph routing model
with first-class **polygon** entities (parking lots, ferry terminals, pedestrian plazas) and **access points**
that connect polygons to the road network and to each other. GPS trajectories are matched to a hybrid path
of links and polygons; polygon traversals are recorded with the access points used to enter and exit.

**Technical approach** (per user input):
- Mirror WeightMatch's file layout, config structure, and CLI argument shape (`src/mm/polymatch/` + `src/app/polymatch.cpp` + `polymatch_test.cpp`).
- Do **not** modify FMM, STMATCH, WEIGHTMATCH, or H3MM source files. All polygon-aware logic is additive:
  new types in `FMM::NETWORK` (polygon + access point loading), a new polygon-aware routing graph in
  `FMM::ROUTING`, and a new matcher in `FMM::MM::POLYMATCH`.
- Reuse `DijkstraState` / `IndexedMinHeap` per Constitution Principle I (no per-query allocation).
- Add a new test executable `polymatch_test` following the `weightmatch_test` pattern.

## Technical Context

**Language/Version**: C++11 (strict; per Constitution Principle V & Technical Constraints).

**Primary Dependencies**: GDAL ≥ 2.2 (shapefile I/O), Boost ≥ 1.56 (Graph, Geometry, Serialization), OpenMP (parallel matching). SWIG is **not** used — polymatch ships no Python bindings (per spec assumption).

**Storage**: ESRI shapefiles (network, polygon layer, access point layer), CSV (GPS trajectory input, matched results output), XML (optional config file). No databases.

**Testing**: Custom C++ test executables built via CMake (`make tests`). New test executable `polymatch_test` mirrors `weightmatch_test` — golden-path comparison against fixture data.

**Target Platform**: Linux + macOS. Windows-via-cygwin best-effort (matches the existing project's support matrix).

**Project Type**: C++ library + CLI executables (single-project, monolithic CMake).

**Performance Goals**: SC-003 — 1,000 GPS points × 100-polygon network in < 10 s single-core, plus ≥ 0.7 × N speedup with N OpenMP threads on independent-trajectory batches. SC-009 — parallel results are bit-identical to single-threaded results. Constitution Principle I — no heap allocation per Dijkstra query; OpenMP-parallel trajectory loop is a **v1 requirement** (FR-017), not deferred. Shared graph data is read-only post-construction; per-trajectory mutable state is per-thread.

**Constraints**:
- No regressions to FMM / STMATCH / WEIGHTMATCH / H3MM behavior (spec FR-012, SC-002, SC-004).
- All four existing test suites (`algorithm_test`, `fmm_test`, `network_test`, `network_graph_test`) must continue to pass without modification.
- No virtual dispatch in HMM inner loop (Constitution Technical Constraints).

**Scale/Scope**: Same envelope as WeightMatch — networks up to ~1M edges, trajectories up to ~10⁵ points, polygon layers up to ~10⁴ polygons with ~10⁵ access point relationships. OpenMP-parallel throughput is required (FR-017, SC-003, SC-009) to handle production batch sizes.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Compliance | Notes |
|-----------|------------|-------|
| **I. Performance-First** | ✅ Compliant | Reuse `DijkstraState` + `IndexedMinHeap` for polygon-augmented routing; document Big-O for new polygon transitions (O((E+P) log (E+P)) where P = polygons); OpenMP outer loop preserved. |
| **II. Algorithmic Correctness** | ✅ Compliant | Hybrid `C_Path` topology invariant restated for polygon transitions (spec FR-015); zero/finite checks remain; Viterbi termination unchanged. Duplicate GPS points handled via existing TransitionGraph logic. |
| **III. Test-Driven (NON-NEGOTIABLE)** | ✅ Compliant — gated by Phase 1 design output | `polymatch_test.cpp` added; covers (a) normal trip through polygon, (b) trip entirely inside polygon, (c) trip beginning/ending mid-polygon, (d) through-routing case, (e) duplicate points, (f) out-of-network points, (g) link-only regression vs WeightMatch (SC-002). One integration test against a shapefile fixture in `test/data/`. |
| **IV. Dual-API Stability** | ✅ N/A | Spec explicitly excludes Python bindings (assumption: "No Python API bindings are provided for polymatch"). `python/fmm.i` is **not** modified. |
| **V. Separation of Concerns** | ✅ Compliant | Polygon + access-point loading lives in `FMM::NETWORK`. Polygon-aware routing lives in `FMM::ROUTING` (extends or sits alongside `LinkGraph`). HMM matching lives in `FMM::MM::POLYMATCH`. `FMM::CONFIG` adds `PolygonConfig` + `AccessPointConfig` leaf types. No back-references introduced. |

**Gate result**: PASS. No violations to track in the Complexity Tracking table.

## Project Structure

### Documentation (this feature)

```text
specs/001-polymatch-algorithm/
├── plan.md              # This file
├── spec.md              # Feature spec (ratified)
├── research.md          # Phase 0 output (technical decisions)
├── data-model.md        # Phase 1 output (entities, fields, relationships)
├── quickstart.md        # Phase 1 output (build + run guide)
├── contracts/           # Phase 1 output (CLI command schema, file formats)
│   ├── polymatch-cli.md
│   ├── polygon-shapefile.md
│   ├── access-point-shapefile.md
│   └── polymatch-output.md
├── checklists/
│   └── requirements.md  # Spec quality checklist (passing)
└── tasks.md             # Phase 2 output (created by /speckit-tasks)
```

### Source Code (repository root)

```text
src/
├── network/
│   ├── network.hpp / .cpp              # (existing — unchanged)
│   ├── network_graph.hpp / .cpp        # (existing — unchanged)
│   ├── link_graph_routing.hpp / .cpp   # (existing — unchanged; new file mirrors pattern)
│   ├── polygon_layer.hpp / .cpp        # NEW — loads polygon shapefile via GDAL
│   ├── access_point_layer.hpp / .cpp   # NEW — loads access point shapefile; resolves link attachment
│   └── poly_link_graph.hpp / .cpp      # NEW — polygon-aware routing graph in FMM::ROUTING
├── config/
│   ├── network_config.hpp / .cpp       # (existing — unchanged)
│   ├── gps_config.hpp / .cpp           # (existing — unchanged)
│   ├── result_config.hpp / .cpp        # (existing — unchanged)
│   ├── polygon_config.hpp / .cpp       # NEW — polygon layer config (file path, ID + cost columns)
│   └── access_point_config.hpp / .cpp  # NEW — access point layer config (file path, node + polygon ID columns)
├── mm/
│   ├── transition_graph.hpp / .cpp     # (existing — unchanged)
│   ├── mm_type.hpp                     # (existing — unchanged)
│   ├── fmm/                            # (existing — unchanged)
│   ├── stmatch/                        # (existing — unchanged)
│   ├── h3mm/                           # (existing — unchanged)
│   ├── weightmatch/                    # (existing — unchanged)
│   └── polymatch/                      # NEW — mirrors weightmatch/
│       ├── polymatch_algorithm.hpp / .cpp   # POLYMATCH class + POLYMATCHConfig
│       ├── polymatch_app.hpp / .cpp         # POLYMATCHApp orchestrator
│       ├── polymatch_app_config.hpp / .cpp  # POLYMATCHAppConfig (CLI + XML parse)
│       └── poly_match_result.hpp / .cpp     # PolyMatchResult (hybrid path with access points)
├── io/
│   ├── mm_writer.hpp                   # (existing — unchanged)
│   └── poly_mm_writer.hpp / .cpp       # NEW — CSV writer for hybrid results (link + polygon + access points)
└── app/
    ├── fmm.cpp, stmatch.cpp, weightmatch.cpp, h3mm.cpp, ubodt_gen_app.cpp  # (existing — unchanged)
    └── polymatch.cpp                   # NEW — main() mirrors weightmatch.cpp

test/
├── algorithm_test.cpp, fmm_test.cpp, network_test.cpp,
│   network_graph_test.cpp, weightmatch_test.cpp   # (existing — unchanged)
├── polymatch_test.cpp                              # NEW — golden-path tests
├── data/
│   ├── (existing fixtures unchanged)
│   └── polymatch/                                  # NEW — polygon + access point + GPS fixtures
└── CMakeLists.txt                                  # MODIFIED — add polymatch_test target only

CMakeLists.txt                                      # MODIFIED — add POLYMATCH_OBJ + polymatch executable
```

**Structure Decision**: Single-project layout mirroring the existing repo. All new files are additive
under new directories (`src/mm/polymatch/`) or new filenames (`polygon_layer.hpp`, `poly_link_graph.hpp`,
etc.). The only modifications to existing files are **additive declarations** in two CMakeLists.txt files
— no source-file edits to FMM, STMATCH, WEIGHTMATCH, H3MM, or shared infrastructure.

## Phase 0 — Research

See [research.md](./research.md). Key decisions resolved:

1. **Polygon-aware routing graph design (R1)**: Augment `LinkGraph` semantics into a new `PolyLinkGraph`
   that adds polygon vertices and access-point arcs to the existing edge-as-vertex graph. Reuses the
   `shortest_edge_to_edges` Dijkstra core via a generalized vertex ID space.
2. **Polygon + access-point shapefile loading (R2, R3)**: GDAL OGR driver, same pattern as
   `Network::read_ogr_file`. Three-condition validation pass on access-point shapefile (FR-005).
3. **Access-point-to-link attachment (R4)**: Determined by **direct node-ID lookup** against the
   network's existing `node_id → NodeIndex` map. The access-point shapefile's `node_id` field shares
   the road network's node ID scheme by design, so no spatial coincidence check is needed. O(1) lookup,
   no EPSILON tolerance for link attachment. `boundary_epsilon` is used only for the polygon-boundary
   geometry validation.
4. **Polygon transition cost computation strategy (R5)**: Costs are split into observation-dependent
   vs. precomputed. Entry, egress, and within-polygon (same-link override) costs are computed at
   routing time because they depend on the matched-point-inside-polygon or per-GPS-pair Euclidean
   distances. **Through-routing** cost is precomputed at graph construction (see #8 below).
5. **Output extension (R6)**: New `PolyMMWriter` produces a CSV with the same columns as
   `CSVMatchResultWriter` plus polygon segment columns (`polygon_ids`, `entry_aps`, `egress_aps`,
   `is_through`, `polygon_distances`). PolyMatch uses `PolyMMWriter` exclusively; existing
   `CSVMatchResultWriter` is unchanged.
6. **Mid-polygon initialization/termination (R7)**: HMM initial/terminal state allows polygon
   candidates; absent access points are represented as empty strings in the output and as a sentinel
   `kNoAccessPoint` value internally.
7. **Result type extension (R8)**: A new `PolyMatchResult` struct wraps the existing `MatchResult` and
   adds polygon segment metadata (polygon ID, entry AP, egress AP, through-routing flag,
   distance-inside). The existing `MatchResult` type is **not** modified.
8. **Boundary point handling and fallback (R9, R10)**: Boost.Geometry `covered_by` treats boundary
   points as inside (emission distance zero). Empty polygon layer falls back to a pure
   `LinkGraph`-equivalent routing graph and link-only matching (FR-012).
9. **Through-routing cost precomputation (R11)**: At `PolyLinkGraph` construction, build a per-polygon
   table of `polygon.weight × distance(AP_i, AP_j)` for every ordered access-point pair. Final
   through-routing cost = `precomputed[i,j] × THROUGH_PENALTY_FACTOR` at lookup time. Storing the
   factor-independent quantity preserves runtime tunability of the factor without rebuilding the graph.
   Space: O(Σ_p n_p²) — negligible for realistic inputs.
10. **Distance-inside-polygon computation (R12)**: Computed at result-assembly time from
    already-cached `eu_dist` values used by the HMM. O(P) extra work per trajectory where P is the
    number of polygon segments matched.
11. **OpenMP-parallel matching design (R13, FR-017)**: Outer trajectory loop is parallelized. Shared
    graph data (Network, PolygonLayer, AccessPointLayer, PolyLinkGraph including the through-routing
    table) is read-only post-construction. Per-thread mutable state (DijkstraState, IndexedMinHeap,
    TransitionGraph) is created inside the parallel region. Writer uses an internal mutex for
    thread-safe append; matched output content is deterministic regardless of thread count (SC-009).

## Phase 1 — Design Artifacts

Generated files (see each for full content):

- **[data-model.md](./data-model.md)** — Polygon, AccessPoint, PolyLinkGraph, PolyMatchCandidate,
  PolyMatchResult, PolygonSegment, POLYMATCHConfig, PolygonConfig, AccessPointConfig.
- **[contracts/polymatch-cli.md](./contracts/polymatch-cli.md)** — `polymatch` CLI argument schema
  (mirrors `weightmatch` plus polygon-specific flags).
- **[contracts/polygon-shapefile.md](./contracts/polygon-shapefile.md)** — polygon shapefile input
  format (geometry, required fields, optional fields).
- **[contracts/access-point-shapefile.md](./contracts/access-point-shapefile.md)** — access point
  shapefile input format (one feature per polygon-AP relationship).
- **[contracts/polymatch-output.md](./contracts/polymatch-output.md)** — output CSV schema (extended
  with polygon segment columns).
- **[quickstart.md](./quickstart.md)** — build + run guide on a worked example.

**Agent context update**: The `<!-- SPECKIT START --> ... <!-- SPECKIT END -->` block in
`/workspace/CLAUDE.md` is updated to point at this plan.

## Post-Design Constitution Re-check

Re-running the gate after design output is available:

| Principle | Compliance after design | Verification |
|-----------|-------------------------|--------------|
| I. Performance-First | ✅ | `PolyLinkGraph` reuses `DijkstraState`/`IndexedMinHeap`; data-model.md documents Big-O of polygon-augmented Dijkstra. |
| II. Algorithmic Correctness | ✅ | data-model.md enumerates invariants: `C_Path` topology with access points, finite costs, Viterbi cycles forbidden. |
| III. Test-Driven Development | ✅ | quickstart.md includes a TDD-first checklist; test fixtures specified in tasks.md. |
| IV. Dual-API Stability | ✅ N/A | No `python/fmm.i` changes; new types are C++-only. |
| V. Separation of Concerns | ✅ | Directory layout above keeps `FMM::NETWORK`, `FMM::ROUTING`, `FMM::MM`, `FMM::CONFIG` cleanly separated; no back-references introduced. |

**Gate result**: PASS. Ready for `/speckit-tasks`.

## Complexity Tracking

No constitution violations to justify. Table intentionally empty.

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|--------------------------------------|
| (none) | (n/a) | (n/a) |
