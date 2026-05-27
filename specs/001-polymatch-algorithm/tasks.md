---

description: "Task list for PolyMatch implementation"
---

# Tasks: PolyMatch — Link-Polygon Map Matching Algorithm

**Input**: Design documents from `/workspace/specs/001-polymatch-algorithm/`

**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/

**Tests**: TDD is mandatory per Constitution Principle III (NON-NEGOTIABLE). Tests are written FIRST and MUST FAIL before implementation.

**Organization**: Tasks are grouped by user story. Note that for this feature the user stories partition the codebase by **concern** (input, matching, output) rather than delivering independent value increments — the MVP requires all three working together. Implementation order is US2 → US1 → US3 because US1's matcher requires US2's loaders.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Different files, no dependencies on incomplete tasks
- **[Story]**: Maps to user stories from spec.md (US1, US2, US3)
- File paths are absolute under `/workspace/`

## Path Conventions

Single C++ project; source under `src/`, tests under `test/`, app entry points under `src/app/`.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Directory creation, CMake additions, test fixture generation

- [X] T001 Create source directory `/workspace/src/mm/polymatch/` for new algorithm files
- [X] T002 Create test fixture directory `/workspace/test/data/polymatch/` for polymatch-specific test data
- [ ] T003 [P] Create Python fixture generator at `/workspace/test/generate_polymatch_test_data.py` that emits a small road network shapefile (~20 edges, 4 intersections), a polygon shapefile (3 polygons including one self-intersecting invalid case and one polygon-with-no-AP case), an access point shapefile (covering link-attached, shared-between-polygons, off-boundary error, orphaned polygon ID error, and contradictory-geometry error cases), and a GPS CSV with trajectories covering all US1 scenarios (link-only, cross-polygon, two-polygons-shared-AP, mid-polygon start, mid-polygon end, fully-inside, through-routing)
- [X] T004 Add `POLYMATCH_OBJ` object library to `/workspace/CMakeLists.txt`: `file(GLOB POLYMATCHGlob src/mm/polymatch/*.cpp)` then `add_library(POLYMATCH_OBJ OBJECT ${POLYMATCHGlob})`; include `$<TARGET_OBJECTS:POLYMATCH_OBJ>` in the `FMMLIB` shared library composition (mirror the `WEIGHTMATCH_OBJ` pattern)
- [X] T005 Add `polymatch` executable target to `/workspace/CMakeLists.txt`: `add_executable(polymatch src/app/polymatch.cpp)` + `target_link_libraries(polymatch FMMLIB)` (mirror the `weightmatch` target)
- [X] T006 Add `polymatch_test` target to `/workspace/test/CMakeLists.txt` mirroring `weightmatch_test`: include `$<TARGET_OBJECTS:POLYMATCH_OBJ>` plus MM_OBJ, CORE, CONFIG, network, io object libs; link GDAL_LIBRARIES, Boost_LIBRARIES, OpenMP_CXX_LIBRARIES, OSMIUM_LIBRARIES

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Pure data types, configuration structs, constants. No business logic — these are header-mostly files that downstream stories depend on.

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [X] T007 [P] Create `PolygonConfig` struct in `/workspace/src/config/polygon_config.hpp` with fields `file` (string), `id_name` (string, default `"id"`), `cost_name` (string, default `"cost"`); declare `validate()` and ptree XML reader
- [X] T008 [P] Implement `PolygonConfig::validate()` and `read_xml()` in `/workspace/src/config/polygon_config.cpp` mirroring `NetworkConfig`
- [X] T009 [P] Create `AccessPointConfig` struct in `/workspace/src/config/access_point_config.hpp` with fields `file` (string), `node_id_name` (string, default `"node_id"`), `polygon_id_name` (string, default `"polygon_id"`); declare `validate()` and ptree XML reader
- [X] T010 [P] Implement `AccessPointConfig::validate()` and `read_xml()` in `/workspace/src/config/access_point_config.cpp`
- [X] T011 [P] Create `Polygon` struct + `PolygonLayer` class declaration in `/workspace/src/network/polygon_layer.hpp` per data-model.md (fields: index, id, geom, weight, bbox; PolygonLayer container with `polygons`, `id_to_index`, `rtree`)
- [X] T012 [P] Create `AccessPointFeature`, `AccessPoint`, `AccessPointLayer` class declarations in `/workspace/src/network/access_point_layer.hpp` per data-model.md (AccessPoint has `index`, `node_id`, `point`, `polygons`, `attached_node`, `attached_edges`; AccessPointLayer container with lookups)
- [X] T013 [P] Declare `PolyLinkGraph` class + `shortest_polylink_to_polylinks()` routing function in `/workspace/src/network/poly_link_graph.hpp` per data-model.md including `n_edges`, `n_polygons`, `adjacency`, `through_cost_tables`, `polygon_local_ap_index`, `vertex_kind()`
- [X] T014 [P] Create `POLYMATCHConfig` struct in `/workspace/src/mm/polymatch/polymatch_algorithm.hpp` mirroring `WEIGHTMATCHConfig` field-for-field plus `through_penalty_factor` (default `1.5`) and `boundary_epsilon` (default `1e-6`); declare `validate()`, `register_arg()`, `load_from_arg()`, `register_help()`, `print()`
- [X] T015 [P] Implement `POLYMATCHConfig` methods (`validate`, `register_arg`, `load_from_arg`, `register_help`, `print`) in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp` mirroring `WEIGHTMATCHConfig` implementation
- [X] T016 [P] Create `PolygonSegment` struct + `PolyMatchResult` struct + `kNoAccessPoint` sentinel constant (`= -1`) in `/workspace/src/mm/polymatch/poly_match_result.hpp` per data-model.md (fields: polygon_id, entry_ap, egress_ap, is_through, distance_inside, position_in_cpath)
- [X] T017 [P] Declare `POLYMATCH` matcher class skeleton in `/workspace/src/mm/polymatch/polymatch_algorithm.hpp`: constructor signature `POLYMATCH(Network&, PolygonLayer&, AccessPointLayer&, PolyLinkGraph&)`, `match_traj()` method signature mirroring `WEIGHTMATCH::match_traj`
- [X] T018 [P] Declare `POLYMATCHAppConfig` class in `/workspace/src/mm/polymatch/polymatch_app_config.hpp` containing `NetworkConfig`, `GPSConfig`, `ResultConfig`, `PolygonConfig`, `AccessPointConfig`, `POLYMATCHConfig`, plus `use_omp`, `help_specified`, `log_level`, `step` (mirror `WEIGHTMATCHAppConfig`)
- [X] T019 [P] Declare `POLYMATCHApp` orchestrator class in `/workspace/src/mm/polymatch/polymatch_app.hpp` taking `POLYMATCHAppConfig const&`; `run()` method declaration

**Checkpoint**: All headers compile. CMake configures without errors. No business logic yet.

---

## Phase 3: User Story 2 (Priority: P2) — Configure Polygon Layer as Routing Input

**Goal**: PolyMatch loads polygon shapefile + access point shapefile, validates them per FR-005, applies FR-013/FR-014 graceful degradation, and falls back to link-only mode per FR-012.

**Independent Test**: `polymatch --polygons X --access_points Y --network N --gps G --output O` either succeeds in loading inputs (printing "loaded N polygons / M access points") or halts with a descriptive error for invalid inputs. Tests in `polymatch_test` exercise each validation path.

### Tests for User Story 2 (TDD — write FIRST, ensure FAIL)

- [ ] T020 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolygonLayer, loads_valid_shapefile)`: load test fixture polygons, assert correct count, IDs, weights, R-tree built
- [ ] T021 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolygonLayer, skips_invalid_geometry)`: fixture includes a self-intersecting polygon; assert it is skipped with warning, valid polygons remain (FR-013)
- [ ] T022 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolygonLayer, point_in_polygon_and_boundary)`: assert `polygons_containing(p)` returns polygon for point inside, on boundary, and empty for point outside; assert `min_boundary_distance` returns 0 inside/boundary, positive outside (FR-006)
- [ ] T023 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, loads_and_deduplicates_shared_node)`: fixture includes a node_id shared across 2 polygons; assert 1 resolved AccessPoint with `polygons.size() == 2`
- [ ] T024 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, rejects_off_boundary_ap)`: feature whose geometry is > `boundary_epsilon` from declared polygon's boundary; assert load throws / halts with descriptive error (FR-005 condition 1)
- [ ] T025 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, rejects_orphaned_polygon_ref)`: feature references polygon_id not in PolygonLayer; assert halt (FR-005 condition 2)
- [ ] T026 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, rejects_contradictory_geometries)`: two features with same node_id but differing geometry; assert halt (FR-005 condition 3)
- [ ] T027 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, link_attachment_via_node_id)`: AP node_id matches a network source/target node ID; assert `attached_node.has_value()` and `attached_edges` is populated (R4)
- [ ] T028 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, polygon_shared_only_no_link_attachment)`: AP node_id does not match any network node; assert `attached_node.has_value() == false` and the AP is still valid because `polygons.size() >= 2`
- [ ] T029 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(AccessPointLayer, polygon_with_no_ap_warned_and_excluded)`: polygon has no features in AP shapefile; assert warning emitted and polygon excluded from candidate generation (FR-014)
- [ ] T030 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(POLYMATCHAppConfig, cli_parses_polygon_flags)`: argc/argv with `--polygons`, `--polygon_id_name`, `--polygon_cost_name`, `--access_points`, `--ap_node_id_name`, `--ap_polygon_id_name`, `--through_penalty_factor`, `--boundary_epsilon`; assert values land in correct config fields
- [ ] T031 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(POLYMATCHAppConfig, xml_parses_polygon_blocks)`: XML config with `<polygon>` and `<access_point>` blocks parses into matching fields
- [ ] T032 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(POLYMATCHApp, fallback_link_only_when_no_polygon_layer)`: instantiate POLYMATCHApp with empty polygon_config; assert matching proceeds as link-only and output equivalent to weightmatch (FR-012, SC-002)
- [ ] T033 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(POLYMATCHApp, fallback_link_only_when_zero_valid_polygons)`: polygon layer present but all polygons invalid (skipped); assert fallback to link-only (FR-012)
- [ ] T086 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolygonLayer, rejects_polygon_id_zero)`: fixture includes a polygon feature with `id == 0`; assert load halts with a descriptive error citing ID 0 (FR-018)
- [ ] T087 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolygonLayer, rejects_duplicate_polygon_id)`: fixture includes two polygons sharing the same ID; assert load halts with a descriptive error listing the offending ID (FR-018)
- [ ] T089 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST(POLYMATCHAppConfig, cli_parity_with_weightmatch)`: parse `polymatch --help` and `weightmatch --help` outputs; assert every flag accepted by weightmatch is also accepted by polymatch (excluding any polygon-specific additions) (FR-011)

### Implementation for User Story 2

- [ ] T034 [P] [US2] Implement `PolygonLayer::load()` in `/workspace/src/network/polygon_layer.cpp`: open shapefile via `GDALOpenEx`, iterate features, parse id + cost fields (using `PolygonConfig` field names), validate via `boost::geometry::is_valid` (skip self-intersect/zero-area with `SPDLOG_WARN` per FR-013), convert to Boost.Geometry polygon, populate `polygons`, `id_to_index`, `rtree`. Hard-reject per FR-018: halt with a descriptive error if any polygon feature has `id == 0`, or if two features share the same `id` (error message MUST list the offending ID).
- [ ] T035 [US2] Implement spatial queries `PolygonLayer::polygons_containing(point)`, `PolygonLayer::polygons_within_radius(point, radius)`, `PolygonLayer::min_boundary_distance(polygon_idx, point)` in `/workspace/src/network/polygon_layer.cpp` using Boost.Geometry `covered_by` (boundary = inside per R9), `distance`, and the R-tree (depends on T034)
- [ ] T036 [P] [US2] Implement `AccessPointLayer::load(AccessPointConfig, PolygonLayer&, Network&)` in `/workspace/src/network/access_point_layer.cpp`: GDAL OGR load of features → build vector of `AccessPointFeature`; run three-condition validation (FR-005); group by `node_id` and deduplicate to `AccessPoint`s; resolve `attached_node` via direct `Network::node_id_to_index` lookup with O(1) hash lookup (R4); precompute `attached_edges` by scanning Network edges; build `node_id_to_index` and `polygon_to_aps` maps; emit per-polygon warnings for polygons absent from `polygon_to_aps` (FR-014)
- [ ] T037 [US2] Implement `POLYMATCHAppConfig::POLYMATCHAppConfig(int argc, char**)` and `load_arg()` in `/workspace/src/mm/polymatch/polymatch_app_config.cpp` mirroring `WEIGHTMATCHAppConfig::load_arg`: register all existing weightmatch flags via shared configs' `register_arg`, then add polygon-specific flags (`--polygons`, `--polygon_id_name`, `--polygon_cost_name`, `--access_points`, `--ap_node_id_name`, `--ap_polygon_id_name`) and POLYMATCHConfig flags (T015)
- [ ] T038 [US2] Implement `POLYMATCHAppConfig` XML parsing in `/workspace/src/mm/polymatch/polymatch_app_config.cpp`: detect `--config` flag, parse XML via ptree, populate sub-configs including new `<polygon>` and `<access_point>` blocks per `contracts/polymatch-cli.md` (depends on T037)
- [ ] T039 [US2] Implement `POLYMATCHAppConfig::validate()` in `/workspace/src/mm/polymatch/polymatch_app_config.cpp`: validate all sub-configs; allow `polygon_config.file` and `access_point_config.file` to be empty (triggers fallback); require both or neither (validate that AP layer is provided iff polygon layer is provided) (depends on T037)
- [ ] T040 [US2] Implement `POLYMATCHAppConfig::print_help()` in `/workspace/src/mm/polymatch/polymatch_app_config.cpp` mirroring `WEIGHTMATCHAppConfig::print_help()` with new polygon-related flag documentation (depends on T037)
- [ ] T041 [US2] Implement `POLYMATCHApp::run()` skeleton in `/workspace/src/mm/polymatch/polymatch_app.cpp`: load Network (existing); if `polygon_config.file` non-empty, load PolygonLayer (T034) and AccessPointLayer (T036); if both empty or zero valid polygons after load, log info and set internal `link_only_mode` flag for fallback (FR-012); stop short of actual matching (matching belongs to US1) (depends on T034, T036)

**Checkpoint**: `polymatch --polygons P.shp --access_points AP.shp --network N.shp --gps G.csv --output O.csv` validates and loads inputs successfully or halts with descriptive errors. All T020-T033 tests pass.

---

## Phase 4: User Story 1 (Priority: P1) — Match GPS Trip Through Mixed Link-Polygon Network 🎯 MVP

**Goal**: Core matcher — given loaded inputs (US2), match GPS trajectories through a hybrid link+polygon network and produce in-memory `PolyMatchResult` with hybrid C_Path and polygon segments carrying entry/egress access points, through-routing flag, and distance_inside.

**Independent Test**: `polymatch_test` golden-path fixtures verify each US1 acceptance scenario by comparing `PolyMatchResult` fields against expected values. Link-only fallback produces results identical to weightmatch on the same input (SC-002).

### Tests for User Story 1 (TDD — write FIRST, ensure FAIL)

- [ ] T042 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyLinkGraph, construction_vertex_and_arc_counts)`: build graph for fixture; assert vertex count `|E|+|P|`, link→link arcs preserved, link↔polygon arcs added for every AP with `attached_node`, polygon↔polygon arcs for every shared-AP pair
- [ ] T043 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyLinkGraph, through_cost_table_correctness)`: polygon with 3 APs at known coords and weight `w`; assert `through_cost_tables[p][i*n+j] == w * dist(AP_i, AP_j)` for all pairs; assert final cost = `precomputed × THROUGH_PENALTY_FACTOR` (SC-006, R11)
- [ ] T044 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyLinkGraph, through_penalty_factor_varies_cost_proportionally)`: same fixture, factor in {1.0, 1.5, 5.0}; assert through-routing cost scales linearly without rebuilding the graph (R11)
- [ ] T045 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(Routing, dijkstra_polylinkgraph_reuses_state_no_allocation)`: run 1000 Dijkstra queries with the same `DijkstraState`/`IndexedMinHeap`; instrument `mallocs` and assert zero allocations during the query phase (Constitution Principle I)
- [ ] T046 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(Emission, distance_zero_for_point_inside_polygon)`: GPS point strictly inside; assert emission distance == 0 (FR-006, SC-007)
- [ ] T047 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(Emission, distance_zero_for_point_on_polygon_boundary)`: GPS point exactly on boundary; assert emission distance == 0 (FR-006, R9)
- [ ] T048 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(Emission, distance_equals_min_boundary_for_point_outside)`: GPS point at known coords outside; assert emission distance == min distance to polygon boundary (FR-006, SC-007)
- [ ] T049 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, link_only_matches_weightmatch_exactly)`: link-only GPS trajectory on shared fixture; assert `PolyMatchResult.base.opath` and `.cpath` equal `WEIGHTMATCH::match_traj` result (US1 scenario 2, SC-002)
- [ ] T050 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, crosses_polygon_records_entry_egress)`: trajectory enters polygon, crosses, exits to different link; assert polygon segment present with correct polygon_id, non-empty entry_ap, non-empty egress_ap, `is_through == false` (US1 scenario 1)
- [ ] T051 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, entirely_inside_polygon_single_segment)`: trajectory fully inside; assert single polygon segment with `entry_ap == kNoAccessPoint && egress_ap == kNoAccessPoint` (US1 scenario 3, FR-010)
- [ ] T052 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, two_polygons_via_shared_ap)`: trajectory passes through two polygons sharing a node_id; assert two consecutive polygon segments where the egress AP of segment 1 equals the entry AP of segment 2 (US1 scenario 4)
- [ ] T053 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, through_routing_no_gps_inside)`: GPS configured so optimal path passes through polygon but no GPS observation falls inside it; assert polygon segment with `is_through == true`, both APs recorded, and applied cost matches `weight × dist(entry_AP, egress_AP) × THROUGH_PENALTY_FACTOR` (US1 scenario 5, FR-008, SC-006)
- [ ] T054 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, mid_polygon_start_no_entry_ap)`: trajectory begins inside polygon; assert first polygon segment has `entry_ap == kNoAccessPoint`, valid `egress_ap` (FR-007, FR-010)
- [ ] T055 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, mid_polygon_end_no_egress_ap)`: trajectory ends inside polygon; assert last polygon segment has valid `entry_ap`, `egress_ap == kNoAccessPoint` (FR-007, FR-010)
- [ ] T056 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, within_polygon_uses_eu_dist_override)`: two consecutive GPS points both inside same polygon; instrument `update_layer` and assert it set `path_distance = eu_dist`, mirroring the same-link override at `weightmatch_algorithm.cpp:323-324` (FR-008)
- [ ] T057 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, duplicate_gps_points_no_nan)`: trajectory with consecutive zero-distance GPS points; assert no NaN/inf in HMM transition probabilities, matcher returns a valid result (Constitution Principle II)
- [ ] T058 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, distance_inside_polygon_correct)`: hand-computed fixtures for each segment type — traversal, through-routing, mid-polygon start, mid-polygon end, fully-inside; assert `PolygonSegment.distance_inside` matches formula within floating-point tolerance (FR-016, SC-008)
- [ ] T059 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(MatchTraj, hybrid_cpath_topology_valid)`: assert every consecutive link-link pair in `C_Path` shares a node; every link-polygon transition uses a shared AP; every polygon-polygon transition uses an AP listed under both polygons (FR-015)
- [ ] T088 [P] [US1] In `/workspace/test/polymatch_test.cpp`, add `TEST(POLYMATCHApp, empty_or_single_point_trajectory_skipped)`: input batch contains one trajectory with 0 GPS points, one with 1 GPS point, and one normal multi-point trajectory; assert the two short trajectories are skipped with per-trajectory warnings identifying their IDs, no output row is produced for them, and the normal trajectory is matched successfully (FR-019)

### Implementation for User Story 1

- [ ] T060 [US1] Implement `PolyLinkGraph::PolyLinkGraph(Network&, LinkGraph&, PolygonLayer&, AccessPointLayer&)` constructor in `/workspace/src/network/poly_link_graph.cpp`: allocate `adjacency` of size `|E|+|P|`; copy link→link arcs from existing `LinkGraph`; for each access point with `attached_node`, add link→polygon arc (from each incident edge) and polygon→link arc (to each incident edge) for every polygon the AP references; for each shared-AP polygon pair, add polygon→polygon arcs (R1) (depends on T013, T034, T036)
- [ ] T061 [US1] Implement the through-routing precomputation in `PolyLinkGraph` constructor in `/workspace/src/network/poly_link_graph.cpp`: for each polygon `p` with access points `[AP_0, ..., AP_n]`, allocate `through_cost_tables[p]` as `n*n` doubles, compute `polygon.weight * distance(AP_i.point, AP_j.point)` for every ordered pair, populate `polygon_local_ap_index[p]` mapping `AccessPointIndex → local index` (R11) (depends on T060)
- [ ] T062 [US1] Implement `shortest_polylink_to_polylinks()` in `/workspace/src/network/poly_link_graph.cpp`: Dijkstra over `PolyLinkGraph` using reused `DijkstraState`/`IndexedMinHeap`; size state to `|E|+|P|` once; during arc relaxation from a polygon vertex, read `parent[P]` to identify the entry-side AP, look up cost in `through_cost_tables[p][entry_local*n + egress_local] * THROUGH_PENALTY_FACTOR`; support `upper_bound_factor` cutoff identical to `shortest_edge_to_edges` (R5, R11) (depends on T061)
- [ ] T063 [US1] Implement `POLYMATCH::POLYMATCH(Network&, PolygonLayer&, AccessPointLayer&, PolyLinkGraph&)` constructor in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp` storing const references (depends on T017)
- [ ] T064 [US1] Implement polygon candidate search in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp` as `search_polygon_candidates(point, POLYMATCHConfig)`: query `PolygonLayer::polygons_containing(point)` for inside candidates (emission distance = 0); query `PolygonLayer::polygons_within_radius(point, radius)` for proximity candidates with emission distance = `min_boundary_distance(p, point)` (FR-006, R9) (depends on T035, T063)
- [ ] T065 [US1] Implement `POLYMATCH::match_traj()` orchestrator in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp` mirroring `WEIGHTMATCH::match_traj`: time-tracked phases — candidate search (link via existing `Network::search_tr_cs_knn` + polygon via T064), build TransitionGraph with hybrid candidates, call `update_tg`, Viterbi backtrack, build C_Path, assemble `PolyMatchResult` (depends on T063, T064)
- [ ] T066 [US1] Implement `POLYMATCH::update_tg()` and `POLYMATCH::update_layer()` in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp` mirroring weightmatch with extensions: (a) same-polygon both-inside override sets `path_distance = eu_dist` (mirror `weightmatch_algorithm.cpp:323-324`); (b) polygon-entry transition cost = `polygon.weight × dist(entry_AP, matched_point_inside_polygon)`; (c) polygon-egress symmetric; (d) other transitions use `shortest_polylink_to_polylinks` (T062) with through-routing cost applied at relaxation time; (e) layer 0 / layer N-1 special case — polygon candidates pay no entry/egress cost (R7) (depends on T062, T065)
- [ ] T067 [US1] Implement `POLYMATCH::build_cpath()` in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp`: walk the Viterbi `TGOpath`, expand to hybrid `C_Path` encoding polygon IDs as negative values (per `contracts/polymatch-output.md`); track per-polygon-segment entry_ap, egress_ap, is_through (true iff no GPS observation matched inside the segment's polygon), and compute distance_inside per R12 (entry-to-first-inside + Σ eu_dist between consecutive inside GPS pairs + last-inside-to-egress; or entry-to-egress for through-routing); populate `PolyMatchResult.polygon_segments` (FR-010, FR-016) (depends on T066)
- [ ] T068 [US1] Complete `POLYMATCHApp::run()` in `/workspace/src/mm/polymatch/polymatch_app.cpp`: construct `PolyLinkGraph` (T060–T061) from loaded inputs; open GPS reader; construct per-thread `DijkstraState`, `IndexedMinHeap`, `TransitionGraph` inside the parallel region; outer trajectory loop annotated `#pragma omp parallel for schedule(dynamic)` if `use_omp` else single-threaded; per trajectory: if the trajectory has 0 or 1 GPS points, emit a per-trajectory warning identifying the trajectory ID, skip it, and continue (no output row produced) per FR-019; otherwise call `POLYMATCH::match_traj` then forward result to writer (writer wiring deferred to US3); honor `link_only_mode` by bypassing polygon candidates (FR-012, FR-017, FR-019, R13) (depends on T041, T067)
- [ ] T069 [US1] Implement `main()` in `/workspace/src/app/polymatch.cpp` mirroring `/workspace/src/app/weightmatch.cpp`: instantiate `POLYMATCHAppConfig(argc, argv)`, handle `--help`, call `validate()`, instantiate `POLYMATCHApp`, call `run()` (depends on T068)

**Checkpoint**: US1 complete. `polymatch` builds, runs end-to-end, produces `PolyMatchResult` matching all US1 acceptance scenarios. Link-only fallback bit-identical to weightmatch (SC-002). All T042–T059 tests pass.

---

## Phase 5: User Story 3 (Priority: P3) — Inspect Matched Path Distinguishing Links from Polygons

**Goal**: `PolyMMWriter` serializes `PolyMatchResult` to CSV with polygon-aware columns (`polygon_ids`, `entry_aps`, `egress_aps`, `is_through`, `polygon_distances`), polygon ID negation in `opath`/`cpath`, thread-safe append, and link-only-mode column omission for SC-002 binary-identical fallback.

**Independent Test**: Verify CSV output column presence, value formatting, polygon ID encoding, empty access point token, through-routing flag, and thread safety against the schema in `contracts/polymatch-output.md`.

### Tests for User Story 3 (TDD — write FIRST, ensure FAIL)

- [ ] T070 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyMMWriter, writes_all_polygon_columns)`: write a `PolyMatchResult` with one polygon segment; parse the output CSV; assert columns `polygon_ids`, `entry_aps`, `egress_aps`, `is_through`, `polygon_distances` are present with correct values (FR-010, FR-016)
- [ ] T071 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyMMWriter, polygon_id_encoding_negative_in_cpath)`: assert polygons appear as negated IDs in `opath`/`cpath`, positive IDs in `polygon_ids` (per `contracts/polymatch-output.md`)
- [ ] T072 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyMMWriter, empty_access_point_token_for_mid_polygon)`: `PolyMatchResult` with `entry_ap == kNoAccessPoint`; assert `entry_aps` column for that segment contains the literal token `-`
- [ ] T073 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyMMWriter, through_routing_flag_value)`: through-routing segment; assert `is_through` value `1`; normal traversal segment; assert `0`
- [ ] T074 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyMMWriter, link_only_fallback_omits_polygon_columns)`: write a `PolyMatchResult` with empty `polygon_segments`; assert output header and rows omit `polygon_ids`/`entry_aps`/`egress_aps`/`is_through`/`polygon_distances` entirely; remaining columns binary-identical to `CSVMatchResultWriter` output (SC-002)
- [ ] T075 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add `TEST(PolyMMWriter, thread_safe_concurrent_writes)`: spawn 2 threads each writing 1000 distinct `PolyMatchResult`s via the same `PolyMMWriter`; parse the resulting CSV and assert all 2000 rows present, no corruption, no interleaving inside a row (FR-017)

### Implementation for User Story 3

- [ ] T076 [P] [US3] Declare `PolyMMWriter` class in `/workspace/src/io/poly_mm_writer.hpp`: constructor `PolyMMWriter(const std::string& filename, const OutputConfig&, bool include_polygon_columns)`; `write_result(const Trajectory&, const PolyMatchResult&)` method; internal `std::mutex write_mutex_` member
- [ ] T077 [US3] Implement `PolyMMWriter::PolyMMWriter()` and header emission in `/workspace/src/io/poly_mm_writer.cpp`: open output file; emit base columns (id, opath, error, offset, spdist, cpath, tpath, mgeom); when `include_polygon_columns`, also emit `polygon_ids;entry_aps;egress_aps;is_through;polygon_distances` (depends on T076)
- [ ] T078 [US3] Implement `PolyMMWriter::write_result()` in `/workspace/src/io/poly_mm_writer.cpp`: acquire `write_mutex_`; serialize base columns from `result.base` with polygon entries negated in opath/cpath per `contracts/polymatch-output.md`; if `include_polygon_columns`, serialize polygon segment columns — `entry_ap`/`egress_ap` use literal token `-` when `kNoAccessPoint`; release mutex on scope exit (FR-017, R13) (depends on T077)
- [ ] T079 [US3] Wire `PolyMMWriter` into `POLYMATCHApp::run()` in `/workspace/src/mm/polymatch/polymatch_app.cpp`: instantiate `PolyMMWriter` with `include_polygon_columns = !link_only_mode`; call `writer.write_result(...)` from inside the OpenMP parallel region after each `match_traj` (depends on T068, T078)

**Checkpoint**: US3 complete. Output CSV format matches `contracts/polymatch-output.md`; link-only fallback omits polygon columns and is byte-for-byte equivalent to weightmatch output. All T070–T075 tests pass.

---

## Phase 6: Polish & Cross-Cutting Concerns

- [ ] T080 [P] In `/workspace/test/polymatch_test.cpp`, add `TEST(Regression, polymatch_vs_weightmatch_link_only_csv_identical)`: run both `polymatch` and `weightmatch` (via the existing weightmatch entry points) on the same link-only fixture trajectory; assert output CSV files are byte-identical after sorting by `id` (SC-002)
- [ ] T081 [P] In `/workspace/test/polymatch_test.cpp`, add `TEST(Performance, baseline_1000_points_100_polygons)` (gated behind a `--bench` flag, default skipped): build fixture with 1000 GPS points and 100 polygons; run single-core; assert wall time < 10 s (SC-003)
- [ ] T082 [P] In `/workspace/test/polymatch_test.cpp`, add `TEST(Parallelism, bit_identical_output_across_thread_counts)`: same input batch matched with OMP_NUM_THREADS in {1, 2, 4, 8}; collect outputs; sort rows by trajectory ID; assert all four outputs byte-identical (SC-009, FR-017)
- [ ] T083 Run all 6 test executables and assert zero failures (Constitution Principle III, SC-004): `cd /workspace/build && make tests && ./test/algorithm_test && ./test/fmm_test && ./test/network_test && ./test/network_graph_test && ./test/weightmatch_test && ./test/polymatch_test`
- [ ] T084 [P] Add doxygen comments documenting Big-O of `PolyLinkGraph` construction (`O(|E| + Σ |AP.polygons| + Σ |AP.attached_links| + Σ_p n_p²)`) and `shortest_polylink_to_polylinks` query (`O((|E|+|P|) log (|E|+|P|))`) directly above the function declarations in `/workspace/src/network/poly_link_graph.hpp` (Constitution Principle I)
- [ ] T085 Manually execute `/workspace/specs/001-polymatch-algorithm/quickstart.md` end-to-end against the generated fixture data; verify the example output rows match the documented schema

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately.
- **Phase 2 (Foundational)**: After Phase 1 (needs directories and CMake declarations). Blocks all user stories.
- **Phase 3 (US2)**: After Phase 2. Implements loaders + configs needed by US1.
- **Phase 4 (US1)**: After Phase 3 (matcher needs loaded inputs). This is the MVP.
- **Phase 5 (US3)**: After Phase 4 (writer wires into `POLYMATCHApp::run()`).
- **Phase 6 (Polish)**: After all user stories.

### Cross-story dependency note

For this feature, US1 depends on US2's loaders (you cannot match without loaded data) and US3 wires into US1's app. The MVP requires all three user stories. Within each phase, the parallelism and TDD discipline still hold.

### Within-Phase Dependencies

- **Phase 2**: T007–T019 mostly [P] (different files); T008 depends on T007, T010 on T009, T015 on T014.
- **Phase 3**: All Phase-3 tests (T020–T033, T086, T087, T089) are [P] (one file, independent functions — written before any implementation per TDD). Implementation T034–T041: T035 depends on T034 (same file); T037 must precede T038/T039/T040 (same file); T041 depends on T034 & T036. T034 covers FR-018 hard-rejects (polygon ID 0, duplicate IDs).
- **Phase 4**: All Phase-4 tests (T042–T059, T088) are [P]. Implementation: T060→T061→T062 sequential (same file). T063→T064→T065→T066→T067 sequential (same file). T068 depends on T041 (US2) + T067 and covers FR-019 (skip zero/single-point trajectories). T069 depends on T068.
- **Phase 5**: All T070–T075 tests are [P]. T076→T077→T078 sequential (poly_mm_writer.cpp). T079 depends on T068 + T078.
- **Phase 6**: T080–T082 and T084 are [P]. T083 and T085 are sequential validations.

### Parallel Opportunities

- **Phase 2**: All 13 header-creation tasks (T007–T019) can run in parallel — different files.
- **Phase 3 tests** (T020–T033, T086, T087, T089): 17 tests in parallel inside `polymatch_test.cpp`.
- **Phase 3 implementations**: T034 + T036 + T037 are in different files — parallel. T035 depends on T034.
- **Phase 4 tests** (T042–T059, T088): 19 tests in parallel.
- **Phase 5 tests** (T070–T075): 6 tests in parallel.
- Across stories: A two-developer team can have one developer driving T034 (PolygonLayer) while the other drives T036 (AccessPointLayer), since they touch different files. T037–T040 require sequencing in the shared `polymatch_app_config.cpp`.

---

## Parallel Example: Phase 4 (US1) Tests

```bash
# Write all US1 tests first (they MUST FAIL before any implementation):
Task: "T042 [P] [US1] Add TEST(PolyLinkGraph, construction_vertex_and_arc_counts)"
Task: "T043 [P] [US1] Add TEST(PolyLinkGraph, through_cost_table_correctness)"
Task: "T046 [P] [US1] Add TEST(Emission, distance_zero_for_point_inside_polygon)"
Task: "T049 [P] [US1] Add TEST(MatchTraj, link_only_matches_weightmatch_exactly)"
Task: "T053 [P] [US1] Add TEST(MatchTraj, through_routing_no_gps_inside)"
# (and the rest of T042–T059 in parallel)
```

After all tests are red, implement T060–T069 in order.

---

## Implementation Strategy

### MVP scope

The MVP is **US2 + US1 + US3** combined — each story implements one third of the pipeline (load, match, output) and the polymatch executable is useless without all three. Strict story-by-story incrementalism would mean US2 can be demoed as "polymatch loads inputs and exits", US1 adds "polymatch produces in-memory matches", US3 adds "polymatch writes CSV output". This staging works for risk control but the production MVP requires all three.

### Incremental Delivery Path

1. Complete Phase 1 + Phase 2 (foundation).
2. Complete Phase 3 (US2): `polymatch --polygons X --access_points Y --network N --gps G` validates inputs end-to-end (no matching yet).
3. Complete Phase 4 (US1): `polymatch` produces in-memory matches; tested via `polymatch_test`. **Internal milestone — MVP candidate**.
4. Complete Phase 5 (US3): `polymatch` writes spec-compliant CSV output. **Production-ready release candidate**.
5. Complete Phase 6 (Polish): regression, performance, parallel-determinism gates → release.

### Single-Developer Sequencing

T001 → T002 → T003 → T004 → T005 → T006 → (T007..T019 in parallel where possible) → (T020..T033 + T086 + T087 + T089) tests → T034 → T035 → T036 → T037 → T038 → T039 → T040 → T041 → (T042..T059 + T088) tests → T060 → T061 → T062 → T063 → T064 → T065 → T066 → T067 → T068 → T069 → T070..T075 tests → T076 → T077 → T078 → T079 → T080 → T081 → T082 → T083 → T084 → T085.

---

## Notes

- **TDD is mandatory** per Constitution Principle III. Every test task (T020–T033, T042–T059, T070–T075, T086, T087, T088, T089) is written **before** the corresponding implementation and must compile + fail before the implementation task is started.
- **No edits to existing matcher source files** (FMM, STMATCH, WEIGHTMATCH, H3MM). Modifications to `CMakeLists.txt` and `test/CMakeLists.txt` are additive only.
- **No Python bindings** are added; `python/fmm.i` is not modified.
- **OpenMP parallelism is a v1 requirement** (FR-017); per-thread state created inside the parallel region, shared graph data treated as const, writer guarded by mutex.
- **Through-routing cost precomputed once** at `PolyLinkGraph` construction (R11); `THROUGH_PENALTY_FACTOR` applied at lookup so factor sweeps don't require a rebuild.
- Commit after each task or at story checkpoints. Don't skip the test-first ordering — Constitution Principle III is non-negotiable.
