---

description: "Task list for Real-Network Validation of Polymatch"
---

# Tasks: Real-Network Validation of Polymatch

**Input**: Design documents from `/workspace/specs/002-real-network-validation/`

**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/, quickstart.md

**Tests**: TDD is mandatory per Constitution Principle III (NON-NEGOTIABLE). Test tasks below MUST be written FIRST and observed to fail before their implementation counterparts begin.

**Organization**: Tasks are grouped by user story so each story can be implemented + verified independently. The MVP is **US2 + US3** (a deterministic CSV the harness can validate); US1 is the developer-facing single-command façade and lands last.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Different files / no dependencies on incomplete tasks.
- **[Story]**: Maps to user stories from spec.md (US1, US2, US3).
- All file paths are absolute under `/workspace/`.

## Path Conventions

Single C++ project. Generator source under `src/app/` + helper under `src/network/`; harness under `test/`; committed fixture artifact under `test/data/polymatch/real_example_area/`.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: New CMake targets + compile-time defines + an empty entry-point shell so the rest of the work can compile incrementally.

- [ ] T001 Add `polymatch_traces_gen` executable target to `/workspace/CMakeLists.txt`: `add_executable(polymatch_traces_gen src/app/polymatch_traces_gen.cpp)` plus `target_link_libraries(polymatch_traces_gen FMMLIB)` (mirror the `polymatch` target a few lines above).
- [ ] T002 Add `FMM_REAL_EXAMPLE_DIR` compile-time define to `/workspace/test/CMakeLists.txt`'s `polymatch_test` target: `target_compile_definitions(polymatch_test PRIVATE FMM_REAL_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/test/data/polymatch/real_example_area")`. Mirrors the existing `POLYMATCH_FIXTURE_DIR` / `FMM_TEST_DATA_DIR` pattern.
- [ ] T003 [P] Create entry-point shell at `/workspace/src/app/polymatch_traces_gen.cpp` containing just `int main(int argc, char**argv){ (void)argc; (void)argv; return 0; }` so T001's target builds cleanly before the generator logic lands.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Pure data types, label constants, header-only utilities that both the generator (Phase 3) and the harness (Phase 4) consume. No business logic yet.

**⚠️ CRITICAL**: No user-story work begins until Phase 2 completes.

- [ ] T004 [P] Add the `kCategoryLabels` constant array and a `TraceCategory` enum (matching the ten labels in spec FR-004) to a new feature-local header `/workspace/src/network/trace_category.hpp`. Used by both the generator and the harness so the label set is single-sourced.
- [ ] T005 [P] Declare the `TraceGenerator` class skeleton in `/workspace/src/network/trace_generator.hpp`: ctor takes `(const Network&, const PolygonLayer&, const AccessPointLayer&, const PolyLinkGraph&, uint64_t seed)`; declares one public method per category (`generate_link_only`, `generate_polygon_traversal`, …, ten total) returning a `std::vector<GeneratedTrace>`; declares `write_csv(const std::string& out_path, const std::vector<GeneratedTrace>& traces) const`. Empty bodies / no .cpp yet.

**Checkpoint**: Headers compile. CMake configures without errors. `polymatch_traces_gen` builds (no-op).

---

## Phase 3: User Story 2 (Priority: P1) — Generate a Large, Varied Synthetic Trace Set

**Goal**: A deterministic committed CSV at `test/data/polymatch/real_example_area/trips.csv` with ≥ 200 traces spanning the ten labeled categories, produced by `polymatch_traces_gen` from the real-area shapefiles + a fixed seed.

**Independent Test**: `./polymatch_traces_gen --network ... --polygons ... --access_points ... --seed 2026 --output /tmp/trips.csv` exits 0, `/tmp/trips.csv` has ≥ 200 unique IDs, every category present has ≥ 20 rows, and `diff` against a re-run with the same seed shows no output.

### Tests for User Story 2 (TDD — write FIRST, observe FAIL)

- [ ] T006 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST_CASE("TraceGenerator emits ten labeled categories deterministically (T006)", "[real_network][us2]")`: construct generator with seed `2026` over the real-area fixtures, call `write_csv()` twice to two different temp paths, assert byte-equality via file-hash check (SC-002 / FR-002).
- [ ] T007 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST_CASE("TraceGenerator covers each category ≥ 20 traces OR exactly 0 (T007)", "[real_network][us2]")`: count rows per `category` column; assert each value in `kCategoryLabels` appears either ≥ 20 times or 0 times (R7).
- [ ] T008 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST_CASE("TraceGenerator IDs unique within batch (T008)", "[real_network][us2]")`: parse the generated CSV, assert no `id` collisions and that each id falls inside its category's 100-wide bucket (contract: `real-network-trips-csv.md`).
- [ ] T009 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST_CASE("TraceGenerator coordinates inside network bounding box (T009)", "[real_network][us2]")`: load the network shapefile, compute its bbox, assert every non-`off-network-noise` trace's coordinates lie inside; assert `off-network-noise` traces' off-network points lie outside by > `4 × radius` from any edge (FR-006).
- [ ] T010 [P] [US2] In `/workspace/test/polymatch_test.cpp`, add `TEST_CASE("TraceGenerator categorical guarantees match category label (T010)", "[real_network][us2]")`: for each emitted trace, run polymatch and verify the matched output exhibits the labeled category's expected shape. Implement one small predicate per category (parallel to the four invariant functions T029 introduces for matcher correctness — see C5 below for shared infrastructure). Predicates: `link-only` → `result.polygon_segments.empty()`; `polygon-traversal` → ≥ 1 polygon segment with `is_through == false` and both APs populated; `polygon-shared-ap` → ≥ 2 consecutive polygon segments where `seg[i].egress_ap == seg[i+1].entry_ap`; `mid-polygon-start` → first polygon segment has `entry_ap == kNoAccessPoint`; `mid-polygon-end` → last polygon segment has `egress_ap == kNoAccessPoint`; `fully-inside` → exactly one polygon segment with both APs `kNoAccessPoint`; `through-routing` → ≥ 1 polygon segment with `is_through == true`; `off-network-noise` → matched without exception (no further shape check); `short-trip` → matched without exception OR skipped per FR-019; `duplicate-points` → matched without exception, no NaN in any output value. Each predicate returns `std::optional<std::string>`. Skipped per-category when the generator emitted zero traces for that category.

  *(C5: T010 and T031 both iterate polymatch over every trace. Share a `for (auto& [traj, cat] : load_real_network_trips(...))` helper introduced in T030 so both call sites consume the same loop.)*

### Implementation for User Story 2

- [ ] T011 [US2] Implement `TraceGenerator::write_csv()` in `/workspace/src/network/trace_generator.cpp`: open `ofstream`, set `imbue(std::locale::classic())` + `std::setprecision(9)`, write header `id;geom;category\n`, then one row per trace sorted by `id` ascending (contract R3 + `real-network-trips-csv.md`). LF line endings only.
- [ ] T012 [US2] Implement the random-walk helper and `TraceGenerator::generate_link_only()` in `/workspace/src/network/trace_generator.cpp`: pick a start edge whose midpoint has no polygon within `radius`, random-walk over `LinkGraph` for 3-8 hops, sample 4-12 GPS points with Gaussian noise scaled to ~`radius/3`, assign IDs 1000-1099 (depends on T011).
- [ ] T013 [US2] Implement `TraceGenerator::generate_polygon_traversal()` in `/workspace/src/network/trace_generator.cpp`: pick a polygon with ≥ 2 link-attached APs, pick incoming/outgoing edges via those APs, sample GPS points across edge→inside-polygon→edge, IDs 1100-1199.
- [ ] T014 [US2] Implement `TraceGenerator::generate_polygon_shared_ap()` in `/workspace/src/network/trace_generator.cpp`: scan `AccessPointLayer` for APs with `polygons.size() >= 2`; if none exist emit zero traces and return empty vector; otherwise pick two polygons sharing an AP and walk endpoint→polygon A→shared AP→polygon B→endpoint, IDs 1200-1299.
- [ ] T015 [US2] Implement `TraceGenerator::generate_mid_polygon_start()` in `/workspace/src/network/trace_generator.cpp`: pick a polygon, sample first GPS point strictly inside via rejection on `polygons_containing`, route out via an AP-incident edge, IDs 1300-1399.
- [ ] T016 [US2] Implement `TraceGenerator::generate_mid_polygon_end()` symmetrically in `/workspace/src/network/trace_generator.cpp`, IDs 1400-1499.
- [ ] T017 [US2] Implement `TraceGenerator::generate_fully_inside()` in `/workspace/src/network/trace_generator.cpp`: pick a polygon, sample 4-10 GPS points uniformly inside its bounding box with rejection-test against `polygons_containing`, IDs 1500-1599.
- [ ] T018 [US2] Implement `TraceGenerator::generate_through_routing()` in `/workspace/src/network/trace_generator.cpp`: pick a polygon with ≥ 2 link-attached APs; pick external edges so a `shortest_polylink_to_polylinks` query crosses the polygon at `through_penalty_factor=0.5`; place GPS points only OUTSIDE the polygon's `radius`-expanded bounding box so no observation falls inside, IDs 1600-1699 (depends on T011, T012 for the random-walk helper).
- [ ] T019 [US2] Implement `TraceGenerator::generate_off_network_noise()` in `/workspace/src/network/trace_generator.cpp`: build a `link-only` base trace, inject 1-3 interior GPS points offset by 4-10 × `radius` perpendicular to the polyline, IDs 1700-1799 (depends on T012).
- [ ] T020 [US2] Implement `TraceGenerator::generate_short_trip()` in `/workspace/src/network/trace_generator.cpp`: pick a single edge, sample 2-3 GPS points along it with light noise, IDs 1800-1899.
- [ ] T021 [US2] Implement `TraceGenerator::generate_duplicate_points()` in `/workspace/src/network/trace_generator.cpp`: build a base `link-only` trace of 5-8 points, replace one consecutive pair with the previous point's exact coordinates, IDs 1900-1999 (depends on T012).
- [ ] T022 [US2] Implement `main()` orchestrator in `/workspace/src/app/polymatch_traces_gen.cpp`: parse CLI flags via `cxxopts` (`--network`, `--polygons`, `--access_points`, `--seed`, `--output`, default seed = `2026`); load `Network`, `PolygonLayer`, `AccessPointLayer`, `LinkGraph`, `PolyLinkGraph(through_penalty_factor=0.5)`; instantiate `TraceGenerator`; call each `generate_*` method with `n_per_category=20`; concatenate, sort by id, write CSV (depends on T011-T021).
- [ ] T023 [US2] Run the generator to produce `/workspace/test/data/polymatch/real_example_area/trips.csv` with the canonical seed `2026`. Commit the resulting CSV. (Manual step — documented in `quickstart.md`. The CSV is the artifact this story delivers.)

**Checkpoint**: `polymatch_traces_gen` runs end-to-end against the real-area shapefiles. The committed `trips.csv` exists, has ≥ 200 rows, satisfies T007-T009 inline checks, and is deterministic across re-runs. All US2 tests pass.

---

## Phase 4: User Story 3 (Priority: P2) — Match Quality Sanity Checks Against the Real Network

**Goal**: A validation harness in `polymatch_test` that loads the committed CSV, runs polymatch + weightmatch over every trace, accumulates per-invariant violation counts via `ViolationLedger`, and fails once at the end if any invariant has any violation.

**Independent Test**: `./polymatch_test '[real_network]'` runs in < 60 s on a single core, prints the aggregate per-invariant summary, exits 0 on a clean run, and on a deliberate matcher regression (e.g., `transition_cost` returning infinity) flips the corresponding invariant's fail count above 0 and exits non-zero.

### Tests for User Story 3 (TDD — write FIRST, observe FAIL)

- [ ] T024 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add a unit test for `ViolationLedger::record_pass / record_fail / any_failures` in isolation: 0 violations → `any_failures() == false`; ≥ 1 violation → `any_failures() == true`; `print_summary` shows the recorded pass/fail counts and up to 10 trace IDs.
- [ ] T025 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add a unit test that the `cpath-topology` invariant function returns `nullopt` for a hand-crafted valid `PolyMatchResult` (single link edge, no polygons) and returns a non-empty string for a `PolyMatchResult` whose `cpath` has two consecutive link IDs whose endpoint nodes don't connect.
- [ ] T026 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add a unit test that the `is-through-has-aps` invariant returns `nullopt` for a hand-crafted `PolygonSegment{ is_through=true, entry_ap=5, egress_ap=9 }` and returns a non-empty string for `{ is_through=true, entry_ap=kNoAccessPoint, egress_ap=9 }`.
- [ ] T027 [P] [US3] In `/workspace/test/polymatch_test.cpp`, add a unit test that the `distance-inside-finite` invariant returns `nullopt` for finite non-negative `distance_inside` and returns a non-empty string for `-1.0` and for `nan`/`inf` values.

### Implementation for User Story 3

- [ ] T028 [US3] Add a `ViolationLedger` struct (per data-model.md) inside the `polymatch_test.cpp` anonymous namespace with `record_pass`, `record_fail`, `print_summary`, `any_failures` methods. ≤ 80 LOC.
- [ ] T029 [US3] Implement the four invariant functions (`check_cpath_topology`, `check_is_through_has_aps`, `check_link_only_eq_weightmatch`, `check_distance_inside_finite`) in the `polymatch_test.cpp` anonymous namespace per data-model.md. Each is pure: takes the relevant references + a category string, returns `std::optional<std::string>`.
- [ ] T030 [US3] Implement a CSV-loader helper `load_real_network_trips(path) -> std::vector<std::pair<Trajectory, std::string /*category*/>>` in `polymatch_test.cpp`: re-parses the file row-by-row to extract both the geometry and the `category` column (the standard `CSVTrajectoryReader` ignores unknown columns; we need the category). Validates header == `id;geom;category`, IDs unique, categories in `kCategoryLabels` (per `real-network-trips-csv.md`). The returned vector is the iteration surface shared by T010 (categorical-guarantees) and T031 (invariant harness) — both TEST_CASEs call this helper to avoid duplicating the parse + per-trace polymatch invocation.
- [ ] T031 [US3] Add the main `TEST_CASE("Real-network validation against committed trace batch (US1-US3)", "[polymatch][real_network]")` to `/workspace/test/polymatch_test.cpp`. Note the tag is bare `[real_network]` (no leading dot) per FR-016 — runs by default. Load the three shapefiles via `FMM_REAL_EXAMPLE_DIR`; load `trips.csv` via T030; construct `Network`, `LinkGraph`, `PolygonLayer`, `AccessPointLayer`, `PolyLinkGraph`, `POLYMATCH`, `WEIGHTMATCH`; explicitly set POLYMATCHConfig per FR-007 (`k=8`, `radius=300`, `gps_error=50`, `boundary_epsilon=1e-6`, `through_penalty_factor=1.5`); iterate every trace via the shared loop helper from T030; per trace run polymatch (and weightmatch for `category == "link-only"` traces), feed results through each of the four invariants from T029 into a `ViolationLedger`; at end-of-test call `ledger.print_summary(std::cerr)` then `CHECK(!ledger.any_failures())` (depends on T028, T029, T030, T023).
- [ ] T032 [US3] In the same TEST_CASE, add per-category-skip handling per R7: tally per-category trace counts on load; for the `link-only-eq-weightmatch` invariant, only apply to traces with `category == "link-only"`; if a category has zero traces, the harness emits `INFO("category=X skipped — no traces produced")` instead of failing (depends on T031).

**Checkpoint**: US3 complete. `./polymatch_test '[real_network]'` runs the real-area fixtures + committed CSV, all invariants pass, ledger prints a structured summary. Suite finishes well under 60 s (SC-001).

---

## Phase 5: User Story 1 (Priority: P1) — Single-Command Smoke Run

**Goal**: A developer reviewing a matcher change runs `make polymatch_test && ./polymatch_test '[real_network]'` and gets a one-line aggregate plus a pass/fail status. Most of the heavy lifting was done in US2 + US3; this story is the visible façade plus documentation.

**Independent Test**: From a clean checkout, the one-line invocation above prints a "Validation summary" block (per the format in `quickstart.md`) and exits with the expected status.

### Tests for User Story 1 (TDD — write FIRST, observe FAIL)

- [ ] T033 [P] [US1] Add an integration assertion to the US3 TEST_CASE that verifies the printed summary lines match the documented `quickstart.md` format: one line per invariant ID, format `"<invariant-id> : N pass / M fail (first failing trace IDs: ...)"`. Use a `std::ostringstream` capture (override `print_summary`'s stream) so the test can grep the output string.

### Implementation for User Story 1

- [ ] T034 [US1] Polish `ViolationLedger::print_summary` formatting in `/workspace/test/polymatch_test.cpp` so the printed lines exactly match the `quickstart.md` schema (depends on T028, T033).
- [ ] T035 [US1] Update `/workspace/specs/002-real-network-validation/quickstart.md` if the actual implemented output drifts from the documented schema (sanity check; document-only).
- [ ] T036 [US1] Update `/workspace/CLAUDE.md` "Test commands" examples (Build Commands section) to mention `./polymatch_test '[real_network]'` as the real-network suite.

**Checkpoint**: US1 complete. A reviewer can run a single command and immediately tell whether a matcher change regressed against the real-area data.

---

## Phase 6: Polish & Cross-Cutting Concerns

- [ ] T037 [P] Run `./polymatch_test` (default invocation, no tag filter) and confirm both the synthetic and real-network suites pass in a single run (Constitution III).
- [ ] T038 [P] Measure single-core wall time for `./polymatch_test '[real_network]'` and assert (manually or via a `[.bench]` perf case) it is < 60 s (SC-001). If it exceeds 30 s, profile the harness — the matcher work is unchanged from 001, so any cost is in the harness loop.
- [ ] T039 [P] Determinism gate: re-run `polymatch_traces_gen --seed 2026` against an unchanged fixture, `diff` against the committed CSV; expect no output (SC-002). Document the command in `quickstart.md` (already present).
- [ ] T040 [P] Run `make tests && ./algorithm_test && ./fmm_test && ./network_test && ./network_graph_test && ./weightmatch_test && ./polymatch_test` from `build/test/`; confirm all six binaries pass cleanly (Constitution III, SC-004 of 001).
- [ ] T041 [P] Add a one-line entry to `/workspace/specs/002-real-network-validation/quickstart.md` "Common issues" pointing developers at this tasks.md when they need to extend the suite.
- [ ] T042 Manually execute `/workspace/specs/002-real-network-validation/quickstart.md` end-to-end against the freshly built executables; verify every step matches the documented behavior (no command typos, no missing flags, output matches the example block).

---

## Dependencies & Execution Order

### Phase dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately.
- **Phase 2 (Foundational)**: After Phase 1 (needs CMake targets and the compile-time define). Blocks both user stories.
- **Phase 3 (US2 — generator)**: After Phase 2. Produces the committed CSV that US3's harness consumes. **Internal milestone — generator-MVP**.
- **Phase 4 (US3 — harness)**: After Phase 3 (the harness needs the CSV from T023). **Internal milestone — validation-MVP**.
- **Phase 5 (US1 — single-command façade)**: After Phase 4 (the façade just polishes the existing US3 summary output). The smallest of the three stories.
- **Phase 6 (Polish)**: After all user stories.

### Cross-story dependency note

This feature's user stories are inverted from the typical priority-equals-leaf pattern: US1's command-line UX depends on US3's harness which depends on US2's CSV. The spec prioritizes US1 + US2 = P1 and US3 = P2 by *value* delivered, but the implementation order is US2 → US3 → US1 by *dependency*.

### Within-phase dependencies

- **Phase 2**: T004 and T005 are [P] — different files.
- **Phase 3 tests** (T006-T010): all [P], same file but independent TEST_CASE blocks.
- **Phase 3 impl**: T011 first (everything else uses `write_csv`); T012 (link-only) is the helper-builder needed by T018 (through-routing) and T019 (off-network-noise) and T021 (duplicate-points); other category generators (T013-T017, T020) are [P] amongst themselves once T011 lands. T022 main() depends on all of T011-T021. T023 (the committed CSV) is the final manual artifact step.
- **Phase 4 tests** (T024-T027): all [P], independent TEST_CASE blocks.
- **Phase 4 impl**: T028 → T029 → T030 in order (each follow-on uses the previous). T031 needs all three plus T023. T032 depends on T031.
- **Phase 5**: T033 (test) before T034 (impl); T035 and T036 are doc-only, parallel with everything else.
- **Phase 6**: T037, T038, T039, T040, T041 are [P]; T042 is the final manual validation.

### Parallel opportunities

- **Phase 2**: 2 header tasks in parallel (T004, T005).
- **Phase 3 tests**: 5 TEST_CASE blocks in parallel (T006-T010).
- **Phase 3 impl**: 6 category-generators in parallel once T011 (writer) and T012 (link-only / random-walk helper) land — T013, T014, T015, T016, T017, T020.
- **Phase 4 tests**: 4 unit tests in parallel (T024-T027).
- **Phase 6**: 5 polish tasks in parallel (T037-T041).

---

## Implementation Strategy

### MVP scope

The MVP is **US2 + US3 combined** — a deterministic CSV plus an invariant-checking harness. US1 (single-command façade) is a small finishing touch on top.

### Incremental delivery path

1. Phase 1 + Phase 2: scaffolding (executable + headers compile).
2. Phase 3 (US2): the generator runs and produces a committed `trips.csv` deterministically. **Internal milestone — generator-MVP** (committable).
3. Phase 4 (US3): the harness validates the CSV against four invariants and aggregates failures. **Internal milestone — validation-MVP**.
4. Phase 5 (US1): polish summary formatting + CLAUDE.md entry. Suite is usable by reviewers.
5. Phase 6: regression + determinism + full-suite gates.

### Single-developer sequencing

T001 → T002 → T003 → (T004, T005 parallel) → (T006-T010 parallel tests) → T011 → T012 → (T013, T014, T015, T016, T017, T020 parallel) → T018 → T019 → T021 → T022 → T023 (run + commit CSV) → (T024-T027 parallel tests) → T028 → T029 → T030 → T031 → T032 → T033 → T034 → T035 → T036 → (T037-T041 parallel polish) → T042.

---

## Notes

- **TDD is mandatory** per Constitution Principle III. Every test task (T006-T010, T024-T027, T033) is written and observed to fail BEFORE the corresponding implementation task starts.
- **The matcher is not modified.** This feature only adds a generator + harness; `polymatch_algorithm.{hpp,cpp}` and friends stay untouched. Any change to the matcher's source as part of this feature is a constraint violation.
- **CSV determinism** depends on `std::mt19937_64` + `std::setprecision(9)` + `std::locale::classic()`. Adding randomized OpenMP work to the generator would break this — keep it single-threaded.
- **Categories the network can't host** are silently skipped per R7 — the suite does not fail on missing trace counts. Document which categories were skipped in the harness's printed summary.
- Commit after each task or at story checkpoints. The generator commit (T022) + the CSV commit (T023) MUST be separate commits — the CSV is an artifact, the code that produces it is reviewable independently.
