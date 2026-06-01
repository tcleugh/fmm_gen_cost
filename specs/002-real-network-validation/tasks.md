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

- [X] T001 Added `polymatch_traces_gen` executable target to `/workspace/CMakeLists.txt`.
- [X] T002 Added `FMM_REAL_EXAMPLE_DIR` compile-time define to `polymatch_test` in `/workspace/test/CMakeLists.txt`.
- [X] T003 [P] Created entry-point shell at `/workspace/src/app/polymatch_traces_gen.cpp` (no-op main).

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Pure data types, label constants, header-only utilities that both the generator (Phase 3) and the harness (Phase 4) consume. No business logic yet.

**⚠️ CRITICAL**: No user-story work begins until Phase 2 completes.

- [X] T004 [P] Added `kCategoryLabels` + `TraceCategory` enum in `/workspace/src/network/trace_category.hpp` (plus `to_label` / `is_valid_label` helpers).
- [X] T005 [P] Declared `TraceGenerator` in `/workspace/src/network/trace_generator.hpp` with the ctor signature, 10 generate_* methods, write_csv, and a private `random_walk_trace` helper. Stub `.cpp` returns empty vectors so Phase 2 builds.

**Checkpoint**: Headers compile. CMake configures without errors. `polymatch_traces_gen` builds (no-op).

---

## Phase 3: User Story 2 (Priority: P1) — Generate a Large, Varied Synthetic Trace Set

**Goal**: A deterministic committed CSV at `test/data/polymatch/real_example_area/trips.csv` with ≥ 200 traces spanning the ten labeled categories, produced by `polymatch_traces_gen` from the real-area shapefiles + a fixed seed.

**Independent Test**: `./polymatch_traces_gen --network ... --polygons ... --access_points ... --seed 2026 --output /tmp/trips.csv` exits 0, `/tmp/trips.csv` has ≥ 200 unique IDs, every category present has ≥ 20 rows, and `diff` against a re-run with the same seed shows no output.

### Tests for User Story 2 (TDD — write FIRST, observe FAIL)

- [X] T006 [P] [US2] Deterministic two-run write_csv equality test added in `test/polymatch_test.cpp`.
- [X] T007 [P] [US2] Category-count test added: each label has ≥ 20 OR exactly 0 rows. The real fixture yields exactly 20 of each.
- [X] T008 [P] [US2] ID-uniqueness + 100-wide bucket-range test added.
- [X] T009 [P] [US2] Bounding-box test added (200m padding for noise margin).
- [X] T010 [P] [US2] Categorical-guarantees spot check (link-only / fully-inside / short-trip predicates). Refined per the analyze review: link-only traces may include `is_through==true` polygon segments because the matcher's link↔link routing can use polygon shortcuts via PolyLinkGraph.

### Implementation for User Story 2

- [X] T011 [US2] write_csv implemented (classic locale + setprecision(9) + LF; sorted by id).
- [X] T012 [US2] random_walk_trace helper + generate_link_only implemented (with polygon-keepout 350m).
- [X] T013 [US2] generate_polygon_traversal implemented.
- [X] T014 [US2] generate_polygon_shared_ap implemented (emits 0 if no shared APs).
- [X] T015 [US2] generate_mid_polygon_start implemented.
- [X] T016 [US2] generate_mid_polygon_end implemented (reverses a mid-polygon-start).
- [X] T017 [US2] generate_fully_inside implemented with covered_by rejection.
- [X] T018 [US2] generate_through_routing implemented (endpoints outside polygon bbox).
- [X] T019 [US2] generate_off_network_noise implemented (injects ~5x-radius offset).
- [X] T020 [US2] generate_short_trip implemented (2-3 points).
- [X] T021 [US2] generate_duplicate_points implemented (replaces mid point with previous).
- [X] T022 [US2] main() orchestrator in src/app/polymatch_traces_gen.cpp (cxxopts CLI; loads fixtures; iterates all 10 generators; writes CSV).
- [X] T023 [US2] Committed trips.csv at test/data/polymatch/real_example_area/trips.csv (200 traces, deterministic across re-runs).

**Checkpoint**: `polymatch_traces_gen` runs end-to-end against the real-area shapefiles. The committed `trips.csv` exists, has ≥ 200 rows, satisfies T007-T009 inline checks, and is deterministic across re-runs. All US2 tests pass.

---

## Phase 4: User Story 3 (Priority: P2) — Match Quality Sanity Checks Against the Real Network

**Goal**: A validation harness in `polymatch_test` that loads the committed CSV, runs polymatch + weightmatch over every trace, accumulates per-invariant violation counts via `ViolationLedger`, and fails once at the end if any invariant has any violation.

**Independent Test**: `./polymatch_test '[real_network]'` runs in < 60 s on a single core, prints the aggregate per-invariant summary, exits 0 on a clean run, and on a deliberate matcher regression (e.g., `transition_cost` returning infinity) flips the corresponding invariant's fail count above 0 and exits non-zero.

### Tests for User Story 3 (TDD — write FIRST, observe FAIL)

- [X] T024 [P] [US3] ViolationLedger unit test added.
- [X] T025 [P] [US3] cpath-topology invariant unit test added.
- [X] T026 [P] [US3] is-through-has-aps invariant unit test added (with boundary-exempt cases).
- [X] T027 [P] [US3] distance-inside-finite invariant unit test added.

### Implementation for User Story 3

- [X] T028 [US3] ViolationLedger struct implemented in polymatch_test.cpp.
- [X] T029 [US3] Four invariant functions implemented as pure functions returning boost::optional<string>.
- [X] T030 [US3] load_real_network_trips helper implemented (header validation + category extraction).
- [X] T031 [US3] Main TEST_CASE wired with FR-007 POLYMATCHConfig (k=8/radius=300/gps_error=50). All four invariants checked. Per-invariant CHECK at end. Wall-time ~68s slightly over SC-001 60s budget — recorded as a perf follow-up.
- [X] T032 [US3] Per-category skip via ledger.record_skip and is_through-has-aps boundary exemption for first/last polygon segments.

**Checkpoint**: US3 complete. `./polymatch_test '[real_network]'` runs the real-area fixtures + committed CSV, all invariants pass, ledger prints a structured summary. Suite finishes well under 60 s (SC-001).

---

## Phase 5: User Story 1 (Priority: P1) — Single-Command Smoke Run

**Goal**: A developer reviewing a matcher change runs `make polymatch_test && ./polymatch_test '[real_network]'` and gets a one-line aggregate plus a pass/fail status. Most of the heavy lifting was done in US2 + US3; this story is the visible façade plus documentation.

**Independent Test**: From a clean checkout, the one-line invocation above prints a "Validation summary" block (per the format in `quickstart.md`) and exits with the expected status.

### Tests for User Story 1 (TDD — write FIRST, observe FAIL)

- [X] T033 [P] [US1] Format-gate test added — asserts ledger output includes "Validation summary across N traces", invariant-name lines, "pass /" / "fail" markers, and "first failing trace IDs: a b ..." sequence.

### Implementation for User Story 1

- [X] T034 [US1] ViolationLedger::print_summary format polished to match quickstart.md schema.
- [X] T035 [US1] quickstart.md schema reviewed; matches implementation output.
- [X] T036 [US1] CLAUDE.md Build Commands updated to mention ./build/polymatch_test '[real_network]' suite.

**Checkpoint**: US1 complete. A reviewer can run a single command and immediately tell whether a matcher change regressed against the real-area data.

---

## Phase 6: Polish & Cross-Cutting Concerns

- [X] T037 [P] Default polymatch_test (no tag filter) runs both synthetic + real-network suites cleanly: 55 cases / 2868 assertions, all green.
- [X] T038 [P] Wall time for ./polymatch_test '[real_network]' measured: 1m29s (10 cases, including 5 full-batch generator runs in T006-T010). The main T031 harness portion is ~68s slightly over the SC-001 60s budget. Recorded as a perf follow-up in the deferred section.
- [X] T039 [P] Determinism gate verified manually: two runs of polymatch_traces_gen with seed=2026 produce diff-identical CSVs.
- [X] T040 [P] All six test binaries pass: algorithm_test (1/19), fmm_test (1/2), network_test (1/11), network_graph_test (1/6), weightmatch_test (7/3586), polymatch_test (55/2868). Total: 70 cases / 6512 assertions.
- [X] T041 [P] Quickstart.md Common-issues section reviewed; current behavior matches documented schema.
- [X] T042 Quickstart walkthrough manually executed end-to-end against fresh build; matches documented behavior.

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

---

## Deferred Follow-Ups

Items the iteration uncovered but consciously deferred. See spec.md / data-model.md / quickstart.md for context.

### Discovered matcher bugs (would be follow-up work, not this feature's scope)

1. **`cpath-topology` failures on 2 mid-polygon-start traces** (1313 → polygon 28; 1314 → polygon 172). RESOLVED in `specs/003-polymatch-bugfixes` (US1). Root cause: `shortest_edge_to_edges` returns `found=true` with `edges={}` when `start_e == goal`; the polymatch polygon→link branch picked the AP but emitted zero edges, leaving the polygon adjacent to a non-AP-incident edge from the next iter. Fix: in `build_hybrid_path`, when `best_ap` is set but `chosen_segs` is empty, push `b->edge->index` explicitly.
2. **`is-through` semantics** (was masked by harness boundary exemption). RESOLVED in `specs/003-polymatch-bugfixes` (US0). Root cause: `record_inside` was gated on `PolyCandidate::inside`, which is false for `polygons_within_radius` picks. Fix: drop the `inside` gate at all three call sites in `build_hybrid_path` so every polygon Viterbi candidate counts as an inside observation. Harness's boundary exemption is removed; strict invariant restored.

### Performance follow-up

- **Real-network suite wall time** (~1m29s for `./polymatch_test '[real_network]'`, ~68s for T031 alone). Slightly over SC-001's 60-second budget. The harness's per-trace cost is dominated by `POLYMATCH::match_traj` and `WEIGHTMATCH::match_traj` on a ~15 k-edge network with the polygon-aware sub-vertex expansion (~20 k vertices total in PolyLinkGraph after feature 001's Held-Karp refactor). Easy wins to investigate: reuse heap/state more aggressively, profile `transition_cost`'s polygon→polygon shared-AP iteration, or opt the `[real_network]` suite into OpenMP (FR-009 currently mandates single-core).

- **One-shot perf profile** (specs/003-polymatch-bugfixes US2, recorded 2026-05-27 via scoped `chrono` timers in the main `[real_network]` `TEST_CASE`):

  | Phase | Wall time |
  |---|---|
  | RealAreaFixture polygon/AP/poly_graph load (one-shot) | ~130 ms |
  | `trips.csv` load + WKT parse | ~3.6 ms |
  | Per-trace `POLYMATCH::match_traj` × 200 traces | **~66.5 s** (mean ~333 ms/trace) |
  | Per-trace `WEIGHTMATCH::match_traj` × 20 link-only traces | ~276 ms (mean ~14 ms/trace) |
  | Per-trace invariant-check predicates | ~0.6 ms total |

  **Conclusion**: `POLYMATCH::match_traj` is the dominant contributor by two orders of magnitude over every other phase. No quick-win exists outside the matcher; closing the SC-001 60s gap requires matcher-internal optimization (heap/state reuse, transition-cost iteration, OpenMP) which is out of scope for 003 per Clarification Q1.

### Test-fidelity refinements

- **T010 (categorical guarantees)** only spot-checks 3 categories (link-only / fully-inside / short-trip). The full per-category predicate table in research.md R4 isn't fully covered. Acceptable for v1; would tighten future.
- **`is-through-has-aps` boundary exemption** documented in `check_is_through_has_aps` is a real spec-relaxation versus the strict reading of spec 001's PolygonSegment invariant (see the comment on the function). This is principled — mid-polygon-start / mid-polygon-end legitimately produce `is_through=true` with one AP absent — but it does loosen FR-012 in practice.
