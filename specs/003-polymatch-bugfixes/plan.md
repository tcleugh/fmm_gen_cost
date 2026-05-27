# Implementation Plan: Polymatch Bug Fixes from Real-Network Validation

**Branch**: `003-polymatch-bugfixes` | **Date**: 2026-05-27 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `/workspace/specs/003-polymatch-bugfixes/spec.md`

## Summary

Two correctness fixes plus a performance diagnostic, all narrowly scoped to `src/mm/polymatch/polymatch_algorithm.cpp` plus the 002 harness assertions:

- **US0 (P1) — `is_through` semantics**: `POLYMATCH::build_hybrid_path` currently sets `has_inside_obs` only when a polygon Viterbi candidate's `inside` flag is `true` (covered_by-strict). The fix sets it whenever a polygon Viterbi candidate for that polygon is matched at any GPS layer, regardless of `inside`. The 002 harness's `check_is_through_has_aps` first/last boundary exemption is then removed.
- **US1 (P1) — cpath-topology violations on traces 1313, 1314**: The matcher emits an edge immediately after a polygon in cpath whose endpoints are not in that polygon's AP set, on two real-network mid-polygon-start traces. Root cause unknown; Phase 0 investigation identifies the buggy code path.
- **US2 (P2) — perf diagnostic**: Time the real-network suite, identify the dominant cost contributor in a 1-2 paragraph summary, record in 002's Deferred Follow-Ups. Apply *only* a single obvious quick-win optimization if one surfaces; defer larger work.

US0 and US1 are likely related — broken `is_through` semantics correlates with mid-polygon-start traces having unexpected polygon-segment shapes, which is exactly where the cpath-topology bug manifests. Phase 0 verifies or refutes the connection.

## Technical Context

**Language/Version**: C++14 (matches existing project; CMake CXX_STANDARD 14, constitution says C++11 but build uses C++14).

**Primary Dependencies**: GDAL ≥ 2.2, Boost ≥ 1.56 (Geometry, Property Tree), OpenMP, Catch2 (vendored). No new dependencies.

**Storage**: No new files. Existing real-area shapefiles + committed `trips.csv` (feature 002) are the regression input.

**Testing**: Catch2 — existing 70 cases across 6 binaries are the regression gate. The 002 [real_network] suite (`polymatch_test '[real_network]'`) is the primary signal for both fixes.

**Target Platform**: Linux + macOS.

**Project Type**: C++ library + CLI executables (single CMake project).

**Performance Goals**: Maintain feature 002's measured ~89s wall time as a *floor* (no perf regression). FR-005 / FR-006 produce a profile but the actual <60s budget restoration is deferred.

**Constraints**:
- The fix lives in `src/mm/polymatch/polymatch_algorithm.cpp` only. No changes to the generator (`src/network/trace_generator.cpp`), the loaders (`src/network/{polygon,access_point}_layer.cpp`), the routing graph (`src/network/poly_link_graph.cpp`), or the harness (`test/polymatch_test.cpp`) other than removing the boundary-exemption block and tightening the `cpath_fail <= 5` tolerance.
- `link-only-eq-weightmatch` is a strictly held invariant (Clarification Q2). A candidate fix that changes any link-only trace's matched cpath is rejected even if the new behavior is also valid.
- The committed `trips.csv` (T023 from 002) MUST NOT change.

**Scale/Scope**: Same as 002 — 15 k-edge real network, 241 polygons, 1.5 k access points, 200 traces in the committed batch.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Compliance | Notes |
|---|---|---|
| **I. Performance-First** | ✅ | The fix is a logic change in `build_hybrid_path` — no new heap allocation per Dijkstra query. US2 (perf diagnostic) explicitly preserves the constitutional non-allocation rule. |
| **II. Algorithmic Correctness** | ✅ | This feature IS Principle II in action — the matcher's `is_through` semantics and cpath topology were silently broken on real data; the fix restores correctness. No new HMM/Viterbi behavior; only the result-assembly phase changes. |
| **III. Test-Driven Development (NON-NEGOTIABLE)** | ✅ | Both fixes are gated by the 002 harness's existing invariants (`is-through-has-aps` strict, `cpath-topology` strict). The TDD discipline: tighten the harness's CHECKs first (so they fail), then fix the matcher until they pass. |
| **IV. Dual-API Stability** | ✅ N/A | No SWIG / Python surface touched. |
| **V. Separation of Concerns** | ✅ | Changes confined to `FMM::MM::POLYMATCH::build_hybrid_path` (and possibly `transition_cost`); no cross-namespace ripple. |

**Gate result**: PASS. No complexity-tracking entries required.

## Project Structure

### Documentation (this feature)

```text
specs/003-polymatch-bugfixes/
├── plan.md                # This file
├── spec.md                # Feature spec (with Clarifications log)
├── research.md            # Phase 0 — root cause analysis of traces 1313 + 1314 + is_through bug
├── data-model.md          # Phase 1 — small; only documents the corrected has_inside_obs semantics
├── quickstart.md          # Phase 1 — how to reproduce / verify / regress
├── checklists/
│   └── requirements.md    # Spec quality checklist (passing)
└── tasks.md               # Phase 2 output (created by /speckit-tasks)
```

No `contracts/` subdirectory — this feature changes no external interface. The matched-output CSV columns, polymatch CLI flags, and validation-harness invariants stay as-is (US0 and US1 fix the *implementation* against the existing contracts).

### Source Code (repository root)

```text
src/
└── mm/polymatch/
    └── polymatch_algorithm.cpp        MODIFIED — build_hybrid_path's record_inside trigger + the
                                        polygon→link edge-emission path for traces 1313/1314. Likely
                                        also a small change to the set_entry / set_egress / record_inside
                                        helpers. Estimated < 50 LOC of net change.

test/
└── polymatch_test.cpp                  MODIFIED — tighten check_is_through_has_aps (remove the
                                        boundary exemption block) + tighten the cpath_topology
                                        tolerance from `<= 5` to `== 0`. Per Clarification Q2:
                                        a handful of synthetic-fixture specific-value assertions
                                        may need updating; each touched assertion gets an inline
                                        comment naming the bug-fix commit.

specs/002-real-network-validation/
└── tasks.md                            MODIFIED — record the perf profile (FR-006) in Deferred
                                        Follow-Ups; update the SC-001 line in 002's tasks/spec
                                        if the profile reveals an obvious quick win that brings
                                        the wall time closer to 60s.
```

**Structure Decision**: Surgical edits only. No new files. No new CMake targets. No new namespace, no new entity types. The feature is intentionally tiny in code-volume terms (likely < 50 lines of matcher code changes + harness assertion tightening) but high in spec-correctness value.

## Complexity Tracking

No constitutional violations to justify.

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|--------------------------------------|
| (none) | (n/a) | (n/a) |
