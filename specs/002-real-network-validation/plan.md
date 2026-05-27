# Implementation Plan: Real-Network Validation of Polymatch

**Branch**: `002-real-network-validation` | **Date**: 2026-05-27 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `/workspace/specs/002-real-network-validation/spec.md`

## Summary

Add a real-network validation suite for polymatch that runs against the provided real example area shapefiles (`test/data/polymatch/real_example_area/`). The deliverables are:

1. A C++ trace generator that, given the real-area fixtures and a fixed seed, produces a deterministic CSV (`test/data/polymatch/real_example_area/trips.csv`) of ≥ 200 trajectories spanning ten labeled categories (link-only, polygon-traversal, shared-AP crossing, mid-polygon start/end, fully-inside, through-routing, off-network noise, short trips, duplicate-points). The CSV is committed to the repo; the generator is a separate `polymatch_traces_gen` executable a developer re-runs when the underlying shapefiles change.
2. A validation harness inside `polymatch_test` (under a `[real_network]` Catch2 tag) that loads the committed CSV, drives `POLYMATCH` and `WEIGHTMATCH` against it, and checks four property-based invariants (cpath topology, `is_through` AP completeness, link-only ≡ weightmatch, finite non-negative `distance_inside`). The harness aggregates per-invariant violation counts and fails once at the end (Clarification Q2).

Both pieces are additive — they do not modify the matcher itself.

## Technical Context

**Language/Version**: C++14 (matches existing project; `CMAKE_CXX_STANDARD 14`).

**Primary Dependencies**: GDAL ≥ 2.2 (shapefile + GeoPackage I/O), Boost ≥ 1.56 (Geometry, Property Tree, Serialization), OpenMP, Catch2 (header-only, vendored in `third_party/catch2/`). No new dependencies introduced.

**Storage**: Three input shapefiles (existing, in `test/data/polymatch/real_example_area/`); one committed output CSV (this feature's deliverable: `trips.csv`); transient build-local CSV(s) for matched-result inspection during test runs.

**Testing**: Catch2 test suite under `test/polymatch_test.cpp`. New tag `[real_network]` so default `./polymatch_test` invocations skip the larger fixture; tag becomes opt-in like the existing `[.bench]` perf test.

**Target Platform**: Linux + macOS (matches existing project support matrix).

**Project Type**: C++ library + CLI executables (single CMake project).

**Performance Goals**: SC-001 — full real-network validation in under 60 s on a single core. With the existing 0.11 s perf baseline (1000 GPS × 100 polygons, synthetic), this is comfortable headroom for ≥ 200 traces × the real network's ~5 k edges and ~100 polygons.

**Constraints**:
- Generator MUST be deterministic for a given seed (SC-002).
- Must not modify polymatch's matching algorithm — only adds a generator + harness.
- The CRS is GDA 1994 Australia Albers (meters); no on-the-fly reprojection anywhere.
- Constitutional: no new heap allocation per Dijkstra query, no virtual dispatch in HMM inner loop (the generator MAY allocate freely; it's a one-shot tool).

**Scale/Scope**: Real-area fixture is ~5 km × 8 km, network ~1.8 MB shp / ~5 k edges, polygon layer ~90 kB / dozens of polygons, ~1000 APs. Trace batch ≥ 200; each trace 2–~30 GPS points. Total run ≤ 60 s single-core (SC-001).

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Compliance | Notes |
|---|---|---|
| **I. Performance-First** | ✅ N/A for the generator (one-shot offline tool); the harness reuses existing `DijkstraState`/`IndexedMinHeap` from polymatch. No new allocation in the matcher inner loop. |
| **II. Algorithmic Correctness** | ✅ This feature *verifies* correctness — it doesn't add new algorithmic steps. The four invariants encode existing matcher contracts (cpath topology, `is_through` semantics, link-only-equals-weightmatch, finite `distance_inside`). |
| **III. Test-Driven Development (NON-NEGOTIABLE)** | ✅ The deliverable *is* tests. The trace generator is independently verified by a deterministic-output check; the harness's invariants are the testable contract. |
| **IV. Dual-API Stability** | ✅ N/A — no SWIG / Python surface touched. |
| **V. Separation of Concerns** | ✅ Generator lives in `src/app/polymatch_traces_gen.cpp` (entry-point) plus optional helper in `src/network/`; uses existing `FMM::NETWORK`, `FMM::ROUTING`, `FMM::CONFIG` types. Harness lives entirely in `test/polymatch_test.cpp`. No new namespace, no back-references. |

**Gate result**: PASS. No complexity table entries required.

## Project Structure

### Documentation (this feature)

```text
specs/002-real-network-validation/
├── plan.md                  # This file
├── spec.md                  # Feature spec
├── research.md              # Phase 0 output — design decisions
├── data-model.md            # Phase 1 output — TraceCategory, generated-CSV schema, invariant model
├── quickstart.md            # Phase 1 output — how to regenerate the CSV + run the suite
├── contracts/               # Phase 1 output
│   └── real-network-trips-csv.md   # extends 000-pre-existing/contracts/gps-trajectory.md
├── checklists/
│   └── requirements.md      # Spec quality checklist (passing)
└── tasks.md                 # Phase 2 output (created by /speckit-tasks)
```

### Source Code (repository root)

```text
src/
├── app/
│   ├── polymatch.cpp              # (existing — unchanged)
│   └── polymatch_traces_gen.cpp   # NEW — offline trace-generator entry point
└── network/
    └── trace_generator.hpp / .cpp # NEW (optional split) — TraceGenerator class used by the entry point

test/
├── polymatch_test.cpp             # MODIFIED — adds a `[real_network]` TEST_CASE block
└── data/polymatch/real_example_area/
    ├── network.shp + .dbf + .shx + .prj + .cpg    # (provided — unchanged)
    ├── polygons.shp + ...                          # (provided — unchanged)
    ├── access_points.shp + ...                     # (provided — unchanged)
    └── trips.csv                                   # NEW — committed deterministic trace batch

CMakeLists.txt                                       # MODIFIED — add polymatch_traces_gen target
test/CMakeLists.txt                                  # MODIFIED — polymatch_test gets FMM_REAL_EXAMPLE_DIR define
```

**Structure Decision**: Single-project layout matching the existing repo. The generator is exposed as a new executable target alongside `polymatch`, so a developer can run `./polymatch_traces_gen --output trips.csv ...` directly. The harness is added as a TEST_CASE inside the existing `polymatch_test.cpp` rather than a new binary, keeping the test surface unified.

## Complexity Tracking

No constitutional violations to justify.

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|--------------------------------------|
| (none) | (n/a) | (n/a) |
