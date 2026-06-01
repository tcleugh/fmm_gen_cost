# Phase 0 Research: Real-Network Validation

This document captures the design decisions behind the trace generator and validation harness. Each section names what was chosen, why, and what alternatives were considered.

---

## R1. Trace generator: in-test C++ vs. separate executable vs. Python

**Decision**: A new C++ executable target `polymatch_traces_gen` (entry point `src/app/polymatch_traces_gen.cpp`). The trace-generation logic itself lives in a small `TraceGenerator` class under `src/network/` (or a feature-local file under `src/app/`), reused by the entry point.

**Rationale**:

- The generator needs `Network`, `PolygonLayer`, `AccessPointLayer`, and `PolyLinkGraph` to make matcher-dependent category decisions (through-routing requires pre-routing via `shortest_polylink_to_polylinks`; polygon-traversal needs a polygon-containment test). All four types are already in `FMM::NETWORK` / `FMM::ROUTING`.
- A separate executable keeps the test suite fast — the committed CSV is the artifact; the generator is only re-run when the input shapefiles change.
- C++ avoids adding a Python toolchain dep for a one-shot offline tool (and matches the convention set in `specs/001-polymatch-algorithm` for in-test fixtures).
- The Spec Clarifications mandate the CSV be committed; that decouples generator runtime from CI cycle time, which removes the "make this fast" pressure from the generator.

**Alternatives considered**:

- *In-test C++ generator that builds the CSV at `polymatch_test` startup* — Rejected. With a committed CSV (Clarification Q1) the test doesn't need to regenerate; doing so would slow every `polymatch_test` invocation by ~5–15 s for no benefit.
- *Python script (`test/generate_polymatch_real_traces.py`)* — Rejected. Reaching the same matcher-dependent decisions (especially through-routing endpoint selection) from Python would require either a SWIG binding for `PolyLinkGraph` (doesn't exist) or a CLI shell-out, both heavier than a 200-line C++ tool.

---

## R2. CSV column layout — extending `gps-trajectory.md` Format 1

**Decision**: Three semicolon-delimited columns `id;geom;category`. `id` is a stable int per trace; `geom` is a WKT `LINESTRING`; `category` is one of the ten labels in spec FR-004.

**Rationale**:

- `CSVTrajectoryReader` matches columns by header name (see `gps-trajectory.md`), so unknown columns are silently ignored. Polymatch consumes the file unchanged.
- The harness needs the category per-trace to pick which invariants apply (Clarification Q3). Storing it in the CSV is the cheapest correct option.
- Deterministic byte-output requires a stable column order; we lock the column order alphabetically: `category,geom,id` — no, **trajectory order** (`id,geom,category`) — to keep grep'ing for a trace ID convenient.

**Alternatives considered**:

- *Side car JSON file `trips_categories.json`* — Rejected. Two files double the artifact count and create a synchronization risk (CSV row 47's id doesn't map to the JSON's "47" via index — needs explicit join).
- *Encode category in the id (`{cat_prefix}{seq}`)* — Rejected. Breaks the int-id convention `polymatch` already enforces; harder to read; loses the ability to renumber.

---

## R3. Determinism strategy

**Decision**: The generator takes a seed (`--seed`, default `2026`), instantiates a `std::mt19937_64` seeded from it, and uses that engine for every random choice: starting edge, end edge, GPS noise, point count, category selection. The CSV is written with `std::ofstream` and an explicit floating-point precision (`std::setprecision(9)`) so platform-specific stream defaults don't cause diffs.

**Rationale**:

- `std::mt19937_64` is reproducible across platforms when seeded identically (it's deterministic per the C++ standard).
- `std::setprecision(9)` is enough to fully round-trip a `double` for typical coordinates in the GDA Albers projection (meter-scale; max range ~10⁶ → 9 significant digits covers full precision).
- Locale matters for `std::ostream` numeric formatting; we explicitly set `out.imbue(std::locale::classic())` to avoid `1,234.5` vs `1234.5` differences.

**Alternatives considered**:

- *Use Boost.Random for determinism* — Rejected. The standard library's `mt19937_64` is deterministic; adding Boost.Random is unnecessary.
- *Hex-encode coordinates to avoid round-trip issues* — Rejected. Overkill; the resulting CSV would be unreadable and diverge from `gps-trajectory.md` Format 1.

---

## R4. Category-selection algorithm

**Decision**: For each category we have a small constructor function in the generator that produces one trace. The generator's main loop calls each constructor `n_per_category` (= 20) times with seeded RNG, accumulating into a single output list. After all categories run, the list is sorted by `id` (id assigned by category, e.g., link-only ids 1000–1019, polygon-traversal 1100–1119, …) and written.

| Category | Constructor approach |
|---|---|
| `link-only` | Pick a random edge with no polygons within `radius`; pick a downstream edge via random walk over `LinkGraph`; sample GPS points along the polyline with Gaussian noise. |
| `polygon-traversal` | Pick a polygon, then pick an entry edge incident to one of its link-attached APs and an exit edge incident to a different AP. Sample points along edge→AP→inside-polygon→AP→edge. |
| `polygon-shared-ap` | Find two polygons sharing an AP; if none exist, emit 0 traces and the harness skips this category. Walk endpoint-edge → AP → polygon A → shared AP → polygon B → AP → endpoint-edge. |
| `mid-polygon-start` | Pick a polygon; sample first GPS point strictly inside; route out via an AP-incident edge. |
| `mid-polygon-end` | Symmetric of above. |
| `fully-inside` | Pick a polygon; sample all GPS points uniformly in its bounding box, rejection-test against `boost::geometry::covered_by`. |
| `through-routing` | Pick a polygon; pick two link-attached APs of it; pick external edges so the optimal route over `PolyLinkGraph` (with the configured `through_penalty_factor`) crosses the polygon, but the sampled GPS points stay outside the polygon's boundary (placed at offsets such that `polygons_within_radius` returns nothing or only outside-distance candidates). |
| `off-network-noise` | Pick a valid network trace, then inject 1–3 GPS points with large Gaussian offset (≥ `4 × radius`) into the middle. |
| `short-trip` | Pick a single edge; sample 2–3 GPS points. |
| `duplicate-points` | Pick a valid trace; replace 1–2 internal consecutive GPS points with the previous point's coordinates. |

**Rationale**: Per-category constructors are small (~30 LOC each), independently testable, and let us inject category-specific guarantees the matcher can't infer post-hoc (e.g., "through-routing implies no GPS inside the polygon" is enforced at generation, not validated by matcher inspection).

**Alternatives considered**:

- *Generic Markov-chain-style random walk with post-hoc category labeling* — Rejected per Clarification Q3 (Option A in the question; we chose Option C, generator-stamped labels).
- *Hand-craft a fixed set of traces* — Rejected as it doesn't meet SC-003 (≥ 200 traces with ≥ 20 per category) at any reasonable authoring cost.

---

## R5. Harness wiring — Catch2 tag, fixture sharing

**Decision**: A new `TEST_CASE("Real-network validation against committed trace batch (US1-US3)", "[polymatch][real_network]")` block in `test/polymatch_test.cpp`. It is **not** hidden behind a `[.bench]`-style "skipped by default" tag — the suite must run by default so polymatch_test catches regressions.

The harness:

1. Loads the real-area shapefiles via the existing `PolygonLayer` / `AccessPointLayer` / `Network` constructors. Path comes from a new compile-time define `FMM_REAL_EXAMPLE_DIR` injected via `target_compile_definitions` (same pattern as `POLYMATCH_FIXTURE_DIR`, established in `CLAUDE.md`'s path convention).
2. Loads `trips.csv` via `CSVTrajectoryReader` (the extra `category` column is silently ignored).
3. Re-reads the same CSV manually to extract the `category` column into a parallel `std::vector<std::string>`.
4. Constructs both `POLYMATCH` and `WEIGHTMATCH` matchers; runs each trace through both.
5. Tallies per-invariant pass/fail counts; at the end, asserts each invariant's fail count is zero. On failure, prints up to 10 failing trace IDs per invariant (Clarification Q2).

**Rationale**:

- Compile-time defines keep the path portable across mount points (Memory: `feedback_path_convention`).
- Catch2's `CHECK`/`INFO` macros let us accumulate violations without short-circuiting — perfect for "aggregate, fail at end."
- Not hiding behind `[.bench]` means CI catches regressions; the projected time budget (≤ 60 s, SC-001) is below the threshold where a CI gate becomes painful.

**Alternatives considered**:

- *Run the suite via a separate test binary `real_network_test`* — Rejected. Adds a CMake target for marginal value; the suite belongs with polymatch's other tests.
- *Tag as `[.real_network]` (hidden by default)* — Rejected. The whole point is to catch regressions in the default flow; hiding it defeats that.

---

## R6. Invariant verification implementation

**Decision**: Each invariant is a small free function in the harness anonymous namespace, taking a `PolyMatchResult` + the original `Trajectory` + the original `Network`/`PolygonLayer`/`AccessPointLayer` references. Returns `std::optional<std::string>` — `std::nullopt` on pass; a one-line failure description on fail.

The harness drives them through a `ViolationLedger` helper that maintains `std::vector<std::pair<int /*trace_id*/, std::string /*reason*/>>` per invariant.

**Rationale**:

- Pure-function invariants are easy to unit-test in isolation (constitution Principle III) and easy to extend (add a new invariant = add a new function).
- The `ViolationLedger` consolidates per-invariant printing and Catch2 assertion at end-of-test, avoiding boilerplate at every check site.

**Alternatives considered**:

- *Use Catch2's REQUIRE per trace per invariant* — Rejected. `REQUIRE` aborts the section on first failure; that's fail-fast, contradicting Clarification Q2.
- *External golden-file comparison* — Rejected. Per spec User Story 3, golden files for 200 real-network traces are prohibitive; property-based invariants are the explicit choice.

---

## R7. Handling categories the network can't support

**Decision**: When the generator can't construct a trace in a given category (e.g., the real fixture has no two polygons sharing an AP, so `polygon-shared-ap` is empty), it emits zero traces for that category. The CSV's `category` column reflects only categories with ≥ 1 trace. The harness counts per-category trace populations on load; for a category with zero traces, the corresponding invariant assertion is replaced with `INFO("category=X skipped — no traces produced")` and the test does not fail on it.

**Rationale**: Already pinned in spec Assumptions. Codified here for the harness implementation.

**Alternatives considered**:

- *Fail the suite when a category is empty* — Rejected. The real fixture is what it is; failing a test because the fixture's topology doesn't support a category penalizes good fixtures.

---

## Summary

| # | Topic | Decision |
|---|---|---|
| R1 | Generator implementation | New C++ executable `polymatch_traces_gen`; one-shot offline |
| R2 | CSV layout | `id;geom;category` — extends `gps-trajectory.md` Format 1 |
| R3 | Determinism | `std::mt19937_64` + explicit precision + `locale::classic()` |
| R4 | Category-selection | Per-category constructor functions, 20 traces each |
| R5 | Harness wiring | TEST_CASE in `polymatch_test.cpp` tagged `[real_network]` (default-on) |
| R6 | Invariant verification | Pure functions returning `optional<string>`; `ViolationLedger` for aggregation |
| R7 | Missing categories | Skip with INFO; don't fail |

All Phase 0 unknowns resolved. No `NEEDS CLARIFICATION` markers introduced.
