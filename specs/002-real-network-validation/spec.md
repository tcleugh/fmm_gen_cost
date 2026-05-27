# Feature Specification: Real-Network Validation of Polymatch

**Feature Branch**: `002-real-network-validation`

**Created**: 2026-05-27

**Status**: Draft

**Input**: User description: "I have provided real example test data for polymatch test/data/polymatch/real_example_area Create a series of test traces for that region (large number of varied data) and verify it is running correctly on a real network"

## Clarifications

### Session 2026-05-27

- Q: Where does the generated trace CSV live in the dev workflow? → A: Committed artifact. The CSV is committed to the repo at `test/data/polymatch/real_example_area/trips.csv`; the trace generator is an offline tool a developer re-runs and commits when the fixture shapefiles change. CI reads the committed CSV directly.
- Q: When an invariant violation is detected, what does the suite do? → A: Aggregate, fail at end. The suite runs every trace, counts per-invariant violations, prints up to 10 failing trace IDs per invariant, then fails the test once at the end if any invariant has any violations.
- Q: How does each trace acquire its category label (especially matcher-dependent categories like through-routing)? → A: Generator stamps a `category` column. The generator picks each trace's category intentionally (e.g., pre-runs Dijkstra over the PolyLinkGraph to find endpoints whose route crosses a polygon) and writes the label as a third CSV column `id;geom;category`. The harness consumes the label to choose which invariants to apply per trace. `CSVTrajectoryReader` silently ignores unknown columns, so polymatch's input parsing is unaffected.

### Post-analyze refinements (also Session 2026-05-27)

These tighten wording in response to a `/speckit-analyze` cross-artifact review:

- FR-016 — the dedicated Catch2 tag is now named explicitly (`[real_network]`) and the spec confirms it is a default-on tag (no leading dot). The tag exists to give developers a filter, not to skip the suite by default. Resolves analyzer findings C1 + C2.
- FR-007 — the POLYMATCHConfig values the harness must set (`k=8`, `radius=300`, `gps_error=50`, `boundary_epsilon=1e-6`, `through_penalty_factor=1.5`) are now pinned in the spec rather than implied via "the defaults are fine." Resolves C3.
- FR-009 — broadened from "matcher invocation under 60 s" to "the full validation under 60 s," matching SC-001's scope. Resolves C6.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Smoke-Run Polymatch on the Real Example Network (Priority: P1)

A developer modifying the polymatch matcher needs a fast, repeatable way to confirm their change still produces sensible matches on a real road network — not just on the 12-edge synthetic fixture the unit tests use today. They drop the provided real-area shapefiles into the repo, run a single command, and get a pass/fail signal plus a sortable summary of any anomalies.

**Why this priority**: Without a real-network signal, polymatch's existing tests only prove the algorithm handles the 12-edge synthetic fixture. A 30-minute reviewer cannot tell whether a change breaks behavior on production-scale data. This story is the smallest unit of value that closes that gap.

**Independent Test**: From a clean checkout, `make polymatch_test && ./polymatch_test '[real_network]'` (or an equivalent invocation that runs only the real-network suite) terminates in under 60 seconds with exit code 0 and prints a one-line summary of how many traces matched, how many were skipped, and how many hit any defined sanity-check failures.

**Acceptance Scenarios**:

1. **Given** the real-network fixtures (network/polygon/access-point shapefiles in `test/data/polymatch/real_example_area/`) are present, **When** the developer runs the new real-network test, **Then** every trace runs to completion (no crash, no exception) and the test reports an aggregate pass/fail.
2. **Given** a developer has just modified `POLYMATCH::transition_cost`, **When** they re-run the suite, **Then** any trace whose matched output now violates a documented sanity invariant (see User Story 3) is identified in the test output by trace ID and the specific invariant violated.

---

### User Story 2 — Generate a Large, Varied Synthetic Trace Set (Priority: P1)

A test-data engineer needs to produce a representative batch of GPS traces that exercises the real network — straight road runs, mid-polygon trips, polygon entry/egress, polygon-crossing routes, mixed link+polygon trips, and edge cases (very short trips, very long trips, off-network noise, duplicate points). The traces must be generated reproducibly so anyone can regenerate identical data, and the resulting batch must be large enough (≥ 200 traces) and varied enough that a single broken cost rule will show up as multiple regressions rather than one.

**Why this priority**: A test suite that runs polymatch is only useful in proportion to the diversity of its input. A handful of hand-crafted traces lets bugs hide; a large, varied, deterministic set catches them. Without this story the suite from US1 has nothing meaningful to run.

**Independent Test**: Running the trace-generation tool against the real-area fixtures produces a CSV with ≥ 200 unique trip IDs covering each of the documented trace categories at least 20 times, deterministically across runs given the same seed.

**Acceptance Scenarios**:

1. **Given** the real-network fixtures are present and the trace generator is run with a fixed seed, **When** the generator completes, **Then** the output CSV contains at least 200 distinct trajectories whose IDs are stable across re-runs of the same seed.
2. **Given** two runs with different seeds, **When** comparing their outputs, **Then** the two CSVs differ in trace coordinates but match in row count and category distribution (within ± 5%).
3. **Given** the output CSV is loaded by `polymatch`, **When** matching completes, **Then** every documented category (link-only, polygon-traversal, polygon-shared-AP crossing, mid-polygon-start, mid-polygon-end, fully-inside, through-routing, off-network noise, short-trip, duplicate-points) is represented in the matched output by at least one trace whose result exhibits the expected shape.

---

### User Story 3 — Match Quality Sanity Checks Against the Real Network (Priority: P2)

A code reviewer needs confidence that the matcher is not just producing *some* output but is producing *plausible* output across the trace set. They want every emitted match to satisfy a small list of invariants that any correct polymatch result must hold — independent of any specific expected output — so regressions show up loudly even without a hand-curated golden file.

**Why this priority**: Hand-authoring expected outputs for ≥ 200 real-network traces is prohibitive. Property-based invariants (see Edge Cases + Functional Requirements below) catch the majority of correctness regressions at a tiny authoring cost. This story turns the matcher output into a verifiable signal without requiring per-trace ground truth.

**Independent Test**: A reviewer can run the same `polymatch_test '[real_network]'` from US1 and see, for each invariant defined in the Functional Requirements, an aggregate pass/fail count plus the first ≤ 10 violations per invariant printed for inspection.

**Acceptance Scenarios**:

1. **Given** the matched-output set for the trace batch, **When** the suite checks the cpath-topology invariant (every consecutive link-link pair shares a network node; every link-polygon transition uses one of that polygon's APs), **Then** the suite either reports zero violations or — if any violation is found — keeps running through every trace, prints up to 10 failing trace IDs at the end, and fails the test.
2. **Given** the matched-output set, **When** the suite checks that every PolygonSegment with `is_through == true` has both `entry_ap` and `egress_ap` set, **Then** the suite passes.
3. **Given** the matched-output set, **When** the suite checks that link-only mode applied to the same trace set produces output bit-identical to `weightmatch` on the same input, **Then** the suite passes (SC-002 from polymatch spec carried forward to the real network).

### Edge Cases

- **Empty / single-point trajectories**: Generator includes a small number of these; matcher MUST skip them with a per-trajectory warning (already enforced by polymatch FR-019).
- **Off-network noise**: Traces whose GPS coords sit far from any network edge or polygon — matcher MUST return either an empty match or a truncated match, never crash.
- **Duplicate consecutive GPS points**: Generator emits some traces with stuttered points; matcher MUST handle them without producing NaN/inf in transition probabilities.
- **Traces wholly outside the polygon layer**: A link-only path through the network, no polygon candidates — the matcher's output for these MUST be byte-identical to what `weightmatch` produces on the same trace.
- **Traces that cross the bounding-box edge of the fixture data**: GPS points outside the loaded network's bounding box — generator MUST avoid these by construction, or the matcher MUST gracefully truncate / skip without crashing.
- **Reproducibility**: With a fixed RNG seed, the generated trace CSV MUST be byte-identical across runs (so CI failures are debuggable).

## Requirements *(mandatory)*

### Functional Requirements

#### Trace generation

- **FR-001**: The trace generator MUST accept a real-example fixture directory (containing `network.shp`, `polygons.shp`, `access_points.shp`) and a random seed; it MUST emit a single CSV at `test/data/polymatch/real_example_area/trips.csv` extending the format documented in [`specs/000-pre-existing/contracts/gps-trajectory.md`](../000-pre-existing/contracts/gps-trajectory.md) (Format 1 — CSV trajectory, semicolon-delimited) with an added `category` column. The columns are `id;geom;category`. `CSVTrajectoryReader` silently ignores unknown columns by name, so polymatch consumes the file unchanged. The CSV is committed to the repository; the harness reads the committed file directly.
- **FR-002**: With the same seed, the generator MUST produce a byte-identical CSV across runs. This guarantees that a developer who re-generates after a fixture change produces a deterministic diff in git.
- **FR-003**: The generator MUST emit at least 200 distinct trace IDs.
- **FR-004**: The generator MUST cover each of the following categories with at least 20 traces each. The category label is written into the CSV's `category` column at generation time; for matcher-dependent categories (e.g., `through-routing`, `polygon-shared-ap`) the generator MUST pre-route via `shortest_polylink_to_polylinks` to pick endpoints that guarantee the matcher will produce the labelled behavior:
  - link-only (start and end on a link; no polygon candidate at any layer)
  - polygon-traversal (link → polygon → link, polygon visible as a candidate)
  - polygon-shared-AP crossing (passes through two polygons sharing an AP)
  - mid-polygon start (trajectory starts inside a polygon)
  - mid-polygon end (trajectory ends inside a polygon)
  - fully-inside (every GPS point inside one polygon)
  - through-routing (polygon is the cheapest routing shortcut between two links but no GPS observation falls inside it)
  - off-network noise (some GPS points placed >> `radius` from any edge)
  - short trips (2-3 points)
  - duplicate-points (consecutive identical GPS coordinates)
- **FR-005**: Every trace ID MUST be unique within the batch.
- **FR-006**: Every trace's coordinates MUST lie within the network's bounding box (with a small documented margin for off-network-noise traces; ≤ 200 m).

#### Matcher invocation

- **FR-007**: The validation harness MUST invoke polymatch against the generated trace set using the real-area network + polygon + access-point layers as input. The CRS units are meters; the harness MUST use the following POLYMATCHConfig values (meters where applicable): `radius = 300`, `gps_error = 50`, `boundary_epsilon = 1e-6`, `through_penalty_factor = 1.5`, `k = 8`. These match POLYMATCHConfig's defaults and are appropriate for a ~5 km × 8 km network with edge lengths in the 50-500 m range. The harness MUST set them explicitly (not rely on the implicit defaults) so a future change to those defaults does not silently move the test goalposts.
- **FR-008**: The harness MUST also invoke `weightmatch` against the same network and trace set in link-only mode for the SC-002 regression check (FR-013 below).
- **FR-009**: The full real-network validation — fixture load + per-trace polymatch matching + per-trace weightmatch matching (for `link-only` traces) + per-invariant verification + summary emission — MUST complete in under 60 seconds for a 200-trace batch on a single core (no OpenMP). This matches the SC-001 budget.
- **FR-010**: The matcher MUST NOT crash or throw an exception for any trace in the batch, including all edge-case categories.

#### Output verification

- **FR-011**: For every matched result, the harness MUST verify the hybrid `cpath` topology: every consecutive link-link pair shares a network node; every link↔polygon transition uses one of that polygon's APs. This invariant applies to *all* traces regardless of `category`.
- **FR-012**: For every `PolygonSegment` with `is_through == true`, the harness MUST verify that both `entry_ap` and `egress_ap` are populated (not `kNoAccessPoint`).
- **FR-013**: For traces whose `category` column is `link-only`, the harness MUST verify that polymatch's `opath`/`cpath` are byte-identical to `weightmatch`'s `opath`/`cpath` on the same trace (SC-002 carried forward). Traces with other category labels are exempt from this invariant.
- **FR-014**: For every `PolygonSegment`, the harness MUST verify that `distance_inside` is finite and non-negative.
- **FR-015**: The harness MUST aggregate violations across all traces (no fail-fast). For each invariant (FR-011 through FR-014) it MUST emit a summary line with the pass count, the fail count, and the first ≤ 10 failing trace IDs. After every trace has been processed, the harness MUST fail the test once if any invariant has any failures; otherwise it MUST pass.
- **FR-016**: The validation suite MUST run as part of `polymatch_test` under the dedicated Catch2 tag `[real_network]`. The tag MUST be a default-on tag (no leading dot), so a plain `./polymatch_test` invocation runs both the existing synthetic-fixture suite and this real-network suite together — the SC-001 60-second budget is set so this single command remains a fast regression gate. The tag exists to give developers a filter to run only the real-network suite during iteration (`./polymatch_test '[real_network]'`), not to skip it by default.

### Key Entities

- **Real-area fixture set**: Three shapefiles in `test/data/polymatch/real_example_area/`:
  - `network.shp` — `wkbLineString` road edges in a projected CRS (GDA 1994 Australia Albers, units: meters). Bounding box approximately 5 km × 8 km.
  - `polygons.shp` — `wkbPolygon` features representing the polygon layer.
  - `access_points.shp` — `wkbPoint` features representing access points.
- **Trace category**: A label (one of the ten in FR-004) that describes the topological / geometric class of a synthetic trip. Stored in the trace CSV's `category` column. Drives both generation distribution (the generator stamps the label) and post-match verification (the harness picks per-trace invariants by category).
- **Trace batch CSV**: The deterministic output of the generator — a single `id;geom` CSV consumed by polymatch and by the verification harness. Committed to the repo as `test/data/polymatch/real_example_area/trips.csv`. Regenerated by hand when the underlying shapefiles change; CI reads the committed file.
- **Match-invariant report**: A structured summary produced by the verification harness, listing each invariant from FR-011 through FR-014 with pass count, fail count, and a short list of failing trace IDs.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A single command (`polymatch_test '[real_network]'` or equivalent) runs the full real-network validation in under 60 seconds on a single core and prints a structured summary of the result.
- **SC-002**: The trace-generation step is fully deterministic: running the generator twice with the same seed produces byte-identical output.
- **SC-003**: At least 200 traces are generated with all ten documented categories represented at least 20 times each.
- **SC-004**: 100% of traces in the batch run through polymatch without crash or exception.
- **SC-005**: 100% of matched outputs satisfy the cpath-topology invariant (FR-011).
- **SC-006**: 100% of `is_through` polygon segments have both entry and egress APs populated (FR-012).
- **SC-007**: For every link-only trace, polymatch's matched path is bit-identical to weightmatch's matched path on the same input (FR-013 / SC-002 of polymatch spec).
- **SC-008**: A reviewer reading the test output can identify which invariant failed (and for which trace IDs) without re-running the suite.

## Assumptions

- The real-example fixtures in `test/data/polymatch/real_example_area/` are well-formed: polygons are valid, access points lie on polygon boundaries within the default `boundary_epsilon = 1e-6` (measured in meters, so this is effectively zero tolerance — the fixtures were prepared to exact precision), and access-point `node_id`s match the network's node ID scheme where link attachment is intended.
- The fixture data is large enough to host all ten trace categories without artificial extrapolation. If a category cannot be produced from the actual network topology (e.g., no two polygons share an AP), the generator MUST emit fewer or zero traces in that category and the harness MUST mark the corresponding invariant test as "skipped — no traces produced" rather than failing.
- The fixture CRS (GDA 1994 Australia Albers, meters) is preserved end-to-end — no on-the-fly reprojection is performed by any tool in this pipeline (consistent with the convention documented in `specs/000-pre-existing/contracts/network-shapefile.md`).
- Polymatch's existing feature spec (`specs/001-polymatch-algorithm/spec.md`) is authoritative for matcher behavior. This feature only adds a validation harness on top; it does NOT change the matcher.
- The trace generator is implemented in C++ (using GDAL/OGR + Boost.Geometry already available to the test target) and lives alongside the existing `polymatch_test.cpp` fixture setup, not as an external Python script — so the entire pipeline runs in CI without extra dependencies. (The Python `test/generate_polymatch_test_data.py` script remains the offline-reference equivalent for the small synthetic fixture; this feature does not extend or replace it.)
- "Real-network validation" is scoped to the provided real-example area only. Testing against other real road networks (e.g., a different city's data) is out of scope for this feature.
