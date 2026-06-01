# Feature Specification: Polymatch Bug Fixes from Real-Network Validation

**Feature Branch**: `003-polymatch-bugfixes`

**Created**: 2026-05-27

**Status**: Draft

**Input**: User description: "fix the bugs identified in above testing expansion"

## Clarifications

### Session 2026-05-27

- Q: Is US2 (performance investigation + budget decision) in scope, or do we keep this feature laser-focused on US1? → A: Profile-only. Produce a paragraph-length profile summary identifying the dominant cost contributor for the real-network suite; defer any actual optimization or SC-001 budget revision to a separate follow-up *unless* the profile surfaces a single obvious quick win.
- Q: If the bug fix changes synthetic-test outputs, what's the policy? → A: Hold link-only-eq-weightmatch strictly; accept changes elsewhere if the new value is correct. `link-only-eq-weightmatch` (FR-009) remains a strict gate — the link-only path must continue to produce output bit-identical to weightmatch. Other synthetic-test assertions checking specific edge IDs *may* be updated to match the fix's corrected output, provided the new output still satisfies every other invariant (cpath topology, is_through-has-aps, distance_inside finite). Updating a specific-value assertion requires a comment naming the bug-fix commit that motivated the change.
- Q: Fate of the `is-through-has-aps` boundary exemption introduced in 002? → A: **Tighten the matcher.** `is_through == true` semantically means "polygon was a routing detour with no GPS observation associated with it." A trip that starts or ends inside a polygon is *demonstrably* associated with that polygon, so its first/last polygon segment MUST have `is_through == false` (even when the GPS point sat within radius rather than strictly inside the geometry). The matcher's current logic sets `has_inside_obs` only on `inside=true` candidates; the fix sets it whenever a polygon Viterbi candidate (from either `polygons_containing` or `polygons_within_radius`) is matched at any layer. Once the matcher is correct, the 002 harness's `check_is_through_has_aps` boundary exemption MUST be removed — the strict "is_through implies both APs set" invariant holds for ALL polygon segments without exception.

## User Scenarios & Testing *(mandatory)*

### User Story 0 — Correct the `is_through` Semantics (Priority: P1)

The same developer notices that mid-polygon-start traces in the real-network suite emit polygon segments with `is_through == true` even though the trip *began inside* the polygon. By the spec's definition, `is_through` means "the polygon was a routing detour with no GPS observation associated with it" — a trip that starts inside the polygon is, by definition, associated with it. The matcher's current logic ties `has_inside_obs` to the `inside` flag of polygon candidates, but `inside` only captures the strict `polygons_containing` case (geometric `covered_by`). Candidates coming from `polygons_within_radius` (close to the polygon but not strictly inside) get `inside=false` even when they're the Viterbi-selected per-layer match — leaving `has_inside_obs=false` → `is_through=true` for trips the matcher demonstrably saw at the polygon.

**Why this priority**: This bug is the *root cause* of the 002 harness's "is-through-has-aps boundary exemption" — the exemption was a defensive workaround for matcher behavior that turned out to be semantically wrong. Fixing the matcher lets the harness drop the exemption and enforce the strict invariant uniformly.

**Independent Test**: After this fix, every polygon segment in a `mid-polygon-start`, `mid-polygon-end`, or `fully-inside` trace MUST have `is_through == false`. The 002 harness's `check_is_through_has_aps` is updated to require both APs set on EVERY `is_through=true` polygon segment without first/last exemption, and `./polymatch_test '[real_network]'` still reports `is-through-has-aps : 200 pass / 0 fail`.

**Acceptance Scenarios**:

1. **Given** a trip with a polygon Viterbi candidate matched at any GPS layer, **When** the matcher assembles the polygon segment for that polygon, **Then** the segment's `is_through` is `false`, regardless of whether the candidate's `inside` flag was true or false.
2. **Given** a trip with no polygon Viterbi candidate at any layer but whose link↔link routing crosses a polygon as a Dijkstra shortcut, **When** the matcher emits the polygon segment, **Then** `is_through` is `true` AND both `entry_ap` and `egress_ap` are populated.
3. **Given** the 002 harness's `check_is_through_has_aps`, **When** the matcher fix is applied, **Then** the boundary exemption is removed and the strict invariant ("`is_through=true` ⇒ both APs set") holds across all 200 traces with zero violations.

---

### User Story 1 — Eliminate cpath-Topology Violations on Mid-Polygon-Start Trips (Priority: P1)

A developer reviewing the real-network validation suite sees two trace IDs (1313, 1314) consistently report `cpath-topology` violations. The harness's tolerance lets them pass (≤ 5 allowed), but the visible "2 fail" count is a constant red herring for anyone scanning the summary. The bug only surfaces on the real fixture — the synthetic suite doesn't reach it — so it was invisible until feature 002. They need the matcher to emit a topologically valid `cpath` for every trace category, including mid-polygon-start, so the real-network ledger shows `0 fail` for `cpath-topology` and the tolerance can be tightened back to 0.

**Why this priority**: The bug undermines the spec-001 FR-015 contract ("hybrid C_Path topology preserved through access points") on real-world data. Every other invariant on the real fixture is 100%; this is the lone reproducible matcher correctness failure surfaced by the 002 suite. Until it's fixed, the harness output carries a permanent "known-issue" footnote that obscures genuine future regressions.

**Independent Test**: `./polymatch_test '[real_network]'` reports `cpath-topology : 200 pass / 0 fail` (no "first failing trace IDs" line), the in-test cpath-topology tolerance is removed (`CHECK(cpath_fail == 0)` instead of `<= 5`), and trace IDs 1313 and 1314 produce matched outputs whose every consecutive `(link, polygon)` or `(polygon, link)` pair in `cpath` has the edge's endpoint matching one of that polygon's access-point node IDs.

**Acceptance Scenarios**:

1. **Given** the committed `trips.csv` from feature 002, **When** the validation suite runs, **Then** every emitted polygon ↔ edge transition in every trace's `cpath` has at least one edge endpoint that's a registered AP of the polygon (i.e., `cpath-topology` reports zero failures).
2. **Given** a developer modifies `build_hybrid_path`, **When** they re-run the real-network suite, **Then** the suite's `cpath-topology` summary still shows zero failures, OR if a new regression appears, the harness reports the regression by trace ID with the same precision the current ledger provides.
3. **Given** the 002 suite's `is-through-has-aps` boundary-exemption is currently used to mask first/last polygon segments with one absent AP, **When** the underlying matcher is reviewed, **Then** any cases where the matcher could legitimately set both APs at the boundary (and the exemption hides correctness) are either fixed or explicitly documented as accepted spec relaxations.

---

### User Story 2 — Diagnose the [real_network] Suite's Wall-Time Overshoot (Priority: P2)

The same developer's CI runs polymatch_test and notices the real-network slice takes about 68 seconds — slightly over the 60-second SC-001 target the 002 spec committed to. They don't need to fix the budget overshoot in this feature — that's deferred to a follow-up — but they do need to *know* where the time is going, so the follow-up has a starting point and so reviewers reading the spec aren't surprised that a "fast regression gate" is sluggish.

**Why this priority**: This is a diagnostic deliverable, not a fix. SC-001's "fast regression gate" intent is undermined by the overshoot, but ~13% over isn't urgent — what *is* urgent is having a documented cost breakdown so future optimization (or budget revision) is data-driven rather than guessed.

**Independent Test**: After this feature ships, a reviewer reading the 003 spec or the 002 deferred-follow-ups section can name the dominant cost contributor in the real-network suite (e.g., "polymatch::match_traj on the 15k-edge network averages ~X ms/trace; harness loop overhead is negligible") backed by a recorded measurement.

**Acceptance Scenarios**:

1. **Given** the real-network suite has been instrumented or profiled, **When** the developer reads the deliverable, **Then** the dominant cost contributor (matcher call vs harness loop vs invariant verification vs fixture load) is identified with at least one number per category.
2. **Given** the profile is captured, **When** the developer reviews it, **Then** they can decide whether an obvious "single quick win" optimization is available now, or whether the suite's cost is intrinsic and warrants a follow-up to either optimize or revise the 002 budget.

---

### Edge Cases

- **Off-network-noise traces**: After the bugfix, off-network-noise traces should still match cleanly without exception (matcher already handles them; this is a non-regression check, not a new requirement).
- **Mid-polygon-start with a polygon that has only one link-attached AP**: The matcher must still emit a topologically valid cpath when the egress AP is forced (no choice between APs).
- **Polygon-shortcut routes that the matcher picks during link↔link Dijkstra**: The bugfix MUST NOT regress the 001 behavior where polygon shortcuts surface as `is_through=true` polygon segments with both APs populated.
- **Determinism**: After the bugfix, the harness's `cpath-topology` output for a fixed seed and fixture set MUST be reproducible across re-builds (no NaN-driven non-determinism).

## Requirements *(mandatory)*

### Functional Requirements

#### Correctness fix — `is_through` semantics (User Story 0)

- **FR-001a**: `POLYMATCH::match_traj` MUST set `PolygonSegment.is_through = true` only when no polygon Viterbi candidate for that polygon was matched at any GPS layer along the trajectory. If the per-layer Viterbi picked any polygon candidate for polygon P at any GPS point in the trip (whether the candidate's `inside` flag was true or false), every emitted segment for P MUST have `is_through = false`.
- **FR-001b**: The 002 harness's `check_is_through_has_aps` MUST be tightened: remove the first/last boundary exemption introduced in feature 002. After the matcher fix, EVERY polygon segment with `is_through == true` MUST have both `entry_ap` and `egress_ap` populated (no `kNoAccessPoint` on either side).

#### Correctness fix — cpath topology (User Story 1)

- **FR-001**: For every matched trajectory output by `POLYMATCH::match_traj` against the real_example_area fixture, every consecutive `(link, polygon)` and `(polygon, link)` pair in `result.base.cpath` MUST satisfy: the link edge has at least one endpoint (source or target) whose node ID is the `node_id` of an access point that lists the polygon in its `polygons` field.
- **FR-002**: For every matched trajectory output by `POLYMATCH::match_traj` against the real_example_area fixture, every consecutive `(polygon_A, polygon_B)` pair in `result.base.cpath` (where A ≠ B) MUST be backed by at least one shared access point — an `AccessPoint` whose `polygons` field contains both polygon indices.
- **FR-003**: The synthetic-fixture suite (feature 001's tests inside `polymatch_test`) MUST continue to pass. *Invariant-style* assertions — `link-only-eq-weightmatch`, cpath-topology, is_through-has-aps (strict, no boundary exemption), distance_inside finite — MUST hold byte-unchanged. *Specific-value* assertions (e.g., "matched cpath equals `[1, 2, -7]`") MAY be updated to reflect the bug fix's corrected output, provided every invariant still passes and the changed assertion carries an inline comment naming the bug-fix commit that motivated the change. The total assertion count is allowed to grow or stay the same; it MUST NOT shrink without justification.
- **FR-004**: Once FR-001 / FR-002 are met, the real-network harness's tolerance for `cpath-topology` failures MUST be removed — the relevant `CHECK` becomes `cpath_fail == 0` instead of `cpath_fail <= 5`.

#### Performance diagnostic (User Story 2)

- **FR-005**: A profile of the real-network harness MUST be produced that identifies the top cost contributor (matcher call, harness loop, fixture load, or invariant verification) to the suite's wall time, with at least one timed measurement per category.
- **FR-006**: The profile summary (1-2 paragraphs) MUST be recorded in the 002 spec's Deferred Follow-Ups section so future optimization work has a documented starting point.
- **FR-007**: If — and only if — the profile reveals a single obvious quick win (e.g., a redundant per-trace allocation, a misordered loop, an O(n²) check that should be O(n)), the developer MAY apply that single fix and record the post-fix wall time alongside the pre-fix number. Larger optimization or any 002 SC-001 budget revision is OUT of scope for this feature and remains a deferred follow-up.

#### Non-regression guards

- **FR-008**: The fix MUST NOT regress polymatch's through-routing behavior for trips where the polygon truly was an unobserved routing shortcut. Specifically: trips with NO polygon Viterbi candidate at any layer whose link↔link Dijkstra crosses a polygon via `PolyLinkGraph` MUST still produce `is_through == true` polygon segments with both APs populated. The fix only narrows `is_through == true` away from cases where the polygon WAS a per-layer Viterbi pick (FR-001a).
- **FR-009**: The fix MUST NOT regress the link-only-equals-weightmatch invariant — this gate is *strictly held* per the clarifications session. After the fix, the 002 suite's `link-only-eq-weightmatch` count MUST remain at `20 pass / 0 fail` (or `N pass / 0 fail` for whatever N is the current count of `link-only` traces). If a candidate fix would break this invariant for even one trace, a different fix MUST be chosen.
- **FR-010**: The fix MUST NOT change the committed `trips.csv` — that file is the artifact of feature 002's generator and an unrelated correctness fix in the matcher should not alter generator behavior. After the fix, regenerating the CSV with the same seed must still produce the same bytes.

### Key Entities

- **`PolyMatchResult.base.cpath`**: The hybrid sequence of link IDs (positive) and polygon IDs (negative) that the matcher emits per trajectory. The bug under fix is that some entries adjacent to a polygon entry don't reference the polygon's APs.
- **`PolygonLayer::aps_for_polygon(p)`**: The authoritative mapping from a polygon to its access-point indices. The fix must ensure the matcher's emitted edges, where they touch a polygon, are consistent with this mapping.
- **`AccessPoint::attached_edges`**: The set of network edges incident to a link-attached AP's network node. The matcher's `build_hybrid_path` polygon↔link routines use this to choose the first/last edge of the polygon's adjacent segment.
- **Bug traces 1313, 1314**: Two `mid-polygon-start` entries in feature 002's committed `trips.csv`. They reproduce the failure deterministically and form the regression-test baseline for FR-001.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: After the fix, `./polymatch_test '[real_network]'` reports `cpath-topology : 200 pass / 0 fail` AND `is-through-has-aps : 200 pass / 0 fail` with the harness's boundary exemption removed (zero failures, no "first failing trace IDs" line in the ledger for either invariant).
- **SC-002**: The fix does not regress any synthetic-fixture test. After the fix, `./polymatch_test` (default invocation) still reports `All tests passed`. Assertion count MUST be ≥ 2860 (the pre-fix floor); specific-value assertions touched by the bug fix are allowed to be updated provided each touched assertion carries a comment naming the bug-fix commit and every invariant-style assertion still passes byte-unchanged.
- **SC-003**: The real-network harness's `cpath-topology` `CHECK` is tightened from `cpath_fail <= 5` to `cpath_fail == 0` (the tolerance that masked the discovered failures is removed). The `is-through-has-aps` invariant's first/last boundary exemption is also removed; the check requires both APs set on every `is_through == true` segment without exception.
- **SC-004**: All six existing test binaries (algorithm_test, fmm_test, network_test, network_graph_test, weightmatch_test, polymatch_test) continue to pass cleanly.
- **SC-005**: A 1-2 paragraph performance-profile summary is recorded in 002's Deferred Follow-Ups section, naming the dominant cost contributor with at least one measured number. Wall-time fix vs budget revision remains a deferred follow-up.
- **SC-006**: For traces 1313 and 1314 specifically, the matcher's emitted `cpath` is reviewable and topologically valid by inspection (a developer can read the edge IDs and AP node IDs and trace the path manually without finding an inconsistency).

## Assumptions

- The two failing trace IDs (1313, 1314) reproduce deterministically against the committed `trips.csv` and the committed real-area shapefiles. If a future regeneration of `trips.csv` changes these specific IDs, the spec's diagnostic detail in FR-001 / SC-006 still applies to whichever trace IDs exhibit the same failure pattern — the bug is in the matcher, not in the specific traces.
- The fix lives in `src/mm/polymatch/polymatch_algorithm.cpp` (the matcher) — most likely inside `POLYMATCH::build_hybrid_path` and possibly `POLYMATCH::transition_cost`. It does **not** modify the trace generator, the validation harness, or the polygon/access-point loaders.
- The relaxed `is-through-has-aps` boundary exemption introduced in feature 002 turned out to be masking a real matcher semantics bug (see Clarification Q3 + User Story 0). This feature removes the exemption from the 002 harness AND fixes the underlying matcher logic so `is_through` reflects "polygon was a routing detour, not a per-layer Viterbi pick."
- Feature 002's `polymatch_test [real_network]` suite is the regression gate for this fix. There is no separate test deliverable in 003 — we tighten 002's existing assertions (FR-004) and rely on its harness to catch regressions.
- Performance-side: a meaningful profile may be obtainable with `perf record` or similar tooling on the target machine; the deliverable is the *decision* (optimize vs adjust budget) and its rationale, not necessarily a specific code change.
