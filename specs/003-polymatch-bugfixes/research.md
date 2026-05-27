# Phase 0 Research: Polymatch Bug Fixes from Real-Network Validation

This phase pins down the root cause of the two failing real-network invariants by reading the matcher code, walking the data flow, and forming a hypothesis the implementation phase can verify in code.

---

## R1. `is_through` semantics bug (US0)

**Decision**: The fix lives in `POLYMATCH::build_hybrid_path` (`src/mm/polymatch/polymatch_algorithm.cpp`). Today the loop that processes opath candidates calls `record_inside(matched_point, layer_idx)` only when `cand.inside == true`. The fix calls it whenever a polygon Viterbi candidate is the per-layer match, regardless of `inside`. The PolyCandidate's `inside` flag is repurposed as advisory only (it still informs the matched_point's geometric origin and the distance_inside formula); it stops being the sole gate on `has_inside_obs`.

**Rationale**:

- `PolyCandidate::inside` is set in `build_candidates` based on which of two queries returned the polygon: `polygons_containing(gps)` (strict `covered_by`, sets `inside=true`) vs `polygons_within_radius(gps, radius)` (within-radius but outside boundary, sets `inside=false`).
- When the matcher's Viterbi picks a polygon candidate at any GPS layer, the trip is *demonstrably* observed at that polygon — the spec's semantic for `is_through == true` ("no GPS observation associated with the polygon") is violated by the current logic if the picked candidate's `inside` flag was false.
- The fix preserves the existing inside-only semantics ONLY for the `record_inside` book-keeping that drives `distance_inside` (R12 of spec 001's research). The two distinct concerns — "was the polygon a Viterbi candidate?" (drives is_through) vs "did the GPS point fall strictly inside?" (drives distance_inside formula choice) — were conflated in feature 001.

**Concrete code change shape**:

```cpp
// Today:
if (b->inside) {
  record_inside(b->matched_point, (int)(i + 1));
}

// After:
if (b->is_polygon()) {
  // Every polygon Viterbi candidate match counts as an inside observation for
  // is_through purposes — see specs/003 clarification Q3.
  record_inside(b->matched_point, (int)(i + 1));
}
// (inside flag remains exposed on PolyCandidate for the distance_inside math)
```

This touches the three call sites of `record_inside` inside `build_hybrid_path`: the initial-element handling, the link→polygon transition, the polygon→polygon (same polygon) transition. Each currently gates on `inside`; each becomes unconditional once we're inside a polygon-candidate branch.

**Alternatives considered**:

- *Set the candidate's `inside=true` for both query types in `build_candidates`* — Rejected. The flag is also used in `distance_inside` calculation to know whether to include AP-to-first-inside segments. Changing its meaning ripples beyond the bug.
- *Add a second flag like `is_through_candidate`* — Rejected. The fix is simpler: a polygon Viterbi candidate is intrinsically an "observed-at-polygon" event for is_through purposes.

---

## R2. cpath-topology violations on traces 1313 + 1314 (US1)

**Decision**: Hypothesis — the bug is in `build_hybrid_path`'s **polygon → link transition** branch (a.is_polygon() && b.is_link()). The matcher iterates the polygon's APs to find a best (lowest-cost) egress route, but the AP it ultimately picks for `set_egress(best_ap)` is *not always* the AP whose attached edges appear at the head of `chosen_segs`. The result: the cpath has `[polygon_id, edge_X]` where `edge_X` is the first emitted edge of `chosen_segs`, but the `best_ap` we recorded as `egress_ap` corresponds to a DIFFERENT AP — leading to the harness's "endpoint not an AP of polygon" failure when the harness scans `cpath` looking for the AP-incident edge.

Actually re-reading the matcher code at `build_hybrid_path`:

```cpp
for (auto ap_idx : aps) {
  ...
  for (auto src_edge : ap.attached_edges) {
    ...
    if (cost < best) {
      best = cost;
      best_ap = ap.node_id;
      chosen_segs = paths[0].edges;  // starts at src_edge of THIS AP
    }
  }
}
...
set_egress(best_ap);   // sets the segment's egress_ap to whichever AP won
for (auto it = chosen_segs.begin(); it != chosen_segs.end(); ++it) {
  push_link(network_.get_edges()[*it].id);  // emits chosen_segs in order
}
```

So `best_ap` and `chosen_segs` are written together — they should be consistent. Then where's the bug?

A more careful hypothesis: the cpath-topology failure isn't necessarily at the polygon→link emission. It could be at a *different* pair — specifically, when the initial-element handling emits a polygon (mid-polygon-start) followed by a `polygon → polygon` (same polygon) transition that, in turn, records candidates differently from what subsequent transitions expect. Or: when iter i is `polygon → polygon (same)` followed by iter i+1 = `polygon → link`, the polygon→link branch uses `a->polygon_index` from the LATER opath candidate (b would be a link). But the cpath at this point already has the polygon pushed; the polygon→link emit appends edges. The first appended edge SHOULD be `chosen_segs[0]` which is an attached_edge of the winning AP.

The actual root cause likely needs a printf-style debug walk through trace 1313 or 1314. **Phase 0 deliverable**: instrument `build_hybrid_path` with debug logging (under `SPDLOG_DEBUG`), match trace 1313, dump the (a, b) transition each iteration emits + the cpath state. Compare to the harness's failure message ("edge 1661 endpoints ... are not APs of polygon 28"). Identify which `build_hybrid_path` branch produced the bad pair.

**Rationale**: We don't yet have enough information to predict the fix without instrumenting; we DO have enough information to write a deterministic reproduction (just run trace 1313 from the committed CSV).

**Hypotheses to verify in implementation**:

1. The polygon→link branch's `best_ap` / `chosen_segs` coupling is fine; the bug is elsewhere. Possible: the LINK→LINK transition emits a polygon sub-vertex via the `shortest_polylink_to_polylinks` Dijkstra, and the sub-vertex's recorded `entry_ap` doesn't match the actual sub-vertex's AP (a mismatched lookup in `polygon_of(v)` vs `ap_of(v)`).
2. The polygon→link branch picks an AP whose `attached_edges[0]` happens to lead to a `chosen_segs[0]` that, due to LinkGraph directionality, is an out-edge of a node that ISN'T the AP. (Recall: AccessPoint.attached_edges includes both incoming and outgoing edges of the AP node. Dijkstra `shortest_edge_to_edges` from an INCOMING edge has unintuitive starting semantics — the "start" edge's path doesn't physically begin at the AP node.)
3. There's a mismatch between which polygon's APs the harness checks vs the polygon the matcher used. The harness extracts `pid = -cpath[i+1]` and checks polygon `pid`'s APs — but if the matcher accidentally emits the polygon id of `a->polygon_index` (not `chosen_segs`'s actual polygon), the harness would look at the wrong AP set.

**Alternatives considered**:

- *Pre-fix without diagnosis* — Rejected. The 2 failures are deterministic; a 20-minute instrumentation walk gives us the actual cause and we fix once.
- *Skip the cpath-topology fix and live with the tolerance* — Rejected; explicit user direction (Clarification Q3 phrasing) treats the violations as real bugs, not acceptable matcher quirks.

---

## R3. Perf diagnostic (US2)

**Decision**: Take a single wall-time measurement of the [real_network] suite, break it down by adding scoped timers to `polymatch_test`'s real-network TEST_CASE around (a) fixture load (one-shot RealAreaFixture construction), (b) trace iteration matching loop, (c) per-trace polymatch::match_traj only, (d) per-trace weightmatch::match_traj only (for link-only traces), (e) the ViolationLedger / invariant verification per trace, (f) trips.csv load. Print the breakdown via `WARN()` so it's visible in CI logs without failing the test.

**Rationale**:

- Instrumentation is reversible and isolated: a few `std::chrono::steady_clock` accumulators around each phase, then a one-time `WARN()` emit at end-of-test.
- The diagnostic is one-shot — we're not building a performance harness, we're answering "where does the time go?" with a single measurement.
- Per Clarification Q1, optimization is out of scope unless an obvious quick win surfaces. The output is a 1-2 paragraph summary committed to 002's tasks.md `Deferred Follow-Ups` section.

**Expected breakdown shape** (predicted; will be replaced with measurements):

| Phase | Time | Note |
|---|---|---|
| RealAreaFixture construction (one-shot) | 1-3 s | LinkGraph + PolyLinkGraph build + R-tree |
| trips.csv load + WKT parsing | < 1 s | 200 rows |
| Per-trace match_traj (polymatch) | ~50-65 s | 200 × ~250-350 ms estimated |
| Per-trace match_traj (weightmatch, link-only only) | < 5 s | 20 × ~250 ms |
| Per-trace invariant checks | < 1 s | mostly index walks |

If the prediction is right, optimization options would target the per-trace polymatch cost, which is matcher-internal work that this feature explicitly defers (see Clarification Q1).

**Alternatives considered**:

- *External profiler (perf record / callgrind)* — Rejected. Adds toolchain dependency; the scoped-timer approach is reproducible without external tools.
- *Bench tag like `[.real_network_perf]` separately* — Rejected. The diagnostic only needs to run once for the deliverable; making it a permanent fixture is over-engineering.

---

## R4. Harness assertion changes

**Decision**:

- Drop the first/last boundary exemption inside `check_is_through_has_aps`. After the matcher fix (R1) takes effect, `is_through == true` should imply both APs set without exception.
- Tighten the `cpath_topology` end-of-test CHECK from `cpath_fail <= 5` to `cpath_fail == 0`.

Both changes go into the same `polymatch_test.cpp` edit. Order of operations:

1. Tighten the harness's check_is_through_has_aps first (drop the exemption). Run the suite — `is-through-has-aps` should now fail on the same mid-polygon-start traces that the matcher mishandles.
2. Apply the matcher fix (R1). Re-run — `is-through-has-aps` returns to 200/0 pass.
3. Investigate trace 1313/1314 per R2. Apply the cpath-topology fix.
4. Tighten the `cpath_fail` tolerance to `== 0`. Re-run — all four invariants pass at 200/0.

**Rationale**: The "tighten harness first, fix matcher second" sequence is the spec-mandated TDD discipline (Constitution III). Each step is reversible.

---

## R5. Determinism + non-regression checks

**Decision**: After each matcher edit, run the full default `./polymatch_test` (49 cases / 2860+ assertions in the synthetic suite plus the 10 real-network cases). Confirm:

- `link-only-eq-weightmatch` stays at 20/0 pass on the real fixture (FR-009, strict per Clarification Q2).
- All synthetic-suite invariant assertions still pass byte-unchanged (FR-003).
- Any synthetic-suite specific-value assertion that changes due to legitimate matcher behavior is updated with an inline comment naming the bug-fix commit (FR-003, Clarification Q2).
- `polymatch_traces_gen --seed 2026` regenerates a CSV byte-identical to the committed `trips.csv` (FR-010 — generator behavior unchanged, since the matcher fix doesn't touch the generator).

**Rationale**: Standard non-regression sweep. The order matters — fix matcher, then run the *full* synthetic suite before touching any synthetic assertion, so we know exactly which ones moved.

---

## Summary

| # | Topic | Decision |
|---|---|---|
| R1 | `is_through` fix | `record_inside` called on every polygon Viterbi candidate, not just `inside=true` ones |
| R2 | cpath-topology fix | Diagnose via debug instrumentation on traces 1313/1314 — three hypotheses to verify |
| R3 | Perf diagnostic | Scoped `std::chrono` timers in the harness's TEST_CASE; `WARN()` emit; record summary in 002's Deferred Follow-Ups |
| R4 | Harness assertion changes | Drop boundary exemption + tighten cpath_topology tolerance; apply in TDD-first order |
| R5 | Non-regression | Full default `polymatch_test` after each matcher edit; verify byte-determinism of `trips.csv` regeneration |

All Phase 0 unknowns identified; the cpath-topology root-cause discovery is itself part of Phase 1 work (the implementation phase will instrument and trace before changing logic).
