# Quickstart: Polymatch Bug Fixes from Real-Network Validation

How to reproduce, verify, and guard against regressing the two bugs fixed in this feature.

## Reproducing the failures (pre-fix baseline)

The two bugs reproduce deterministically on `master`/pre-fix code against the committed real-network fixtures.

```bash
cd /workspace
mkdir -p build && cd build && cmake .. >/dev/null
make polymatch_test -j$(nproc)
./polymatch_test '[real_network]' 2>&1 | tail -20
```

Pre-fix ledger output (recorded in feature 002):

```text
[real_network] Validation summary across 200 traces:
  cpath-topology            :  198 pass /    2 fail (first failing trace IDs: 1313 1314)
      sample reasons:
        - trace 1313: edge 1661 endpoints (src=454311032 tgt=12655159115) are not APs of polygon 28
        - trace 1314: edge 1484 endpoints (src=267204714 tgt=267204715) are not APs of polygon 172
  is-through-has-aps        :  200 pass /    0 fail   (← masked by harness boundary exemption)
  link-only-eq-weightmatch  :   20 pass /    0 fail
  distance-inside-finite    :  200 pass /    0 fail
```

The cpath-topology line shows 2 known failures (within the tolerated `<= 5` budget). The is-through-has-aps line shows 0 failures only because of the boundary exemption the harness applies; under the strict invariant this feature restores, it would show some mid-polygon-start / mid-polygon-end failures too.

## Diagnosing the cpath-topology bug

The implementation phase will instrument `polymatch::build_hybrid_path` with `SPDLOG_DEBUG` traces around each `push_link` / `push_polygon` call, then re-run the suite with that one trace isolated:

```bash
# Hand-extract trace 1313's geometry to a one-row CSV:
head -1 ../test/data/polymatch/real_example_area/trips.csv > /tmp/one.csv
grep '^1313;' ../test/data/polymatch/real_example_area/trips.csv >> /tmp/one.csv

# Then invoke polymatch directly with verbose logging:
./polymatch \
  --network ../test/data/polymatch/real_example_area/network.shp \
  --polygons ../test/data/polymatch/real_example_area/polygons.shp \
  --access_points ../test/data/polymatch/real_example_area/access_points.shp \
  --gps /tmp/one.csv --output /tmp/one_out.csv \
  --weight cost --log_level 1 \
  -k 8 -r 300 -e 50 --through_penalty_factor 1.5
```

Inspect the printed cpath against polygon 28's AP node IDs in `access_points.shp` to identify which transition produces the bad edge.

## Verifying the fix (post-fix gate)

After the matcher fix lands, the same `[real_network]` invocation MUST report all four invariants at 0 fail:

```text
[real_network] Validation summary across 200 traces:
  cpath-topology            : 200 pass /   0 fail
  is-through-has-aps        : 200 pass /   0 fail
  link-only-eq-weightmatch  :  20 pass /   0 fail
  distance-inside-finite    : 200 pass /   0 fail
```

The harness's `check_is_through_has_aps` will have the boundary-exemption block removed, and the trailing `CHECK(cpath_fail <= 5)` will be tightened to `CHECK(cpath_fail == 0)`.

## Non-regression sweep

After any candidate fix:

```bash
# Default synthetic suite still green (feature 001 + 002):
cd build && ./polymatch_test
# Expected: All tests passed (>= 2860 assertions in >= 54 test cases)

# All six binaries:
cd test
../algorithm_test && ../fmm_test && ../network_test && \
  ../network_graph_test && ../weightmatch_test && ../polymatch_test
# Expected: all six print "All tests passed"

# Generator is unchanged → committed trips.csv MUST still be byte-identical:
cd ../..
./build/polymatch_traces_gen \
  --network test/data/polymatch/real_example_area/network.shp \
  --polygons test/data/polymatch/real_example_area/polygons.shp \
  --access_points test/data/polymatch/real_example_area/access_points.shp \
  --seed 2026 --output /tmp/trips_verify.csv >/dev/null 2>&1
diff test/data/polymatch/real_example_area/trips.csv /tmp/trips_verify.csv && \
  echo "TRIPS.CSV DETERMINISM: PASS"
```

## Performance diagnostic walkthrough (US2)

The implementation will add scoped `std::chrono::steady_clock` timers to the harness's main `TEST_CASE("Real-network validation against committed trace batch ...")` block, around five phases:

1. RealAreaFixture construction (one-shot at first `real_fixture()` call).
2. `trips.csv` load + WKT parsing.
3. Per-trace `POLYMATCH::match_traj` accumulated.
4. Per-trace `WEIGHTMATCH::match_traj` accumulated (link-only traces only).
5. Per-trace invariant checks accumulated.

Output goes to `WARN()` so it appears in CI logs without failing the test. After the suite runs once, the resulting numbers are recorded in `specs/002-real-network-validation/tasks.md` under "Deferred Follow-Ups → Performance follow-up" with at least one number per category — that is the deliverable.

## What this feature does NOT change

- `trips.csv` — bit-identical pre- and post-fix.
- The trace generator (`polymatch_traces_gen`) — no source changes.
- The real-area shapefiles (network / polygons / access_points).
- Polymatch's CLI flags or output schema (the polymatch-output.md contract is honored as before).
- Spec 001 or 002's SC budgets (002's SC-001 wall-time remains deferred; this feature doesn't move it).

## Common issues during the fix

| Symptom | Likely cause | Action |
|---|---|---|
| `is-through-has-aps` shows new failures after removing the boundary exemption | The matcher fix (R1) hasn't landed yet, or `record_inside` isn't being called on all polygon Viterbi candidate branches. | Re-check that ALL three call sites — initial element, link→polygon, polygon→polygon (same) — call `record_inside` unconditionally for polygon candidates. |
| `link-only-eq-weightmatch` regresses | The matcher fix changed link→link transition behavior. | Strict gate per Clarification Q2 — discard this candidate fix; try a different one that touches only polygon-related code paths. |
| `trips.csv` regeneration diff | The generator was inadvertently touched. | The fix must not edit `src/network/trace_generator.{hpp,cpp}` or `src/app/polymatch_traces_gen.cpp`. Revert any generator changes. |
| Synthetic-suite specific-value assertion fails | Legitimate behavior change from the fix (per Clarification Q2). | Verify the new value still satisfies every invariant; update the assertion with an inline comment naming the bug-fix commit. |
