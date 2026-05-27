# Quickstart: Real-Network Validation

A developer's-eye view of regenerating the trace batch and running the validation suite. Refer to [spec.md](./spec.md) for the feature scope and [research.md](./research.md) for design choices.

## Prerequisites

Same as the rest of the project: C++11 toolchain, CMake ≥ 3.5, GDAL ≥ 2.2, Boost ≥ 1.56, OpenMP. The real-area input shapefiles are committed at `test/data/polymatch/real_example_area/`.

## Build

```bash
cd /workspace
mkdir -p build && cd build
cmake ..
make -j$(nproc) polymatch polymatch_traces_gen polymatch_test
```

Three new things compared to feature 001:

- `polymatch_traces_gen` — new executable for offline trace generation.
- `polymatch_test` — gains a `[real_network]` Catch2 tag (already-built target).
- The CMake configuration injects `FMM_REAL_EXAMPLE_DIR` into the test binary (mirrors the existing `POLYMATCH_FIXTURE_DIR` / `FMM_TEST_DATA_DIR` pattern).

## Regenerate the committed trace CSV

You only need to do this when the underlying shapefiles change (different polygon set, different access points). Otherwise the committed CSV is authoritative.

```bash
./polymatch_traces_gen \
  --network        ../test/data/polymatch/real_example_area/network.shp \
  --polygons       ../test/data/polymatch/real_example_area/polygons.shp \
  --access_points  ../test/data/polymatch/real_example_area/access_points.shp \
  --seed           2026 \
  --output         ../test/data/polymatch/real_example_area/trips.csv

# Verify determinism (optional, but recommended after generator changes):
./polymatch_traces_gen \
  --network ../test/data/polymatch/real_example_area/network.shp \
  --polygons ../test/data/polymatch/real_example_area/polygons.shp \
  --access_points ../test/data/polymatch/real_example_area/access_points.shp \
  --seed 2026 \
  --output /tmp/trips_verify.csv
diff ../test/data/polymatch/real_example_area/trips.csv /tmp/trips_verify.csv
# Expected: no output. Any diff = a determinism bug.
```

Commit the updated CSV alongside whichever shapefile change motivated it.

## Run the validation suite

The default `./polymatch_test` invocation runs all 44 existing test cases plus the new real-network suite — no extra flag needed.

```bash
cd build
./polymatch_test
```

To run only the real-network suite (e.g., for fast iteration during development):

```bash
./polymatch_test '[real_network]'
```

Expected output on a clean pass:

```text
===============================================================================
All tests passed (NNN assertions in NN test cases)
```

## Reading the output on failure

When an invariant fails, the harness prints a structured summary like this just before the Catch2 final-status line:

```text
[real_network] Validation summary across 240 traces:
  cpath-topology           : 240 pass /   0 fail
  is-through-has-aps       : 240 pass /   0 fail
  link-only-eq-weightmatch :  20 pass /   2 fail  (first failing trace IDs: 1003, 1017)
  distance-inside-finite   : 240 pass /   0 fail
```

The first 10 failing trace IDs per invariant are listed inline. Re-run the suite with verbose output for one trace to see the specifics:

```bash
./polymatch_test '[real_network]' -s -d yes
```

To debug a specific failing trace, grep the committed CSV:

```bash
grep "^1003;" ../test/data/polymatch/real_example_area/trips.csv
```

…then feed that single trace into a one-off polymatch invocation for inspection.

## Acceptance walkthrough

| Spec requirement | Where to verify |
|---|---|
| SC-001 — full suite < 60 s single-core | `./polymatch_test '[real_network]'` from `build/`; check Catch2's reported wall-time. |
| SC-002 — deterministic generator | Two runs of `polymatch_traces_gen` with the same seed → byte-identical CSV (verified via `diff` above). |
| SC-003 — ≥ 200 traces, ≥ 20 per (supported) category | Harness's "Validation summary" lists the trace count and per-category populations. Categories the network can't host are explicitly tagged "skipped — no traces produced." |
| SC-004 — 100% completion, no crashes | `./polymatch_test` exits 0; no `[critical]` log lines. |
| SC-005 — cpath-topology 100% pass | "cpath-topology" line shows fail count 0. |
| SC-006 — `is_through` APs 100% present | "is-through-has-aps" line shows fail count 0. |
| SC-007 — link-only ≡ weightmatch | "link-only-eq-weightmatch" line shows fail count 0. |
| SC-008 — reviewer can identify failures by trace ID | Failing trace IDs are inlined in the per-invariant summary. |

## Common issues

| Symptom | Likely cause | Fix |
|---|---|---|
| `polymatch_traces_gen` silently produces zero traces in one category | Real fixture's topology can't host that category (e.g., no two polygons share an AP). | This is by design (research.md R7); the harness skips the category. Verify with `cut -d';' -f3 trips.csv \| sort \| uniq -c`. |
| Determinism diff after generator code change | RNG call sequence changed; expected. | Re-commit the updated CSV. |
| Determinism diff with no generator changes | Compiler optimization is reordering randomized work in parallel sections, or `std::mt19937_64`'s SO-level implementation drifted. | Verify the generator is single-threaded; check libstdc++ version parity across machines. |
| Real-network suite reports `link-only-eq-weightmatch` failures | Polymatch's link-only fallback regressed against weightmatch. | Inspect the failing trace IDs; check that link-only mode still delegates to `WEIGHTMATCH::match_traj` byte-for-byte. |
