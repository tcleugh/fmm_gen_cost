---

description: "Task list for Polymatch Bug Fixes from Real-Network Validation"
---

# Tasks: Polymatch Bug Fixes from Real-Network Validation

**Input**: Design documents from `/workspace/specs/003-polymatch-bugfixes/`

**Prerequisites**: plan.md, spec.md, research.md, data-model.md, quickstart.md

**Tests**: TDD is mandatory per Constitution Principle III. The discipline this feature uses is "tighten harness first (assertions go red), then fix matcher (assertions return to green)." Every fix has a pre-existing 002 harness check it must pass — no new test files are created.

**Organization**: Three user stories from spec.md plus a small Polish phase. US0 + US1 are P1 correctness fixes; US2 is a P2 diagnostic.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Different files / no dependencies on incomplete tasks.
- **[Story]**: Maps to user stories from spec.md (US0, US1, US2).
- All file paths are absolute under `/workspace/`.

## Path Conventions

Surgical edits only. Source under `src/mm/polymatch/polymatch_algorithm.cpp`; harness under `test/polymatch_test.cpp`; 002 follow-ups under `specs/002-real-network-validation/tasks.md`. No new files.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Confirm the pre-fix baseline before changing anything, so we have a known starting point.

- [X] T001 Run `./polymatch_test '[real_network]'` on the current branch's HEAD and capture the full ledger summary (4 invariant lines + "first failing trace IDs" for cpath-topology) into `/tmp/003_baseline.txt`. This file is the diagnostic anchor for the rest of the work — do not commit it; it lives in `/tmp` for the duration of the iteration.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Tighten the harness's assertions so the bugs are visible as failures. Per Constitution III + research.md R4's TDD-first sequence: harness goes red BEFORE any matcher change.

**⚠️ CRITICAL**: After Phase 2, `./polymatch_test '[real_network]'` MUST be RED. That's the proof the assertions are doing their job. US0 + US1 turn them green.

- [X] T002 Tighten `check_is_through_has_aps` in `/workspace/test/polymatch_test.cpp`: remove the first/last boundary-exemption block. Every `is_through == true` segment, regardless of position in `polygon_segments`, MUST require both `entry_ap` and `egress_ap` set. Update the function's doc-comment to drop the "EXCEPT at trajectory boundaries" paragraph and explain that the strict invariant is restored by 003.
- [X] T003 Tighten the cpath_topology tolerance at the end of the main `[real_network]` `TEST_CASE` in `/workspace/test/polymatch_test.cpp` from `CHECK(cpath_fail <= 5)` to `CHECK(cpath_fail == 0)`. Drop the "Allow up to 5 ..." comment block above the check; replace it with a comment naming the 003 spec.
- [X] T004 Update unit test `check_is_through_has_aps invariant function (T026)` in `/workspace/test/polymatch_test.cpp` to reflect the strict invariant. The boundary-exemption test cases (single-segment trajectory, first/last with one AP absent) MUST now FAIL the invariant; assertions inverted accordingly. Keep a separate sub-case verifying the strict "both APs ⇒ pass" case is still accepted.
- [X] T005 Build + run `./polymatch_test '[real_network]'`. Verify it goes RED — at minimum the cpath_topology `CHECK == 0` fails, and `check_is_through_has_aps` reports a non-zero failure count on mid-polygon-start / mid-polygon-end traces that were previously masked by the exemption. Capture the post-tighten ledger output to `/tmp/003_red.txt`.

**Checkpoint**: harness is strict, real-network suite is red on cpath_topology + is-through-has-aps. US0 + US1 work begins.

---

## Phase 3: User Story 0 (Priority: P1) — Correct the `is_through` Semantics

**Goal**: `POLYMATCH::match_traj` sets `is_through=true` only when no polygon Viterbi candidate was matched at any GPS layer. Per research.md R1, change `record_inside` call sites in `build_hybrid_path` to fire on every polygon Viterbi candidate, not only `inside==true` ones.

**Independent Test**: After this story, the `[real_network]` ledger's `is-through-has-aps` line returns to `200 pass / 0 fail` even with the boundary exemption removed (T002).

### Implementation for User Story 0

- [X] T006 [US0] In `/workspace/src/mm/polymatch/polymatch_algorithm.cpp`, locate the three call sites of `record_inside(matched_point, layer_idx)` inside `build_hybrid_path`: the initial-element branch, the link→polygon transition branch, and the polygon→polygon same-polygon transition branch. Today each call is gated on the candidate's `inside` flag. Remove the `if (b->inside)` / `if (first_c->inside)` gates so `record_inside` is called unconditionally when the corresponding polygon Viterbi branch is taken. Each touched line gets an inline comment naming this commit (or the 003 feature in general) per Clarification Q2's bug-fix-citation rule.
- [X] T007 [US0] Rebuild `polymatch_test`. Re-run `./polymatch_test '[real_network]'` and verify the `is-through-has-aps` line is now `200 pass / 0 fail` with the boundary exemption already removed (T002). Capture the diff between `/tmp/003_red.txt` and the new output to `/tmp/003_us0_green.txt`.
- [X] T008 [US0] Re-run the default `./polymatch_test` invocation (no tag filter). Verify all synthetic-fixture tests still pass. Any specific-value assertion that now fails because the matcher's polygon-segment shape changed legitimately gets updated per Clarification Q2 — each touched assertion gets an inline comment like `// updated by specs/003-polymatch-bugfixes: matcher now sets is_through=false for any polygon Viterbi candidate`. `link-only-eq-weightmatch` MUST stay strict-green.

**Checkpoint**: US0 complete. `is-through-has-aps : 200 pass / 0 fail` with the strict invariant. Full polymatch_test still green. cpath-topology still red (US1 work).

---

## Phase 4: User Story 1 (Priority: P1) — Eliminate cpath-Topology Violations on Mid-Polygon-Start Trips

**Goal**: Trace 1313 and 1314's mid-polygon-start cpath emissions stop violating the polygon ↔ edge AP-incidence rule. The `cpath-topology` line returns to `200 pass / 0 fail` with the tightened `CHECK == 0` from T003.

**Independent Test**: After this story, the `[real_network]` ledger's `cpath-topology` line shows `200 pass / 0 fail` and the unit-test assertion `CHECK(cpath_fail == 0)` passes.

### Diagnose then fix (research.md R2)

- [X] T009 [US1] In `/workspace/src/mm/polymatch/polymatch_algorithm.cpp`, instrument `build_hybrid_path` with `SPDLOG_DEBUG` traces around every `push_link`, `push_polygon`, `set_entry`, `set_egress` call. Each trace logs the current `cpath.size()`, the value being pushed/recorded, and the source branch (initial / link→link / link→polygon / polygon→link / polygon→polygon). Used only for this diagnostic step; will be removed in T012.
- [X] T010 [US1] Hand-extract trace 1313 from the committed CSV into a one-row test CSV at `/tmp/one_1313.csv` and run `./polymatch --network ... --polygons ... --access_points ... --gps /tmp/one_1313.csv --output /tmp/out_1313.csv --log_level 1 -k 8 -r 300 -e 50 --weight cost --through_penalty_factor 1.5` (full args per quickstart.md "Diagnosing the cpath-topology bug"). Capture the matcher's debug log to `/tmp/003_t1313.log`. Cross-reference with polygon 28's AP node IDs in `access_points.shp` to identify the bad transition. Repeat for trace 1314 + polygon 172. Pick whichever hypothesis (R2 #1, #2, or #3) matches the evidence.
- [X] T011 [US1] Apply the fix indicated by T010's diagnosis to `/workspace/src/mm/polymatch/polymatch_algorithm.cpp`. The fix MUST be the minimum change that makes the bad transition emit an AP-incident edge; if the diagnosis points at a different branch than the polygon→link emit, document that in an inline comment.
- [X] T012 [US1] Remove the `SPDLOG_DEBUG` traces added in T009 — they were diagnostic only. Confirm the matcher's release-build log surface is unchanged (no leftover noise at INFO level).
- [X] T013 [US1] Re-run `./polymatch_test '[real_network]'`. Verify the `cpath-topology` line is now `200 pass / 0 fail` and the final `CHECK(cpath_fail == 0)` passes. Capture output to `/tmp/003_us1_green.txt`.
- [X] T014 [US1] Re-run the full default `./polymatch_test` (no tag). All synthetic tests still pass; per Clarification Q2 any specific-value assertions affected by the fix get inline-comment updates. `link-only-eq-weightmatch` stays strict-green.

**Checkpoint**: US1 complete. All four `[real_network]` invariants pass at 200/0. Default `polymatch_test` green. Traces 1313 + 1314 produce topologically valid cpath.

---

## Phase 5: User Story 2 (Priority: P2) — Diagnose the [real_network] Suite's Wall-Time Overshoot

**Goal**: Produce a 1-2 paragraph profile summary identifying the dominant cost contributor in the real-network suite. Apply *only* a single obvious quick win if one surfaces — defer larger optimization per Clarification Q1.

**Independent Test**: After this story, `specs/002-real-network-validation/tasks.md` contains an updated "Performance follow-up" entry with at least one measured number per phase (fixture load / per-trace match / invariant verification / etc.).

### Instrument + measure

- [X] T015 [US2] In `/workspace/test/polymatch_test.cpp`'s main `TEST_CASE("Real-network validation against committed trace batch ...")`, add scoped `std::chrono::steady_clock` accumulators around five phases per research.md R3: (a) RealAreaFixture construction (one-shot, measured once), (b) `load_real_network_trips()` call, (c) per-trace `POLYMATCH::match_traj` total, (d) per-trace `WEIGHTMATCH::match_traj` total, (e) per-trace invariant-check total. Emit the breakdown via `WARN(...)` so it appears in CI logs without failing the test.
- [X] T016 [US2] Run `./polymatch_test '[real_network]'` once with the instrumentation in place. Capture the breakdown output to `/tmp/003_perf.txt`. Form a 1-2 paragraph summary naming the dominant contributor (e.g., "per-trace POLYMATCH::match_traj averages X ms × 200 = Y s, dominating the wall time; harness loop overhead is Z s; fixture load is W s").
- [X] T017 [US2] Decide based on T016: does the profile reveal a SINGLE OBVIOUS quick win (e.g., a redundant per-iteration allocation, a misordered loop, an O(n²) check)? If yes, apply it now in `/workspace/src/mm/polymatch/polymatch_algorithm.cpp` or wherever the cost lives, then re-measure. If no, leave the matcher untouched.
- [X] T018 [US2] Append the profile summary (from T016, plus any quick-win delta from T017) to `/workspace/specs/002-real-network-validation/tasks.md` under "## Deferred Follow-Ups → ### Performance follow-up". Include the pre-fix and post-fix wall times for the suite if T017 applied a fix.
- [X] T019 [US2] Remove the `WARN` instrumentation timers from `/workspace/test/polymatch_test.cpp` — they were for the diagnostic measurement, not a permanent fixture (per quickstart.md "Performance diagnostic walkthrough").

**Checkpoint**: US2 complete. Perf profile recorded in 002's Deferred Follow-Ups. The wall-time decision (optimize vs revise budget) remains a deferred follow-up unless T017 applied a single fix.

---

## Phase 6: Polish & Non-Regression Sweep

- [X] T020 [P] Regenerate the committed real-network trace CSV with the canonical seed to prove the matcher fix doesn't bleed into the generator: `./polymatch_traces_gen --network ... --polygons ... --access_points ... --seed 2026 --output /tmp/trips_verify.csv`. `diff` against the committed `trips.csv`. The diff MUST be empty (FR-010 — generator unaffected).
- [X] T021 [P] Run all six test binaries from `build/test/`: `./algorithm_test && ./fmm_test && ./network_test && ./network_graph_test && ./weightmatch_test && ./polymatch_test`. All six MUST print "All tests passed" (constitution Principle III).
- [X] T022 [P] Update the `[real_network]` ledger expected-output block in `/workspace/specs/002-real-network-validation/quickstart.md` to show all four invariants at 200/0 pass (or whatever the post-fix steady state actually is), removing the pre-fix examples that showed `2 fail` for `cpath-topology`.
- [X] T023 [P] In `/workspace/specs/002-real-network-validation/tasks.md`, update the "Discovered matcher bugs" entry under Deferred Follow-Ups to record that traces 1313 + 1314 + the is-through-boundary-exemption issue are RESOLVED in 003, with a one-line pointer to this feature's spec.
- [X] T024 Manually walk `/workspace/specs/003-polymatch-bugfixes/quickstart.md` end-to-end against the freshly built executables. Confirm "Reproducing the failures" no longer reproduces them; "Verifying the fix (post-fix gate)" shows the expected green ledger; "Non-regression sweep" all passes.

---

## Dependencies & Execution Order

### Phase dependencies

- **Phase 1 (Setup)** → no dependencies; capture the baseline before anything changes.
- **Phase 2 (Foundational — harness tightening)** → after Phase 1. Blocks both user stories. **Suite goes red after Phase 2.**
- **Phase 3 (US0 — is_through fix)** → after Phase 2. Turns `is-through-has-aps` back green. cpath-topology still red.
- **Phase 4 (US1 — cpath-topology fix)** → after Phase 3. Turns `cpath-topology` back green. All four invariants at 200/0.
- **Phase 5 (US2 — perf diagnostic)** → after Phase 4 (matcher must be correct before timing it).
- **Phase 6 (Polish)** → after all user stories.

### Within-phase dependencies

- **Phase 2**: T002 (`check_is_through_has_aps` strict) and T003 (cpath tolerance) are different concerns but same file; sequence them T002 → T003 → T004 (unit-test update for the strict invariant) → T005 (build + verify red). T004 depends on T002. T005 depends on all three.
- **Phase 3**: T006 (matcher fix) → T007 (verify is-through-has-aps green) → T008 (default-suite non-regression). Strictly sequential.
- **Phase 4**: T009 (instrument) → T010 (diagnose) → T011 (fix) → T012 (un-instrument) → T013 (verify) → T014 (non-regression). Strictly sequential.
- **Phase 5**: T015 → T016 → T017 (conditional) → T018 → T019. Strictly sequential.
- **Phase 6**: T020, T021, T022, T023 are [P]. T024 is the final manual walkthrough.

### Parallel opportunities

- **Phase 6**: 4 polish tasks in parallel (T020-T023). T024 is the final serial walk.

---

## Implementation Strategy

### MVP scope

The MVP is **US0 + US1 combined** — fixing the two correctness bugs the real-network suite identified. US2 (perf diagnostic) is a P2 cleanliness deliverable that can land later if scope pressure arises.

### Incremental delivery path

1. Phase 1 (baseline capture).
2. Phase 2 (harness goes red — TDD anchor).
3. Phase 3 (US0 — is_through green).
4. Phase 4 (US1 — cpath-topology green). **Internal milestone — correctness MVP**.
5. Phase 5 (US2 — perf diagnostic).
6. Phase 6 (polish + non-regression).

### Single-developer sequencing

T001 → T002 → T003 → T004 → T005 (Phase 2 red anchor) → T006 → T007 → T008 (Phase 3 US0 green) → T009 → T010 → T011 → T012 → T013 → T014 (Phase 4 US1 green) → T015 → T016 → T017 → T018 → T019 (Phase 5 US2 diagnostic) → (T020, T021, T022, T023 parallel) → T024.

---

## Notes

- **TDD discipline** per Constitution III: T002-T005 deliberately go red before any matcher change. Skipping this would let a buggy fix appear to "pass" against the previously-loose assertions.
- **No new files.** Every task edits an existing source file or appends to an existing markdown file.
- **`trips.csv` stays untouched** (FR-010). T020 verifies this by re-running the generator and diffing.
- **Synthetic-test specific-value assertions** that legitimately change due to the fix get inline-comment updates per Clarification Q2. `link-only-eq-weightmatch` is the one gate that stays strict — a fix breaking it is rejected outright.
- **US2 perf optimization is OUT of scope** unless T016's profile reveals a single obvious quick win. Otherwise the perf decision (optimize vs revise SC-001) stays a deferred follow-up in 002.
- Commit by phase. Phase 2 produces a "harness goes red" commit (intentionally failing CI — note this in the commit message). Phase 3 and Phase 4 each produce a "matcher fix" commit that returns CI to green.
