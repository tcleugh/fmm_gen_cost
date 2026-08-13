# Brainstorm: Event-level GPS accuracy in map matching

**Date:** 2026-06-01
**Status:** Direction converged (2026-06-01); **implementation-detail revisit 2026-08-13** — see
[Implementation addendum](#implementation-addendum--2026-08-13-code-grounded-revisit) for the
code-grounded plan (Task 1 / Branch C of
[scopestorms/260623-0117-road-traffic-update-revisit.md](../scopestorms/260623-0117-road-traffic-update-revisit.md)).
**All revisit open decisions resolved 2026-08-13**; only the concrete parameter *values*
(`sigma_floor`/`sigma_cap`/`k_max`) remain, pending the car accuracy-distribution parquet — see the
[parameter-setting protocol](#accuracy-distribution-data-car--parameter-setting-protocol-2026-08-13).

## Goal

Replace (or augment) the single global `gps_error` parameter in map matching with
**per-event (per-point) reported GPS accuracy**, so the emission probability reflects
the actual reported uncertainty of each individual GPS fix rather than one fixed value
for the whole dataset.

**Motivation (confirmed — framing A):** Heterogeneous trace quality. Traces mix clean
open-sky fixes with noisy urban-canyon/tunnel fixes. A single global sigma is a bad
compromise: set it low and the noisy points snap confidently to the wrong road; set it
high and the clearly-good points lose their pull, letting transition probability drag the
match off correct snaps. Per-point sigma lets each fix contribute weight proportional to
its actual reported reliability.

**Target matcher (behaviour):** WeightMatch (`src/mm/weightmatch/`). The *I/O + data model*
(accuracy field on `Trajectory`, reader parsing) is **shared infrastructure** and is added
across all readers for consistency. But the *use* of per-point accuracy in matching stays
WeightMatch-only — see "Shared-structure nuance" below.

### Shared-structure nuance
`TransitionGraph` and the `Network` candidate-search functions are shared by FMM, STMatch
and WeightMatch. To keep FMM/STMatch behaviour unchanged:
- `TransitionGraph` constructor changes from a scalar `gps_error` to a **per-point sigma
  vector**; FMM/STMatch fill it with the global value replicated (no behaviour change),
  WeightMatch fills it from reported accuracy.
- The sigma-band candidate admission is wired into **WeightMatch's** search path only;
  shared search functions stay backward-compatible (current callers reproduce today's
  behaviour).

### How it works today (baseline)
- Emission probability: `calc_ep(dist, error) = exp(-0.5 * (dist/error)^2)` —
  a zero-mean Gaussian where `error` = global `gps_error` (sigma, in map units),
  `dist` = point-to-candidate-road distance. (`transition_graph.cpp:39`)
- `gps_error` is set once in config and passed to the `TransitionGraph` constructor
  (`transition_graph.cpp:19`), applied uniformly to every point of every trajectory.
- `Trajectory` struct (`core/gps.hpp:26`) holds `id`, `geom`, `timestamps` only —
  **no per-point accuracy field**, and the GPS reader does not parse one.
- WeightMatch builds the graph at `weightmatch_algorithm.cpp:170`.

### Layers this touches
1. **Data model** — add per-point accuracy to `Trajectory` / candidate structures.
2. **I/O** — parse reported accuracy from GPS input (CSV / shapefile / API).
3. **Math** — `calc_ep` uses a per-point sigma instead of the global constant.

## Options Explored

The plumbing (CSV accuracy column → parallel `std::vector<double>` on `Trajectory` →
threaded into candidate search + `TransitionGraph`) is common to every approach and is
mechanical. The real design choices are (a) how reported metres become each point's
emission sigma, and (b) how to stop the correct candidate being pruned for noisy points.

### Emission-sigma mapping (how reported metres → sigma in `calc_ep`)

- **1 — Direct drop-in:** `sigma_i = reported_i`. Smallest change. *Risk:* trusts the
  reported metres as literally 1σ despite unknown confidence → possible systematic bias;
  no guard against degenerate values. **Verdict:** Discarded as the whole answer — too
  fragile given known ep numerical edge cases.
- **2 — Calibrated sigma:** `sigma_i = clamp(scale · reported_i, floor, cap)`. Tunable
  `scale` absorbs the unknown-confidence problem; floor/cap guard degenerate feed values.
  **Verdict:** Kept — forms the emission half of the chosen direction.
- **4 — Quadrature blend:** `sigma_i = sqrt(base² + reported²)`. Never collapses to zero.
  **Verdict:** Parked — viable floor mechanism, but `clamp` floor is simpler/clearer.
- **5 — Outlier down-weighting only:** keep global sigma, just clamp bad points.
  **Verdict:** Discarded — leaves most of the signal unused; doesn't deliver framing A.

### Candidate-gating fix (the binding constraint — confirmed by user)

The real failure isn't search radius — it's **`k`**. `filter_candidates`
(`network.cpp:417`) keeps only the `k` nearest edges within radius. In a dense network a
noisy point's true edge is *within radius* but ranked past `k` by distance, so it is
pruned before the HMM sees it.

- **3 — Coupled sigma + adaptive search (CHOSEN family):** derive per-point candidate
  admission from reported accuracy so the true edge survives filtering. Two sub-options
  were considered:
  - **(i) scale `k`** — keep "k-nearest", vary `k` per point. *Discarded:* still pure
    distance ranking, needs an invented `g(sigma)`, and `k` is the expensive lever.
  - **(ii) sigma-scaled distance band (CHOSEN):** admit edges within `band_i = c·sigma_i`,
    so the *same* per-point sigma drives both emission and admission.

**Candidate rule (final form of (ii), incorporating the min-count guard):**
Among all edges within an outer bound `R = max_radius`, the candidate set is **every edge
within `band_i = c·sigma_i`** (a distance *threshold*, NOT a top-k — keep them all), then
bounded only by count: **never fewer than `k` nearest** (expand search out to `R` to meet
it) and **never more than `k_max` nearest**. Inside the band there is no distance ranking;
ranking applies only when the floor/ceiling bounds bite.

| Situation | band | Result |
|---|---|---|
| Noisy point | wide | all within band, capped at `k_max` nearest |
| Clean point, dense net | tiny, ≥ `k` inside | keep those (tight, fast) |
| Clean point, sparse net | tiny, < `k` inside | expand to `k` nearest within `R` (the 2 m-accuracy case) |
| Noisy + very dense | very wide, > `k_max` | cap at `k_max` nearest (must log — silent truncation) |

This **supersedes** the earlier `r_i = max(r, c·sigma)`: the floor is now *count-based*
(`k` candidates) not *distance-based* (`r` metres), which is the only thing that
guarantees coverage for a clean point in a sparse network. It generalizes FMM's existing
`search_tr_cs_knn_with_fallback` (`network.cpp:335`) from "expand if empty" to
"expand until ≥ `k`".

**Verdict:** Kept — Approach 3 / sub-option (ii) with count-based floor is the design.

## Risks & Pitfalls

- **`k` drives cost super-linearly.** WeightMatch routes shortest paths between candidate
  pairs across consecutive layers → roughly O(k²) shortest-path computations per layer
  transition. Raising `k` is more expensive than raising `r` (which only affects the
  rtree query + filtering). A point at `k=32` is ~16× the transition cost of `k=8`.
  → **A hard `k_max` cap is load-bearing**, and per-point scaling is good *because* it
  concentrates the cost only on the genuinely noisy points.
- **How fast must `k` grow with sigma?** To "reach" a true edge at distance ~sigma in a
  network of edge-density ρ, the number of nearer decoy edges scales with *area* ~ ρ·sigma².
  So in the worst case the required `k` grows with **sigma², not sigma** — which is exactly
  why an unbounded scaling is dangerous and a cap is mandatory. Open question: linear vs.
  quadratic `g(sigma)` in practice.
- **Degenerate reported values** (0, negative, NaN, absurdly large) will appear in real
  feeds — recent commits already fought ep `inf`/crash cases. Floor/cap on sigma and cap
  on `k`/`r` must defend against these.
- **Threading parallel vectors.** accuracy must stay index-aligned with points through the
  `matchable` trimming (`first_index`/`last_index` slicing at `weightmatch_algorithm.cpp`)
  — same hazard the timestamp handling already has.

## Preferred Direction

**Approach 3 / sub-option (ii)**, driven entirely by the per-point reported accuracy
column. One calibrated `sigma_i` per point drives **both** candidate admission and
emission weight:

1. **Ingest (all readers, for consistency):** add a per-point reported-accuracy
   `std::vector<double>` on `Trajectory`, **optional/empty when absent** — mirroring
   exactly how optional `timestamps` already work. Populate it in:
   - **CSVTrajectoryReader** — parallel comma-list column (e.g. one `accuracy` value per point).
   - **CSVPointReader** — a per-row accuracy column, grouped per trajectory.
   - **GDALTrajectoryReader (shapefile)** — a new field (comma-list per feature, like its timestamp field).
   - **Programmatic / Python** — optional accuracy arg on the `Trajectory` constructor.
   - **`match_wkt`** — geometry-only; stays as-is and uses the global-`gps_error` fallback.
   - **`GPSConfig`** gains an accuracy column/field name + presence detection.
   Keep index alignment through the `matchable` trimming.

   **Naming (resolved 2026-08-13):** the per-point field holds the *raw reported value* (always in
   metres, by contract), distinct from the global `gps_error` param and from the derived `sigma_i`.
   Column/field name is **`accuracy`** (not `gps_error`, to avoid colliding with the existing param).
2. **Per-point sigma:** `sigma_i = clamp(scale · reported_i, sigma_floor, sigma_cap)`.
3. **Candidate search (the fix that matters):** generalize
   `search_tr_cs_knn_with_fallback` to admit **all** edges within `band_i = c·sigma_i`,
   bounded by count to `[k, k_max]` (expand to `R = max_radius` to meet the `k` floor;
   trim to `k_max` nearest at the ceiling).
4. **Emission:** `calc_ep(dist, sigma_i)` per point instead of global `gps_error`.
5. **Fallback / backward-compat:** a point with missing/invalid accuracy uses the global
   `gps_error` for sigma and the global `k`/`radius` for search → **old CSVs without the
   column behave exactly as today.**
6. **Precondition (metric CRS):** input data is required to be in a metric CRS (metres).
   Document this explicitly, and add a startup guard that warns when coordinates look like
   lat/lon degrees (|x|≤180, |y|≤90) — a mistake that would silently corrupt every sigma.

### Config surface (proposed)
| Param | Role | Proposed default |
|---|---|---|
| `gps_error` (existing) | emission sigma *fallback* when accuracy missing | 50 (unchanged) |
| `radius` + `backup_radius` (existing) | consolidate into `R` — single hard outer search bound (`R = max(radius, backup_radius)`) | 300 (unchanged) |
| `k` (existing) | **min** candidate count (floor; now triggers graded expansion to `R`) | 8 (unchanged) |
| `backup_k` (existing) | **dropped** — folds into `[k, k_max]` | — |
| `allow_truncation` (existing) | terminal policy when even `R` finds nothing | unchanged |
| `scale` (new) | reported metres → sigma calibration | 1.0 (treat reported as 1σ) |
| `c` (new) | band width in sigmas | 3.0 (≈ Gaussian 3σ ⇒ emission ≥ e^−4.5; tighter clips plausible edges) |
| `k_max` (new) | **max** candidate count (cost backstop) | e.g. 24 — tune |
| `sigma_floor`/`sigma_cap` (new) | guard degenerate reported values (0/neg/huge) | tune |

Defaults chosen so a dataset with **no** accuracy column reproduces today's behaviour.

## Open Questions

- [x] Motivation = framing A (heterogeneous trace quality).
- [x] Accuracy is present in the user's raw source data (but NOT yet in FMM's data model / readers).
- [x] Reported in **metres, no stated confidence level** → mapping to sigma is a design choice.
- [x] Ingestion path = **CSVTrajectoryReader** (add parallel accuracy comma-list column).
- [x] Sigma mapping = `clamp(scale · reported, sigma_floor, sigma_cap)` (calibrated, guarded).
- [x] Backward-compat = missing accuracy → global `gps_error`/`k`/`radius` (today's behaviour).
- [x] Candidate gating = sigma-scaled distance band, count-bounded `[k, k_max]`, not k-nearest.
- [x] Radius adaptation = replaced by count-based floor (`k` candidates up to `R`).
- [x] I/O scope = **all readers** (CSV trajectory, CSV point, shapefile, programmatic) for
      consistency; `match_wkt` stays geometry-only with global fallback.
- [x] Behaviour scope = WeightMatch only; FMM/STMatch unaffected via backward-compatible
      shared `TransitionGraph` / search signatures.
- [x] **Units / CRS (resolved):** input CRS is always metres → reported metres are directly
      comparable to map units, `scale≈1` is meaningful, no conversion needed.
      **New requirement:** the metric-CRS assumption (previously implicit in the
      `gps_error=50`/`radius=300` defaults) is now load-bearing. Must be (a) clearly
      documented as a precondition, and (b) ideally guarded with a startup warning if
      coordinates fall in lat/lon ranges (|x|≤180, |y|≤90), which signals degrees by mistake.
- [x] **Min-`k` expansion vs. existing fallback (resolved):** two distinct concepts, not
      three. (1) **Outward expansion** — `radius`/`backup_radius`/`max_radius` are all the
      same role; consolidate to a single outer bound `R` (impl: `R = max(radius,
      backup_radius)` for config compat). The band + min-`k` graded expansion *subsumes*
      the old binary backup tier: "expand when `< k`" generalizes "expand when `== 0`".
      (2) **Terminal policy** — `allow_truncation` is orthogonal and unchanged; it governs
      what happens when even `R` finds nothing (drop unmatchable ends for long-driveway /
      off-network cases). Expansion feeds into it, never replaces it.
      Knock-ons: `backup_k` folds into the `[k, k_max]` count bounds (drop it);
      `backup_radius` folds into `R`.
      **Behaviour change to test:** expansion now triggers on "fewer than `k`", not only
      on "zero candidates" — a point that previously accepted 3 candidates now expands to
      reach `k`.
- [ ] Concrete defaults for `c`, `scale`, `k_max`, `sigma_floor`/`sigma_cap` (tune empirically).
- [ ] **Residual risk (accepted):** when `k_max` cap bites (noisy point in very dense net),
      the true edge can still be pruned. Mitigate with generous `k_max` + a log line when it bites.
- [ ] Field/column naming: `accuracy` vs `reported_accuracy` (avoid clash with `gps_error`).
- [x] Verify nothing else reads the global `gps_error` in the WeightMatch generalized-cost path.
      **Resolved (2026-08-13):** in WeightMatch, `config.gps_error` is consumed **only** at the
      `TransitionGraph` construction (`weightmatch_algorithm.cpp:170`). `update_layer` reads the
      pre-computed `ep` off each node, never `gps_error`; the generalized cost uses edge
      `weight`/`length`, not `gps_error`. So per-point sigma need only reach the TG constructor.

## Implementation addendum — 2026-08-13 (code-grounded revisit)

Revisited to bring the converged design to implementation-grade detail. Verified against the live
tree; the design holds. Confirmed line references: `calc_ep` (`transition_graph.cpp:39`),
scalar-`gps_error` TG constructor (`:19`), `Trajectory` + optional `timestamps` (`gps.hpp:26`),
`filter_candidates` top-k (`network.cpp:417`), `search_tr_cs_knn_with_fallback` (`network.cpp:335`),
graph build (`weightmatch_algorithm.cpp:170`).

### Locked decisions (2026-08-13)

- **D1 — Strict gate (backward-compat).** No-accuracy WeightMatch stays **byte-identical**. A
  trajectory with no/empty accuracy vector keeps calling today's `search_tr_cs_knn_with_fallback`
  and today's scalar-`gps_error` emission — existing WeightMatch goldens are frozen and must not
  move. Only accuracy-bearing trajectories route into the new sigma-band search + per-point
  emission. Two search paths coexist, branched in `match_traj`. This resolves the doc's (a)/(b)
  tension: the "expand on fewer-than-`k`" behaviour change is confined to the accuracy path and
  never touches no-accuracy trips. FMM/STMATCH untouched (separate `search_tr_cs_knn`).
- **D2 — I/O scope = CSV only.** Car matcher input is CSV. Wire accuracy into `CSVTrajectoryReader`
  (mode 1, comma-list column mirroring `string2time`) and `CSVPointReader` (mode 2, per-row
  column). **Deprioritise `GDALTrajectoryReader`** — it does not even wire `timestamps` today (it
  resolves the index only to answer `has_timestamp()`, then returns a 2-arg `Trajectory`), so
  there is no timestamp path to mirror; leave shapefile accuracy as a documented gap. `match_wkt`
  stays geometry-only (empty accuracy → global fallback). Supersedes the doc's "all readers" scope.

### Refined architecture (consequences of D1)

- **Two-path candidate search.** New WeightMatch-only `search_tr_cs_knn_sigma_band(...)`
  implementing the band + `[k, k_max]` count-bound rule; existing `search_tr_cs_knn_with_fallback`
  is left untouched (serves no-accuracy trips; FMM/STMATCH use their own function). Branch in
  `match_traj` on `traj.accuracy` presence.
- **Overload the TG constructor, don't replace it.** Keep `TransitionGraph(tc, double gps_error)`
  verbatim (single shared definition at `:19`, also used by FMM/STMATCH); add
  `TransitionGraph(tc, const std::vector<double>& sigmas)`. FMM/STMATCH literally unchanged.
- **Single slice, no re-expansion.** Emission sigma is consumed only at TG construction, never
  surfaced in `MatchResult`. So compute `sigma_full[0..N-1]` once (the band search consumes it over
  the full `traj.geom`, pre-trim), then slice `sigma_trimmed[i] = sigma_full[i + first_index]`
  alongside the `tc` slice at `:162-165`. The re-expansion block (`:236-260`) needs no change —
  simpler than the doc's "thread it through everything" framing.
- **Config = additive (supersedes the doc's consolidation).** D1 makes the proposed
  `R = max(radius, backup_radius)` consolidation and `backup_k` drop counterproductive: the
  untouched old path still needs `backup_k`/`backup_radius` verbatim. Keep all existing params;
  add `scale`/`c`/`k_max`/`sigma_floor`/`sigma_cap`. The band search's outer bound reuses
  `max(radius, backup_radius)` **without renaming config**.

### Adjacent finding — `eu_dists` misalignment through trimming

`match_traj` trims candidates to `[first_index,last_index]` (`:162-165`) so
`layers.size()==matchable_count`, but `update_tg` computes `eu_dists = cal_eu_dist(traj.geom)` over
the **full** geometry (`:273`) and indexes it with the *trimmed* layer index. When `first_index>0`
(leading points dropped under `allow_truncation`), `eu_dists[i]` is the gap for full-point `i`, not
`first_index+i` — misaligned. Pre-existing, independent of this work, reachable only with
`allow_truncation` + missing leading candidates. The sigma vector must slice as
`sigma_full[i+first_index]` to be correct; the same discipline would fix `eu_dists`.
**Resolved 2026-08-13: fix in this pass** — thread `first_index` into `update_tg` and index
`eu_dists[i + first_index]` (the same slice-at-`first_index` discipline as the sigma vector).
Regenerate any golden that exercises a leading trim (`first_index > 0`).

### Per-point fallback — heterogeneous admission (resolved 2026-08-13)

Decision #2 (user): a bad per-point value falls back to **global `k`/`radius` admission + global
`gps_error` emission** for that point — *not* clamp-only. So the accuracy path is **per-point
heterogeneous** within one trajectory:

| Point | Valid (`finite(reported) && reported > 0`) | Invalid (sentinel `-1`, `0`, `<0`, `NaN`/`inf`) |
|---|---|---|
| Admission | band `c·sigma_i` within `R`, count-bounded `[k, k_max]`, graded expansion to reach `k` | today's `filter_candidates` top-`k` within `radius` (+ backup tier) — **exact legacy** |
| Emission `sigma` | `clamp(scale·reported_i, floor, cap)` | global `gps_error` |

Candidate search is independent per point, so an invalid point gets a **byte-identical** candidate
set and emission to today regardless of its neighbours' rules. Consequences:
- `search_tr_cs_knn_sigma_band` takes the raw per-point `reported` vector (or pre-computed `sigma` +
  a validity flag) and branches per point: valid→band, invalid→the existing per-point legacy routine
  (reuse `get_all_candidates` + `filter_candidates` + backup, already factored in `network.cpp`).
- The emission sigma vector passed to the overloaded TG constructor is
  `sigma_emission[i] = valid_i ? clamp(scale·reported_i, floor, cap) : gps_error`, then sliced
  `[first_index,last_index]`.
- **Degenerate-but-correct:** an accuracy-bearing trajectory whose points are *all* invalid
  reproduces today's behaviour exactly (every point on the legacy branch), just routed through the
  new function.

**Sentinel = `-1` (recommended; user sets it upstream).** Physically impossible as a real accuracy
(metres ≥ 0) so unambiguous; parses trivially with the existing `>> double` comma-list reader (no
`NaN`-text fragility, no empty-field desync); matches the repo's `-1 = unset` idiom
(`backup_k`/`backup_radius`, fake-candidate `dist`/`offset`); and unifies the sentinel with every
degenerate value under the single `finite && > 0` predicate. Inferred readings are replaced with
`-1` upstream to put them on the fallback path.

### Open decisions (this revisit)

- [x] **Sigma validity + per-point fallback (resolved 2026-08-13).** Sentinel `-1`; predicate
      `finite(reported) && reported > 0`; invalid → global `gps_error` emission + global `k`/`radius`
      admission (see *Per-point fallback* above).
- [~] Concrete defaults for `scale`/`c`/`k_max`/`sigma_floor`/`sigma_cap` + tuning approach
      (car global error 150 m vs walk 25 m). **Protocol defined 2026-08-13** — see
      [Accuracy-distribution data (car)](#accuracy-distribution-data-car--parameter-setting-protocol-2026-08-13);
      values pending the parquet bundle.
- [x] Column naming (resolved 2026-08-13): **`accuracy`**, value **always in metres** by contract → no
      column-unit config or conversion; directly comparable to map units under the metric-CRS precondition,
      so `scale ≈ 1` is meaningful. (Absent-column fallback still covers legacy/non-car inputs per D1.)
- [x] Tunnel-portal events (resolved 2026-08-13): the pipeline writes a **small explicit `accuracy` value
      (metres)** on synthetic portal points → they take the valid band + tight-emission path; **no C++
      confidence flag, no Task 1 special-casing**. Owned by **Task 2**. Task 2 notes: set the value
      ≥ `sigma_floor` (else it's silently clamped — still tight); the count-floor `k` guarantees the portal
      edge is in the candidate set, so pinning comes from emission weight, not a starved candidate set.
- [x] `eu_dists` misalignment (resolved 2026-08-13): **fix in this pass** — `eu_dists[i + first_index]`
      via `first_index` threaded into `update_tg`; regenerate any golden with a leading trim.
- [x] Metric-CRS guard (resolved 2026-08-13): at `network.cpp:~135` inside the existing
      `ogrsr != nullptr` block — **`ogrsr->IsGeographic()` → `SPDLOG_WARN`** (authoritative CRS-metadata
      check; **supersedes** the doc's coordinate-range heuristic). Range check (`|x|≤180, |y|≤90` over the
      network bbox) kept only as the `ogrsr == nullptr` fallback. Warn, not hard-error. Residual (unguarded):
      a projected non-metre CRS, e.g. US survey feet — rare vs the degrees mistake.
- [x] `has_accuracy()` predicate (resolved 2026-08-13): new `has_accuracy()` uses `accuracy_idx >= 0`.
      **Also fix the existing `has_timestamp()` `> 0` → `>= 0`** in all three readers
      (`gps_reader.cpp:90, 161, 331`) — same latent bug (a column at index 0 is silently dropped;
      reachable in GDAL where field 0 is arbitrary). Done as part of Task 1.

## Accuracy-distribution data (car) — parameter-setting protocol (2026-08-13)

The car trip-extraction repo is producing an 8-output parquet bundle (car-only for now; 2 runs) as
the empirical basis for the tuning open question. It is a **marginal + provenance + spatial +
engagement** view of reported accuracy — deliberately **not paired with true error** — which is
exactly why it settles the *guard rails and cost backstop* but not the *calibration* (`scale`).

### The eight outputs (as delivered)
1. **`trip_windows`** — reconciliation baseline (`n_trips`, `n_trip_points`, first/last start). *Gate.*
2. **`point_provenance`** — per `point_source` × run: `n_points`/`n_trips`/`n_valid`/`n_null`/`median_m`/`p90_m`.
   Inferred-vs-measured composition + correctness gate (`unmatched`, `provenance_unset` must be ~0).
3. **`accuracy_percentiles`** — per run × {`all`,`measured`,`implied_slow_seg`}: full `p1…p99_9` ladder +
   min/mean/max. **Valid only** (`NOT NULL && NOT NaN && >0`). The marginal distribution.
4. **`degenerate_counts`** — raw counts: `n_null`, `n_null_on_matched_event`, `n_nan`, `n_zero`,
   `n_negative`, `n_integer_valued`, `n_under_5m`, `n_over_150m`/`300m`/`1000m`, `n_measured`,
   `n_implied_slow_seg`, `n_tunnel_synthetic`, `n_unmatched`, `n_provenance_unset`, `n_duplicate_event_rows`.
5. **`value_spikes`** — top-30 exact values × {all,measured,implied} (**keeps 0/neg/NaN**). A spike at a
   round constant = a filled-in default → sets sentinel policy.
6. **`per_trip_invalid_summary`** — `point_weighted_invalid_frac`, `mean_trip_invalid_frac`, per-trip
   invalid-frac percentiles, `n_trips_fully_valid`/`majority_invalid`/`fully_invalid`. The engagement measure.
7. **`per_trip_invalid_histogram`** — 10 buckets of per-trip invalid frac (top bucket `[0.9,1.0]` closed).
8. **`accuracy_by_geography`** — per run × {SOSR, GCCSA, SA3, SA2}: `median_m`/`p90_m`/`p95_m`/
   `invalid_frac`/`share_over_150m`. Per-fix geography by meshblock. SA3/SA2 restricted to Greater
   Melbourne + Greater Sydney, ≥5,000 pts/area.

### Which output sets which default
| Default | Source(s) | Rule |
|---|---|---|
| `sigma_floor` | #3 `measured` `p1_m`/`p5_m`/`min_m` · #4 `n_under_5m` · #5 low spikes | Floor at/near measured `p1_m` (or a physical few-m floor); `n_under_5m` sizes the lift |
| `sigma_cap` | #3 `measured` `p99_m`/`p99_9_m`/`max_m` · #4 `n_over_300m`/`n_over_1000m` · #5 high spikes | Cap ≈ measured `p99_9_m`; over-1000 m count = tail clipped |
| `k_max` | #3 tail · #4 `n_over_150m`/`300m` · #8 `share_over_150m` (SA2/SA3) | Distribution scopes *where/how often* bands go wide; **final cap = a network band-crowding run** on dense inner-city, not a percentile |
| validity predicate | #4 degenerate breakdown · #5 `value_spikes` (**QA only**) | `finite && >0`, sentinel `-1`; the fills are **inferred points sentineled upstream by provenance** — no C++ magnitude/share rule; #5 just verifies the `measured` class is clean |
| engagement | #6 `point_weighted_invalid_frac`,`n_trips_fully_invalid` · #7 histogram · #2 provenance | High invalid frac / many fully-invalid trips → per-point mostly collapses to the global fallback |
| tunnel (Task 2) | #2 provenance `n_tunnel_synthetic` | Synthetic points carry no accuracy → explicit tight sigma, not loose fallback |

Notes: floor/cap read off the **`measured`** split *only* — every `implied_slow_seg` point is sentineled to
`-1` upstream (§ *Inferred points* below), so that split never contributes a sigma and is diagnostic-only.
#8 is the empirical test of framing A
(does accuracy degrade in the dense inner-city where the parallel-road error bites) *and* localizes where
`k_max` cost concentrates. #1/#2 are gates — don't trust a percentile until `unmatched`/`provenance_unset` ~0.

### What this data cannot set (residual gaps)
- **`scale`** — needs (reported, true-error) pairs; the marginal has none. In-house proxy = reported vs
  point-to-matched-road residual from a baseline match. **Resolved 2026-08-13: accept `scale = 1.0`**
  (reported treated as 1σ) and tune it downstream on match quality (Task 2, per-state, no-regression);
  the paired-residual extract is **not** commissioned now.
- **`k_max` final value** — the distribution bounds *how often* the band is wide; the actual candidate
  count within `c·sigma` is a network property → a matcher-side experiment on dense SA2s (Task 2 tuning).
- **`c`** — stays an a-priori 3σ coverage choice; cost-validated by the `k_max` run, not the marginal.

### Open (this data)
- [x] **`scale` calibration (resolved 2026-08-13):** accept `scale = 1.0` for now, tune on match quality
      downstream (Task 2). Paired `reported_accuracy` vs `baseline_match_residual_m` extract declined —
      not worth commissioning now.
- [x] **Sentinel predicate vs default-fill (resolved 2026-08-13):** the "default-fill" population *is* the
      inferred points (slow-seg duplicated/inherited accuracy), known by provenance and **sentineled to `-1`
      upstream** (reaffirms the locked "inferred → `-1`" decision). C++ predicate stays `finite && >0` with no
      magnitude/share heuristic. #5 `value_spikes` demotes to a **QA check** on the `measured` class — a
      residual round-constant spike there = an unknown *device-vendor* default (distinct from our inference)
      → revisit the heuristic only then. Consequence: slow-seg points route to the global fallback (car
      `gps_error` ≈ 150 m, loose) — the conservative choice vs trusting one seeding fix's tight value across a dwell.

## Discarded Ideas

| Idea | Reason discarded |
|---|---|
| _(none yet)_ | |
