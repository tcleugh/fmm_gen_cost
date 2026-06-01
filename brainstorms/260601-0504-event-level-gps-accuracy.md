# Brainstorm: Event-level GPS accuracy in map matching

**Date:** 2026-06-01
**Status:** Complete (direction converged; defaults to tune in implementation)

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

   **Naming:** the per-point field holds the *raw reported value* (metres), distinct from
   the global `gps_error` param and from the derived `sigma_i`. Name it `accuracy` /
   `reported_accuracy` (not `gps_error`) to avoid colliding with the existing param.
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
- [ ] Verify nothing else reads the global `gps_error` in the WeightMatch generalized-cost path.

## Discarded Ideas

| Idea | Reason discarded |
|---|---|
| _(none yet)_ | |
