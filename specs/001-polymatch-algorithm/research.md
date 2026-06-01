# Phase 0 Research: PolyMatch

This document resolves the open technical questions in the implementation plan. Each section captures
a **Decision**, the **Rationale**, and **Alternatives considered**.

---

## R1. Polygon-aware routing graph design

**Decision**: Introduce a new class `FMM::ROUTING::PolyLinkGraph` that owns a generalized vertex ID
space:
- Vertices `[0, |E|)` map to road `EdgeIndex` (same as existing `LinkGraph`).
- Vertices `[|E|, |E|+|P|)` map to polygon vertices (one per polygon).

Adjacency is built once at construction:
- **link→link** edges copied from the existing `LinkGraph`.
- **link→polygon** edges: from a road edge whose endpoint coincides with an access point of polygon P,
  add an arc to P's polygon vertex with cost determined at routing time (depends on the matched
  point inside P).
- **polygon→link** edges: symmetric to the above.
- **polygon→polygon** edges: when an access point's node number references both P and Q in the access
  point shapefile, add arcs P↔Q.

The existing Dijkstra core (`shortest_edge_to_edges` in `link_graph_routing.cpp`) is generalized into a
template over the vertex type, or copied to a sibling function `shortest_polylink_to_polylinks` that
operates on `PolyLinkGraph`. The same `DijkstraState` / `IndexedMinHeap` are reused (must be sized for
`|E|+|P|`).

**Rationale**: Treating polygons as additional vertices in the same routing graph keeps the algorithmic
core (Dijkstra with reusable state) intact, mirrors the existing edge-as-vertex convention, and lets
polygon traversal costs be expressed as arc weights — no special-case logic in the inner loop.

**Alternatives considered**:
- *Modify the existing `LinkGraph` in place*: Rejected — violates user constraint "Existing functionality
  of other matching types should not be changed" and risks breaking WeightMatch / STMATCH.
- *Two-pass matching (link Dijkstra, then polygon postprocessing)*: Rejected — would produce incorrect
  HMM transition probabilities because polygon segments must compete with link-only paths during
  Viterbi, not after.
- *Polygon as a "fat" edge in `LinkGraph`*: Rejected — a polygon may have many access points, so it
  cannot be modeled as a single edge; modeling as a complete clique of access-point pairs would inflate
  edge count to O(AP²).

---

## R2. Polygon shapefile loading via GDAL

**Decision**: Mirror `FMM::NETWORK::Network::read_ogr_file()` pattern. New class
`FMM::NETWORK::PolygonLayer` opens the shapefile via `GDALOpenEx(..., GDAL_OF_VECTOR | GDAL_OF_READONLY)`,
iterates features, parses `id` (configurable field name) and optional `cost` field, validates
`wkbPolygon` / `wkbMultiPolygon` geometry, and stores polygons as Boost.Geometry `polygon<point_2d>`
for fast point-in-polygon and distance-to-boundary queries via Boost.Geometry's R-tree.

**Rationale**: Same library and idiom already used everywhere in `FMM::NETWORK`. Boost.Geometry
polygons provide `boost::geometry::within()` (point-in-polygon), `boost::geometry::distance()` (point to
boundary), and `boost::geometry::index::rtree` (spatial index) — all O(log N) for typical cases. No new
dependencies.

**Alternatives considered**:
- *Use raw GDAL geometry types throughout*: Rejected — GDAL's geometry primitives are heavier and
  inconsistent with the rest of the codebase, which uses Boost.Geometry internally.
- *Custom polygon class*: Rejected — Boost.Geometry is already a dependency and is battle-tested.

---

## R3. Access point shapefile loading

**Decision**: New class `FMM::NETWORK::AccessPointLayer` loads the access point shapefile, expecting
features with:
- Point geometry (`wkbPoint`).
- A node-number field (configurable, default `node_id`).
- A polygon-ID field (configurable, default `polygon_id`).

Loading produces:
- A `std::vector<AccessPointFeature>` (one entry per feature in the shapefile).
- A map `polygon_id → vector<access_point_node>` for fast lookup during routing-graph construction.
- A map `node_id → vector<polygon_id>` for resolving polygon-to-polygon shared access points.
- A spatial R-tree of unique access point points for candidate-radius queries.

**Validation pass** (executed during load; halts with descriptive error per FR-005):
1. For each feature, check that its point geometry lies on the boundary of its declared polygon
   (`boost::geometry::distance(point, polygon_boundary) <= EPSILON`, default `1e-6` map units).
2. For each feature, check that its polygon_id exists in the loaded `PolygonLayer`. Reject orphans.
3. Group features by `node_id`; within each group, verify all geometries are identical (within EPSILON).
   Reject contradictory shared access point definitions.

**Rationale**: Aligns with the spec's three FR-005 validation conditions. The validation work happens
once at load time, before any matching — cheap to do, expensive to skip.

**Alternatives considered**:
- *Lazy validation during matching*: Rejected — fails the spec's "halt with descriptive error" requirement
  and makes errors trajectory-dependent.

---

## R4. Access-point-to-link attachment

**Decision**: An access point's `node_id` field directly matches the `source` or `target` node ID of
road links in the network — the shapefiles share a common node ID scheme by design. Link attachment is
determined by **direct ID lookup**, not spatial matching:

1. During `AccessPointLayer` construction, look up the access point's `node_id` in the network's
   `node_id → NodeIndex` map (already built by `Network` during shapefile load).
2. If the node ID exists in the network, the access point is attached to that node, and connects to
   all road edges incident to it (queryable via `Network::get_edges_at_node()` or by scanning the
   edge list for `source == node_id || target == node_id`).
3. If the node ID does **not** exist in the network, the access point is **polygon-shared only** (no
   link attachment). Valid per spec — must be attached to ≥ 2 polygons.

**Rationale**: The user clarified that access points and network links share a node ID scheme. Direct
ID-based attachment is simpler, faster (O(1) lookup vs. O(log N) R-tree query), and avoids any
EPSILON-tolerance ambiguity. No geometric coincidence check is needed for link attachment.

**Implication**: `boundary_epsilon` is now used only for the polygon-boundary validation in
`AccessPointLayer` (FR-005 condition 1) — it is no longer applied to link attachment.

**Alternatives considered**:
- *Spatial coincidence check between access point geometry and network node geometry*: Rejected — the
  user explicitly stated the IDs match, making spatial matching redundant and slower.
- *Make link attachment an explicit field in the shapefile*: Rejected — the user described two
  attachment modes (link OR multiple polygons) as inherent to the data; the node ID match is
  sufficient signal.

---

## R5. Polygon transition cost computation strategy

**Decision**: Split the four cost types into **precomputed** vs. **observation-dependent**:

| Cost type | When computed | Why |
|-----------|--------------|-----|
| **Entry** (outside → inside) | At routing time, during arc relaxation | Depends on the matched-point-inside-polygon (observation-dependent). |
| **Egress** (inside → outside) | At routing time, during arc relaxation | Symmetric to entry. |
| **Within-polygon** (inside → inside, same polygon) | At HMM transition update time | Uses `path_distance = eu_dist`, a per-GPS-pair value mirroring the same-link override at `src/mm/weightmatch/weightmatch_algorithm.cpp:323-324`. Yields the most favorable transition probability via `calc_tp(eu_dist, eu_dist)`. |
| **Through-routing** (outside → polygon → outside) | **Precomputed** at graph build time | See R11. All inputs (polygon weight, access-point geometries, `THROUGH_PENALTY_FACTOR`) are known at startup. |

`THROUGH_PENALTY_FACTOR` is loaded once during config parsing.

**Rationale**: Entry, egress, and within-polygon costs all depend on per-trajectory data (the GPS
matched point, or the GPS-to-GPS Euclidean distance), so they MUST be computed during the HMM/Dijkstra
inner loop. Through-routing cost depends only on inputs known at startup, so it should be precomputed
per Constitution Principle I (no per-query computation in hot paths).

**Alternatives considered**:
- *Precompute everything*: Rejected — entry/egress depend on observation-dependent data.
- *Compute everything on-the-fly*: Rejected — through-routing is the most common polygon-relaxation
  case during Dijkstra (any polygon passed through without a GPS match contributes one through arc
  relaxation per neighbor). Precomputing it saves a `sqrt` + two multiplications per relaxation,
  reducing inner-loop arithmetic and improving cache locality.
- *Apply through-penalty at backtracking*: Rejected — would not affect Viterbi probabilities and would
  let the algorithm prefer through-routing when it shouldn't.

---

## R6. Output extension

**Decision**: A new writer `FMM::IO::PolyMMWriter` produces a CSV file with the same columns as
`CSVMatchResultWriter` **plus** four new columns:
- `polygon_ids` — semicolon-separated list of polygon IDs in the path order.
- `entry_aps` — semicolon-separated list of entry access point node numbers (empty token where absent).
- `egress_aps` — semicolon-separated list of egress access point node numbers (empty token where
  absent).
- `is_through` — semicolon-separated list of `0`/`1` flags per polygon segment.

The existing `CSVMatchResultWriter` is **not** modified. PolyMatch uses `PolyMMWriter` exclusively.

**Rationale**: User constraint — "Existing functionality of other matching types should not be changed."
A new writer satisfies FR-010 (output identifies link vs polygon segments, records access points, marks
through-routing) without touching shared I/O code.

**Alternatives considered**:
- *Extend `CSVMatchResultWriter` with optional polygon columns*: Rejected — adds risk of regression to
  existing matchers' output.
- *Add a second polygon-only CSV file alongside the existing CSV*: Rejected — splits the result across
  two files, harder for downstream consumers; semicolon-separated columns in one CSV is the existing
  project convention (e.g., `cpath`).

---

## R7. Mid-polygon initialization and termination

**Decision**:
- During candidate generation for the **first** GPS point, polygon candidates have a transition cost of
  zero (no entry access point applied) — i.e., the HMM initial distribution treats polygon candidates
  the same as link candidates. The matched first-point's polygon segment records `entry_ap = empty`.
- Symmetric for the **last** GPS point — polygon candidates pay no egress cost, and `egress_ap = empty`.

This is implemented as a flag passed into the transition update for layer 0 and layer N-1.

**Rationale**: Matches spec FR-007 ("MUST additionally support trajectories that begin or end
mid-polygon") and FR-010 ("absent access point MUST be represented as empty").

**Alternatives considered**:
- *Require a synthetic entry/egress access point for every polygon-initial trip*: Rejected — would
  contradict the spec edge case and require fabricating data.

---

## R8. Result type extension

**Decision**: A new struct `FMM::MM::PolyMatchResult` wraps a standard `MatchResult` and adds:
```cpp
struct PolygonSegment {
  PolygonID polygon_id;
  NodeID entry_ap;      // kNoAccessPoint if absent
  NodeID egress_ap;     // kNoAccessPoint if absent
  bool is_through;
};

struct PolyMatchResult {
  MatchResult base;                          // unchanged from existing pipeline
  std::vector<PolygonSegment> polygon_segments;  // empty if pure link-only result
};
```

The base `MatchResult` type is **not** modified.

**Rationale**: Additive composition keeps existing `MatchResult` consumers (FMM, STMATCH, WEIGHTMATCH,
H3MM, the existing `CSVMatchResultWriter`) untouched.

**Alternatives considered**:
- *Add a `polygon_segments` field directly to `MatchResult`*: Rejected — touches the shared type used
  by every algorithm, violating the no-regression constraint.
- *Use std::variant<MatchResult, PolyMatchResult>*: Rejected — overkill; composition is cleaner.

---

## R9. Boundary point handling

**Decision**: For point-in-polygon classification (used by FR-006), use
`boost::geometry::covered_by(point, polygon)` which returns true for points **inside or on the
boundary**. This satisfies the spec's "boundary points are treated as inside (distance zero)" edge
case without extra code.

**Rationale**: Direct support in Boost.Geometry; no custom code needed.

---

## R13. OpenMP-parallel matching (v1 requirement)

**Decision**: Parallelize at the **outer trajectory loop**, mirroring the existing WeightMatch pattern.
Each thread processes one trajectory at a time end-to-end (candidate generation → HMM → Viterbi →
result assembly).

**Thread model**:

| Resource | Sharing model |
|----------|---------------|
| `Network`, `LinkGraph` | Read-only, shared across all threads (already true in WeightMatch). |
| `PolygonLayer`, `AccessPointLayer` | Read-only, shared. Built once at startup, never mutated. |
| `PolyLinkGraph` (incl. through-routing precomputed table) | Read-only, shared. Built once at startup. |
| `DijkstraState` | **Per-thread** instance. Sized to `|E|+|P|` once, reused across trajectories on the same thread. |
| `IndexedMinHeap` | **Per-thread** instance. |
| `TransitionGraph` | **Per-thread** instance. Reused across trajectories. |
| `PolyMMWriter` (CSV output) | Shared, with an internal mutex guarding `write_result()` so completed rows are appended atomically. Mirrors the existing `CSVMatchResultWriter` thread-safety pattern. |

The outer loop is annotated with `#pragma omp parallel for schedule(dynamic)` matching the existing
WeightMatch pattern in `weightmatch_app.cpp`. Per-thread state is created inside the parallel region as
firstprivate or as thread-local variables, exactly as WeightMatch does.

**Determinism**: Result rows may be written in non-deterministic order (because trajectories complete
in non-deterministic order), but the **content** of each row is deterministic — matching is a pure
function of trajectory + shared graph data. SC-009 verification sorts output rows by trajectory ID
before diffing.

**Pitfalls avoided**:
- No `std::shared_ptr` reference counting in hot paths (would cause cache-line contention across
  threads).
- No lazy initialization in shared structures (would race).
- Logging from inside the parallel region uses spdlog's thread-safe sinks (already the project default).

**Rationale**: Production inputs are large batches of trajectories. The natural parallelism is
embarrassingly parallel at the trajectory level, and the existing WeightMatch already uses this
pattern. The through-routing precomputation (R11) and per-polygon static lookups (R4) are read-only
and safe to share. The only synchronization required is at the writer.

**Alternatives considered**:
- *Inner-loop parallelism (Dijkstra)*: Rejected — Dijkstra is hard to parallelize correctly and the
  per-query workload is small. Outer-loop parallelism gives much better scaling.
- *Per-thread output files merged at end*: Rejected — adds complexity and complicates streaming
  consumers. A mutex-protected single writer is sufficient and matches existing project convention.

---

## R12. Distance-inside-polygon computation

**Decision**: For each polygon segment in the matched output, compute the distance inside the polygon
during result assembly (the same phase that builds `PolygonSegment` records). Computation is per-segment
and uses values already available at that point:

```text
For each polygon segment:
  d = 0
  if entry_ap is present and there is ≥ 1 inside GPS point:
    d += distance(entry_AP coords, first inside GPS coords)
  d += Σ eu_dist(consecutive inside GPS coords)           # already cached as eu_dist used by HMM
  if egress_ap is present and there is ≥ 1 inside GPS point:
    d += distance(last inside GPS coords, egress_AP coords)
  if is_through (no inside GPS points):
    d = distance(entry_AP coords, egress_AP coords)
```

**Reuse**: The `Σ eu_dist(consecutive inside GPS coords)` term overlaps with values already computed
during HMM layer construction (the same `eu_dist` used in WeightMatch's same-link override at
`weightmatch_algorithm.cpp:323-324`). Cache these per consecutive-pair so the distance computation is
O(number of polygon segments) extra work, not O(GPS points × candidates).

**Cost**: Negligible — O(P) per trajectory where P is the number of polygon segments in the matched
path (typically a small constant).

**Rationale**: Direct, observable, and exactly testable against hand-computed fixtures (SC-008).
Reusing `eu_dist` values avoids recomputing distances already used by the HMM.

**Alternatives considered**:
- *Use geometric "intersect matched LineString with polygon, sum length"*: Rejected — adds a heavy
  Boost.Geometry operation per segment, and the matched LineString itself is built later in the
  pipeline. The Euclidean-path formulation is consistent with the rest of the cost model.
- *Use polygon weight × distance as the output*: Rejected — the user asked for distance, not cost. Weight
  is already an input; users can multiply on their own if needed.

---

## R11. Through-routing cost precomputation

**Decision**: At `PolyLinkGraph` construction time, precompute a per-polygon table of
`weight × distance(AP_i, AP_j)` for every ordered pair of access points on that polygon. Store as
`std::vector<double>` indexed by `(local_ap_i * n_aps + local_ap_j)` where `local_ap_i` indexes into
the polygon's local access-point list.

Final through-routing cost = `precomputed[i,j] * THROUGH_PENALTY_FACTOR` (one multiplication at lookup
time). Storing the factor-independent quantity (not `weight × dist × factor`) preserves the ability to
vary `THROUGH_PENALTY_FACTOR` for testing without rebuilding the graph.

**Lookup pattern in Dijkstra**: When the search relaxes an outgoing arc from polygon vertex P, the
relaxation routine reads `parent[P]` to identify the entry access point (which by invariant must be
one of P's access points or a link node attached to one of them), computes the egress access point
from the arc's target, and indexes into P's precomputed table.

**Space cost**: `Σ_p (n_p × n_p)` doubles where `n_p` is the number of access points on polygon p.
For typical inputs (10⁴ polygons × ~5 APs each), this is ~200 KB — negligible.

**Time cost**: `O(Σ_p n_p²)` distance computations once at load. For the same inputs, ~250,000
distances — well under a second.

**Rationale**: All inputs are known at startup; computation is deterministic and observation-
independent. Precomputation saves `sqrt + 2 muls` per through-routing arc relaxation, aligning with
Constitution Principle I's no-per-query-computation rule. Also makes through-routing costs
unit-testable in isolation against a known fixture (supporting spec SC-006).

**Alternatives considered**:
- *Precompute the full cost including `THROUGH_PENALTY_FACTOR`*: Rejected — couples the table to the
  config value, requiring rebuild on factor change. The one extra multiplication is negligible.
- *Expand polygon vertices into n entered-from-AP_i vertices (Held-Karp-style state expansion) so arc
  costs become fully static*: Rejected — adds O(n) vertices per polygon and complicates the goal-set
  representation in `shortest_edge_to_edges`. The lookup-via-parent design is simpler and equally fast.
- *Store the table on the `Polygon` struct instead of `PolyLinkGraph`*: Rejected — keeps `Polygon` a
  pure data carrier in `FMM::NETWORK`; the routing-specific precomputation belongs in `FMM::ROUTING`.

---

## R10. Empty / no-polygon-layer fallback

**Decision**: At app-config validation, if `polygon_layer` is empty (not specified), `PolyLinkGraph`
construction skips the polygon-vertex range, leaving a pure `LinkGraph`-equivalent graph. The matcher
behaves identically to WeightMatch in this mode. SC-002 (zero regression vs WeightMatch) is verified by
unit test.

If the polygon layer is specified but contains zero valid polygons (e.g., all skipped by FR-013), the
same fallback applies, plus a warning is emitted.

**Rationale**: Satisfies FR-012 "fall back to link-only matching when no polygon layer is provided or
when the polygon layer contains zero valid polygons" with minimal code complexity.

---

## Summary of Resolved Decisions

| # | Topic | Decision |
|---|-------|----------|
| R1 | Routing graph design | New `PolyLinkGraph` with generalized vertex ID space |
| R2 | Polygon loading | GDAL OGR + Boost.Geometry, mirror existing pattern |
| R3 | Access point loading | New `AccessPointLayer` with three-condition validation pass |
| R4 | Link attachment | Geometric coincidence with road network node (EPSILON tolerance) |
| R5 | Cost computation | Entry/egress/within at routing time; through-routing precomputed (see R11) |
| R6 | Output extension | New `PolyMMWriter` with 4 added columns; existing writer untouched |
| R7 | Mid-polygon HMM | Initial/terminal layers skip entry/egress cost; empty AP in output |
| R8 | Result type | `PolyMatchResult` composes `MatchResult` + polygon segments |
| R9 | Boundary points | `boost::geometry::covered_by` handles boundary = inside |
| R10 | Fallback | `PolyLinkGraph` collapses to `LinkGraph`-equivalent when no polygons |
| R11 | Through-routing precomputation | Per-polygon table of `weight × dist(AP_i, AP_j)`; multiply by factor at lookup |
| R12 | Distance-inside-polygon computation | At result-assembly time, sum entry-to-first + between-inside + last-to-egress (or entry-to-egress for through-routing) |
| R13 | OpenMP-parallel matching (v1) | Outer trajectory loop parallelized; shared graph state read-only; per-thread Dijkstra/heap/transition-graph; mutex-guarded writer |

All Phase 0 unknowns resolved. No NEEDS CLARIFICATION markers remain.
