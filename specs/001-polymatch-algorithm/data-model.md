# Phase 1 Data Model: PolyMatch

This document defines the entities introduced by PolyMatch, their fields and relationships, and the
invariants that hold across them. Existing entities (`Network`, `Edge`, `LinkGraph`, `Trajectory`,
`Candidate`, `TransitionGraph`, `MatchResult`, `C_Path`) are referenced but **not modified**.

---

## Entity: `Polygon`

**Namespace**: `FMM::NETWORK`
**Header**: `src/network/polygon_layer.hpp`

| Field | Type | Description |
|-------|------|-------------|
| `index` | `PolygonIndex` (uint32_t) | Contiguous internal index `[0, |P|)`. |
| `id` | `PolygonID` (int64_t) | Original ID from shapefile; may be discontinuous. |
| `geom` | `boost::geometry::model::polygon<point_2d>` | Boundary geometry. |
| `weight` | `double` | Generalized cost weight (default `1.0` if absent in shapefile). |
| `bbox` | `boost::geometry::model::box<point_2d>` | Cached bounding box for R-tree. |

**Invariants**:
- `index` is contiguous starting at zero.
- `geom` is a valid simple polygon (no self-intersection, non-zero area). Validated at load — invalid
  polygons are dropped per FR-013.
- `weight >= 0`. (Negative weights are rejected at load.)

**Loaded from**: ESRI shapefile via GDAL OGR; field names configurable through `PolygonConfig`.

---

## Entity: `AccessPointFeature`

**Namespace**: `FMM::NETWORK`
**Header**: `src/network/access_point_layer.hpp`

Represents one row in the access point shapefile (a single polygon-AP relationship). Multiple features
with the same `node_id` represent a single shared access point.

| Field | Type | Description |
|-------|------|-------------|
| `node_id` | `NodeID` (int64_t) | Access point node number (shared key across rows). |
| `polygon_id` | `PolygonID` (int64_t) | The polygon this feature associates the access point with. |
| `point` | `point_2d` (Boost.Geometry) | Geometric location on the polygon boundary. |

**Invariants** (enforced at load — load halts on any failure per FR-005):
- `point` is within EPSILON (`1e-6` map units, configurable) of the boundary of the polygon identified
  by `polygon_id`.
- `polygon_id` corresponds to an existing polygon in the loaded `PolygonLayer`.
- All features sharing the same `node_id` have identical `point` geometries (within EPSILON).

---

## Entity: `AccessPoint` (resolved)

**Namespace**: `FMM::NETWORK`
**Header**: `src/network/access_point_layer.hpp`

The resolved, deduplicated access point — one per unique `node_id`. Built from `AccessPointFeature`s
after validation.

| Field | Type | Description |
|-------|------|-------------|
| `index` | `AccessPointIndex` (uint32_t) | Contiguous internal index. |
| `node_id` | `NodeID` | Node number from shapefile; identical to the road network's node ID scheme. |
| `point` | `point_2d` | Geometric location (on polygon boundary). |
| `polygons` | `std::vector<PolygonIndex>` | All polygons this access point connects to (size ≥ 1). |
| `attached_node` | `optional<NodeIndex>` | Road network `NodeIndex` for `node_id`, if `node_id` exists in the network's node map. Determined by direct ID lookup — not spatial matching. |
| `attached_edges` | `std::vector<EdgeIndex>` | Cached list of road edges incident to `attached_node` (empty if `attached_node` is absent). |

**Relationships**:
- An `AccessPoint` references ≥ 1 `Polygon` via the `polygons` field.
- An `AccessPoint` either has `attached_node.has_value()` (link-attached — `node_id` matched a network
  node) or `polygons.size() >= 2` (polygon-shared). The spec guarantees at least one of these
  (FR-004 / user clarification).

---

## Entity: `PolygonLayer`

**Namespace**: `FMM::NETWORK`
**Header**: `src/network/polygon_layer.hpp`

Container for all loaded polygons + spatial index.

| Field | Type | Description |
|-------|------|-------------|
| `polygons` | `std::vector<Polygon>` | Indexed by `PolygonIndex`. |
| `id_to_index` | `std::unordered_map<PolygonID, PolygonIndex>` | Lookup by original ID. |
| `rtree` | `boost::geometry::index::rtree<...>` | Polygons indexed by bbox for spatial search. |

**Operations**:
- `polygons_containing(point) -> std::vector<PolygonIndex>` — candidate polygons for a GPS point.
- `polygons_within_radius(point, radius) -> std::vector<PolygonIndex>` — candidates within radius.
- `min_boundary_distance(polygon_idx, point) -> double` — for emission probability (FR-006).

---

## Entity: `AccessPointLayer`

**Namespace**: `FMM::NETWORK`
**Header**: `src/network/access_point_layer.hpp`

Container for resolved access points + lookups.

| Field | Type | Description |
|-------|------|-------------|
| `access_points` | `std::vector<AccessPoint>` | Indexed by `AccessPointIndex`. |
| `node_id_to_index` | `std::unordered_map<NodeID, AccessPointIndex>` | Lookup by node number. |
| `polygon_to_aps` | `std::unordered_map<PolygonIndex, std::vector<AccessPointIndex>>` | Reverse map. |

---

## Entity: `PolyLinkGraph`

**Namespace**: `FMM::ROUTING`
**Header**: `src/network/poly_link_graph.hpp`

The polygon-aware routing graph. Vertices represent both road edges and polygons in a single ID space.

| Field | Type | Description |
|-------|------|-------------|
| `n_edges` | `size_t` | `|E|` — number of road edges. |
| `n_polygons` | `size_t` | `|P|` — number of polygons. |
| `adjacency` | `std::vector<std::vector<Arc>>` | Outgoing arcs per vertex; size `|E| + |P|`. |
| `vertex_kind(v)` | inline | Returns `Edge` if `v < n_edges`, else `Polygon`. |
| `through_cost_tables` | `std::vector<std::vector<double>>` | Per polygon `p`, a flattened `n_p × n_p` matrix of factor-independent through-routing weights `polygon.weight × dist(AP_i, AP_j)`. Indexed by `local_ap_i * n_p + local_ap_j` where `local_ap_*` indexes into the polygon's local AP list. Final cost = `through_cost_tables[p][i*n+j] * THROUGH_PENALTY_FACTOR` at lookup. |
| `polygon_local_ap_index` | `std::vector<std::unordered_map<AccessPointIndex, uint16_t>>` | Maps a global `AccessPointIndex` to its local index in polygon `p`'s AP list, for O(1) table lookups during Dijkstra relaxation. |

**Arc types**:
- **link→link**: copied from existing `LinkGraph` adjacency.
- **link→polygon**: from any edge incident to access point AP, to polygon vertex for each polygon in
  `AP.polygons`. Arc weight = entry-cost placeholder; final cost computed at relaxation time using the
  matched-point-inside-polygon.
- **polygon→link**: symmetric.
- **polygon→polygon**: when AP has ≥ 2 polygons, every pair (P_i, P_j) in `AP.polygons` gets bidirectional
  arcs. Arc weight = through-routing cost placeholder.

**Big-O**:
- Construction: O(|E| + Σ |AP.polygons| + Σ |AP.attached_links| + Σ_p n_p²) where the last term is the
  through-routing precomputation (R11).
- Dijkstra query: O((|E|+|P|) log (|E|+|P|)) using reused `IndexedMinHeap`. Each through-routing arc
  relaxation is O(1) (table lookup + one multiplication by `THROUGH_PENALTY_FACTOR`).

**Invariant**: Vertex IDs `[0, n_edges)` exactly correspond to `EdgeIndex` values from the underlying
`Network`. Polygon vertex `n_edges + p_idx` corresponds to `PolygonIndex` `p_idx`.

---

## Entity: `PolygonSegment`

**Namespace**: `FMM::MM`
**Header**: `src/mm/polymatch/poly_match_result.hpp`

One segment in a hybrid matched path that corresponds to a polygon traversal.

| Field | Type | Description |
|-------|------|-------------|
| `polygon_id` | `PolygonID` | Polygon traversed. |
| `entry_ap` | `NodeID` | Entry access point node number, or `kNoAccessPoint` (sentinel `-1`) if absent (mid-polygon start). |
| `egress_ap` | `NodeID` | Egress access point node number, or `kNoAccessPoint` if absent (mid-polygon end). |
| `is_through` | `bool` | True if no GPS observations fell inside this polygon segment (through-routing). |
| `distance_inside` | `double` | Euclidean path length spent inside the polygon (FR-016). For traversal: `dist(entry_AP, first inside GPS) + Σ dist(consecutive inside GPS) + dist(last inside GPS, egress_AP)`, dropping terms whose AP is absent. For through-routing: `dist(entry_AP, egress_AP)`. |
| `position_in_cpath` | `size_t` | Index into the hybrid `C_Path` where this segment appears. |

**Invariants**:
- At most one of `entry_ap` / `egress_ap` may be `kNoAccessPoint` for the **first** and **last**
  polygon segments (mid-polygon start/end). Interior polygon segments must have both AP fields set.
- If `is_through == true`, both `entry_ap` and `egress_ap` MUST be set (not `kNoAccessPoint`).

---

## Entity: `PolyMatchResult`

**Namespace**: `FMM::MM`
**Header**: `src/mm/polymatch/poly_match_result.hpp`

Result of matching a single trajectory with polymatch.

| Field | Type | Description |
|-------|------|-------------|
| `base` | `MatchResult` | Standard result (id, opath, cpath, indices, mgeom). `cpath` may include both edge IDs (positive) and polygon IDs (negative, distinguished by a separate flag). |
| `polygon_segments` | `std::vector<PolygonSegment>` | One entry per polygon traversal in the hybrid path. Empty for pure link-only matches. |

**Invariant**: For every entry in `polygon_segments`, `position_in_cpath` indexes a valid position in
`base.cpath`. The set of polygon positions in `base.cpath` exactly equals the `position_in_cpath` values
of `polygon_segments`.

---

## Entity: `POLYMATCHConfig`

**Namespace**: `FMM::MM`
**Header**: `src/mm/polymatch/polymatch_algorithm.hpp`

Algorithm-level configuration. Mirrors `WEIGHTMATCHConfig` field-for-field, plus polygon-specific
fields.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `k` | `int` | `8` | Number of candidates per GPS point. |
| `radius` | `double` | `300` | Search radius (map units). |
| `gps_error` | `double` | `50` | GPS measurement error (map units). |
| `backup_k` | `int` | `-1` | Backup candidate count if first search empty. |
| `backup_radius` | `double` | `-1` | Expanded backup search radius. |
| `upper_bound_factor` | `double` | `10.0` | Dijkstra upper-bound multiplier. |
| `allow_truncation` | `bool` | `false` | Allow truncating ends if no candidates. |
| `through_penalty_factor` | `double` | `1.5` | Multiplier for through-routing cost (FR-008, FR-009). |
| `boundary_epsilon` | `double` | `1e-6` | Tolerance for access-point-to-polygon-boundary distance check (FR-005). Not used for link attachment — that is by direct node ID match. |

**Validation**: `validate()` checks `k > 0`, `radius > 0`, `gps_error > 0`, `through_penalty_factor >= 0`.

---

## Entity: `PolygonConfig`

**Namespace**: `FMM::CONFIG`
**Header**: `src/config/polygon_config.hpp`

Specifies how to load the polygon layer.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `file` | `std::string` | (required) | Path to polygon shapefile. |
| `id_name` | `std::string` | `"id"` | Field name for polygon ID. |
| `cost_name` | `std::string` | `"cost"` | Field name for polygon weight; optional. |

---

## Entity: `AccessPointConfig`

**Namespace**: `FMM::CONFIG`
**Header**: `src/config/access_point_config.hpp`

Specifies how to load the access point shapefile.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `file` | `std::string` | (required) | Path to access point shapefile. |
| `node_id_name` | `std::string` | `"node_id"` | Field name for the access point node number. |
| `polygon_id_name` | `std::string` | `"polygon_id"` | Field name for the polygon association. |

---

## Entity: `POLYMATCHAppConfig`

**Namespace**: `FMM::MM`
**Header**: `src/mm/polymatch/polymatch_app_config.hpp`

Top-level CLI config. Mirrors `WEIGHTMATCHAppConfig` plus polygon configs.

| Field | Type | Description |
|-------|------|-------------|
| `network_config` | `CONFIG::NetworkConfig` | (existing) Road network input. |
| `gps_config` | `CONFIG::GPSConfig` | (existing) GPS trajectory input. |
| `result_config` | `CONFIG::ResultConfig` | (existing) Result output config. |
| `polygon_config` | `CONFIG::PolygonConfig` | NEW — polygon layer input. |
| `access_point_config` | `CONFIG::AccessPointConfig` | NEW — access point layer input. |
| `polymatch_config` | `POLYMATCHConfig` | NEW — algorithm config (mirror of WEIGHTMATCHConfig + extras). |
| `use_omp` | `bool` | (existing pattern) |
| `help_specified` | `bool` | (existing pattern) |
| `log_level` | `int` | (existing pattern) |
| `step` | `int` | (existing pattern) |

---

## Relationship Diagram

```text
PolygonLayer  ─────owns─────►  Polygon
      ▲                          ▲
      │                          │ refs by index
      │                          │
AccessPointLayer ──owns──► AccessPointFeature ──validated to─► AccessPoint
      ▲                                                            │
      │                                                            │ refs Polygon by index
      │                                                            ▼
      │                                                       PolygonIndex
      │
      └──────used to build──────► PolyLinkGraph (FMM::ROUTING)
                                       │
                                       ▼
                                  POLYMATCH::match_traj()
                                       │
                                       ▼
                                  PolyMatchResult { base: MatchResult, polygon_segments[] }
                                       │
                                       ▼
                                  PolyMMWriter → CSV
```

---

## State Transitions (HMM)

The HMM state-transition model is unchanged from WeightMatch except for two extensions:

1. **Layer 0 (initial)**: polygon candidates are admitted with the same emission probability rule as
   link candidates (FR-006). No entry cost is charged (Phase-0 R7).
2. **Layer N-1 (terminal)**: symmetric — no egress cost charged for polygon candidates.

Transition costs between consecutive layers use the four-case cost model from FR-008:
- Same polygon, both inside → `path_distance = eu_dist` (override mirroring the existing same-link
  override at `weightmatch_algorithm.cpp:323-324`).
- Different routing elements → Dijkstra on `PolyLinkGraph`, with entry/egress/through costs computed at
  arc-relaxation time.

---

## Invariant Summary (cross-references to spec FRs)

| Invariant | FR |
|-----------|-----|
| Each polygon has unique `PolygonIndex` and unique `PolygonID`. | FR-003 |
| Every `AccessPoint` references ≥ 1 polygon and either attaches to a network node OR ≥ 2 polygons. | FR-004 |
| Every `AccessPointFeature` geometry lies on the boundary of its declared polygon. | FR-005 |
| Polygon emission distance = 0 for inside/boundary points; min boundary distance for outside points. | FR-006 |
| Routing transitions to/from polygons occur via access points only. | FR-007, FR-015 |
| Cost model is applied per arc-relaxation, with through-routing distinguished. | FR-008 |
| `THROUGH_PENALTY_FACTOR` is configurable. | FR-009 |
| Output records polygon ID + entry AP + egress AP per segment; through-routing marked. | FR-010 |
| Empty access point output where mid-polygon start/end. | FR-010 |
| Fallback to link-only matching when no polygons. | FR-012 |
| Invalid polygon geometry skipped with warning. | FR-013 |
| Polygons without access points are excluded. | FR-014 |
| `C_Path` topology preserved through access points. | FR-015 |
| Each polygon segment reports a numeric distance-inside-polygon. | FR-016 |
| Shared graph state read-only post-construction; mutable state per-thread; output serialized; results bit-identical across thread counts. | FR-017 |
| Polygon shapefile load halts with descriptive error if any polygon has `id == 0` or if two polygons share the same `id`. | FR-018 |
| Trajectories with 0 or 1 GPS points are skipped with per-trajectory warning; remaining trajectories continue; no output row produced for skipped trajectories. | FR-019 |
