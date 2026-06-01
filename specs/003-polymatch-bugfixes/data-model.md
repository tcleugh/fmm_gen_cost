# Phase 1 Data Model: Polymatch Bug Fixes

This feature adds no new entities. It clarifies the semantics of two existing fields whose meaning was muddy in feature 001 and surfaced as bugs in feature 002.

---

## Clarified semantics: `PolygonSegment::is_through`

**Before this feature**: "True if no GPS observations fell inside this polygon segment (through-routing)" — with `inside` interpreted strictly as `polygons_containing` / `covered_by` semantics.

**After this feature**: "True if and only if no polygon Viterbi candidate for this segment's polygon was matched at any GPS layer along the trajectory." Equivalently: the polygon is in the output cpath *purely* because it was a routing detour discovered by the link↔link Dijkstra over `PolyLinkGraph`, not because the HMM picked the polygon as the best per-layer candidate for any GPS point.

**Concrete contract** (post-fix):

- A trip whose Viterbi path has *any* PolyCandidate-of-Polygon-P at any layer → every emitted PolygonSegment for P has `is_through = false`.
- A trip whose Viterbi path has NO PolyCandidate of Polygon P but whose link→link Dijkstra crosses Polygon P → the emitted PolygonSegment for P has `is_through = true` AND both `entry_ap` and `egress_ap` populated.

**Implication for the invariant** (`check_is_through_has_aps` in the harness):

```
is_through == true   ⇒   entry_ap != kNoAccessPoint  AND  egress_ap != kNoAccessPoint
```

Strict, applies to ALL polygon segments. No first/last-segment exemption.

---

## Clarified semantics: `PolyCandidate::inside`

**Before this feature**: "True if GPS observation is inside / on boundary of polygon." Was the sole gate on `record_inside()` (which drove `has_inside_obs` which drove `is_through`).

**After this feature**: Still "true if `polygons_containing(gps)` returned the polygon"; still informs `distance_inside` math (the `if (has_inside_obs && entry_ap != kNoAccessPoint) d += point_distance(entry_point, first_inside)` branch). But NO LONGER the gate on `record_inside` — every polygon Viterbi candidate match calls `record_inside` regardless of `inside`.

Internal field; no change to its declared type or location.

---

## Clarified contract: `build_hybrid_path` polygon→link emission

(Pending Phase 0 / Phase 1 diagnostic on traces 1313 + 1314)

The contract this feature MUST restore is:

> For every consecutive `(polygon, edge)` or `(edge, polygon)` pair in the emitted `cpath`, the edge has at least one endpoint (`source` or `target` node ID) that is the `node_id` of an `AccessPoint` whose `polygons` field contains the polygon's index.

Equivalently: the matcher must NEVER emit a polygon adjacent to an edge that has no AP-incidence with that polygon. Today the matcher mostly does this correctly (198/200 real-network traces pass) but fails on traces 1313 + 1314. The Phase 1 implementation work identifies the buggy code path and corrects it.

---

## No new entities

No new structs, no new enums, no new entity files. The feature is purely a behavioral correction to existing data flow.
