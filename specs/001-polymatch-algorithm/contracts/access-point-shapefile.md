# Contract: Access Point Shapefile Input

The access point layer is provided as an ESRI shapefile. **One feature per polygon-access-point
relationship** — a single physical access point shared between multiple polygons appears as multiple
features with the same `node_id`.

## Geometry

- Type: `wkbPoint`.
- Coordinate reference system: MUST match the road network and polygon layer CRS.
- Each point geometry MUST lie on the boundary of its declared polygon (within `boundary_epsilon`,
  default `1e-6` map units).

## Required attribute fields

| Field | Type | Default name | Description |
|-------|------|--------------|-------------|
| Node ID | Integer or Int64 | `node_id` | Access point node number. Multiple features may share this value to indicate a shared access point. |
| Polygon ID | Integer or Int64 | `polygon_id` | The polygon this feature associates the access point with. Must reference an existing polygon. |

## Field name overrides

CLI / XML config may override the default field names — see `polymatch-cli.md`
(`--ap_node_id_name`, `--ap_polygon_id_name`).

## Validation rules at load (FR-005)

Load **halts with a descriptive error** if any of the following conditions are found:

1. **Boundary mismatch**: an access point feature's geometry is not on the boundary of its declared
   polygon (distance > `boundary_epsilon`).
2. **Orphaned polygon reference**: an access point feature's `polygon_id` does not exist in the polygon
   shapefile.
3. **Contradictory shared geometry**: two features sharing the same `node_id` have differing geometries
   (distance > `boundary_epsilon`).

## Warnings (not errors)

- A polygon with no access point features referencing it is **excluded** from candidate generation; a
  warning is emitted (one per excluded polygon, per FR-014).

## Examples

### Single access point on one polygon

| node_id | polygon_id | geometry |
|---------|------------|----------|
| 1001 | 42 | POINT (500.0 200.0) |

### Shared access point between two polygons

| node_id | polygon_id | geometry |
|---------|------------|----------|
| 2002 | 42 | POINT (510.0 200.0) |
| 2002 | 43 | POINT (510.0 200.0) |  ← identical geometry; required by validation rule 3 |

### Access point on three polygons

| node_id | polygon_id | geometry |
|---------|------------|----------|
| 3003 | 50 | POINT (700.0 300.0) |
| 3003 | 51 | POINT (700.0 300.0) |
| 3003 | 52 | POINT (700.0 300.0) |

## Link attachment (by node ID match, not spatial matching)

The `node_id` field shares the road network's node ID scheme. The loader determines link attachment by
**direct ID lookup**:

- If `node_id` matches a `source` or `target` node ID of any road link in the network, the access point
  is link-attached to that node — it connects to all road edges incident to that node.
- If `node_id` does not match any network node ID, the access point is polygon-shared only (no link
  attachment); it must reference ≥ 2 polygons (validated implicitly by FR-007 routing rules).

No spatial coincidence check is performed for link attachment; access point geometry is used only for
the polygon-boundary validation (FR-005 condition 1).
