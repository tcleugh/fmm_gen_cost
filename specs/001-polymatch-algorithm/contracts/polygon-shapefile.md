# Contract: Polygon Shapefile Input

The polygon layer is provided as an ESRI shapefile (`.shp` + `.dbf` + `.shx` + optional `.prj`).

## Geometry

- Type: `wkbPolygon` or `wkbMultiPolygon`.
- Coordinate reference system: MUST match the road network's CRS. No on-the-fly reprojection is
  performed.

## Required attribute fields

| Field | Type | Default name | Description |
|-------|------|--------------|-------------|
| Polygon ID | Integer or Int64 | `id` | Unique identifier per polygon. May be discontinuous or negative. |

## Optional attribute fields

| Field | Type | Default name | Default value | Description |
|-------|------|--------------|---------------|-------------|
| Weight / cost | Real | `cost` | `1.0` | Generalized cost weight (FR-008). |

## Field name overrides

CLI / XML config may override the default field names — see `polymatch-cli.md`
(`--polygon_id_name`, `--polygon_cost_name`).

## Validation rules at load

### Hard-reject rules (FR-018) — load halts with descriptive error

Load **halts** with a descriptive error if any of the following conditions are found:

- A polygon feature has `id == 0`. Polygon ID 0 is reserved by the output format (polygons are encoded
  as negated IDs in `opath`/`cpath`; ID 0 cannot be disambiguated from edge ID 0) and cannot appear in
  input data.
- Two or more polygon features share the same `id`. The error message MUST list the offending ID.

### Soft-skip rules (FR-013) — warning, matching continues

A polygon feature is **skipped with a warning** (matching continues with the remaining valid polygons)
if:
- Its geometry is empty, self-intersecting, or zero-area.
- Its weight value is negative or NaN.

A polygon feature is **accepted** otherwise.

## Example layout

```text
shp:  polygons.shp        (geometry)
dbf:  polygons.dbf        (attributes)
       fields: id (int64), cost (real)
shx:  polygons.shx        (index)
prj:  polygons.prj        (CRS — same as network)
```
