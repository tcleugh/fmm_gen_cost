# Contract: PolyMatch Output CSV

Output produced by `FMM::IO::PolyMMWriter`. Mirrors the existing `CSVMatchResultWriter` output with
four additional columns for polygon traversal information.

## Format

Semicolon-separated CSV (`;` delimiter — matching existing FMM convention). Sub-lists within a field
use `,` separators.

## Columns

| Column | Source | Type | Description |
|--------|--------|------|-------------|
| `id` | trajectory ID | int64 | Carried from input GPS data. |
| `opath` | per-point edge ID | list<int64> | Existing — one matched edge ID per GPS point. May contain polygon IDs as negative numbers when the matched element is a polygon (see polygon ID encoding below). |
| `error` | per-point distance | list<double> | Existing — emission distance per GPS point. |
| `offset` | per-point offset | list<double> | Existing — projection offset per GPS point. |
| `spdist` | inter-point distance | list<double> | Existing — straight-line distance between consecutive points. |
| `cpath` | complete path | list<int64> | Existing — sequence of traversed elements. Polygons appear as negative IDs (see encoding). |
| `tpath` | per-point cpath index | list<int> | Existing — index into cpath for each GPS point. |
| `mgeom` | matched geometry | WKT | Existing — matched path geometry. |
| `polygon_ids` | **NEW** | list<int64> | Polygon IDs in cpath order (one entry per polygon segment). |
| `entry_aps` | **NEW** | list<int64\|empty> | Entry access point node numbers per polygon segment. Empty token (`-`) where absent (mid-polygon start). |
| `egress_aps` | **NEW** | list<int64\|empty> | Egress access point node numbers per polygon segment. Empty token (`-`) where absent (mid-polygon end). |
| `is_through` | **NEW** | list<int> | Per-polygon-segment flag: `1` if through-routing (no GPS inside), `0` otherwise. |
| `polygon_distances` | **NEW** | list<double> | Per-polygon-segment Euclidean distance spent inside the polygon (FR-016). Includes entry/egress access-point contributions when present; equals entry-to-egress straight-line distance for through-routing segments. |

## Polygon ID encoding in `opath` / `cpath`

To distinguish polygons from road edges in the existing `opath`/`cpath` integer columns without
breaking existing parsers, polymatch encodes polygon IDs as **negated values**:

- Road edge `EdgeID = 42` → appears as `42` in cpath.
- Polygon `PolygonID = 7` → appears as `-7` in cpath.
- `PolygonID = 0` is reserved and rejected at load (to avoid ambiguity with edge ID 0).

The `polygon_ids` column independently lists polygon IDs in their **positive** form, in the order they
appear in cpath. Consumers can use either representation.

## Empty access point representation

Where an access point is absent (mid-polygon start or end, per FR-010), the field uses the literal
token `-` (single hyphen). Example:

```text
polygon_ids: 7,42
entry_aps:   -,1001
egress_aps:  1001,-
is_through:  0,0
```

This represents a trip that began inside polygon 7 (no entry AP), exited via AP 1001 into polygon 42
(also AP 1001 as entry of polygon 42 — i.e., the shared access point), and ended inside polygon 42 (no
egress AP).

## Output field selection

The `--output_fields` flag (existing `ResultConfig` feature) accepts the new column names as additional
options. Default output includes all four new columns when polygons are configured; the columns are
omitted entirely when running in link-only fallback mode (no polygon layer provided), preserving
binary-identical CSV output for SC-002 regression checks.

## Example row

```text
id;opath;cpath;polygon_ids;entry_aps;egress_aps;is_through;polygon_distances;mgeom
1001;1,2,-7,-7,3,4;1,2,-7,3,4;7;1001;2002;0;47.83;LINESTRING(...)
```

This represents trajectory 1001 matching to edges 1,2, then polygon 7 (twice — two GPS points inside),
then edges 3,4. The polygon was entered via AP 1001 and exited via AP 2002, is not a through-route, and
the matched path spent 47.83 (map units) of Euclidean distance inside polygon 7 (entry-AP-to-first-GPS
+ between-GPS + last-GPS-to-egress-AP).
