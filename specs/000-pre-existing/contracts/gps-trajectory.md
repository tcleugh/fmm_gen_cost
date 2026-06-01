# Contract: GPS Trajectory Input

Every matcher consumes a batch of GPS trajectories via `IO::GPSReader`. Three input formats are supported, auto-selected by file extension + the `gps_point` flag.

## Format selection

| Extension | `gps_point` flag | Reader | Internal format ID |
|---|---|---|---|
| `.csv` / `.txt` | unset (default) | `CSVTrajectoryReader` | `1` — one trajectory per row, geometry in a WKT column |
| `.csv` / `.txt` | set | `CSVPointReader` | `2` — one point per row, grouped by id |
| `.gpkg` / `.shp` | (ignored) | `GDALTrajectoryReader` | `0` — one trajectory per feature, geometry from OGR |

Any other extension yields a critical-level log and a -1 sentinel from `GPSConfig::get_gps_format()`.

## Format 1 — CSV trajectory (most common)

Semicolon-delimited (`;`). Header row required; the reader matches header tokens against the configured `id` / `geom` / `timestamp` column names.

```text
id;geom;timestamp
1;LINESTRING(0 0,1 0,2 0);1000,1010,1020
2;LINESTRING(5 5,5 6,5 7);
3;LINESTRING(2 2,3 3);1100,1110
```

Columns:

- **`id`** (required, `int`): trajectory identifier. Default column name `id`.
- **`geom`** (required, WKT): the trajectory polyline. Default column name `geom`. Parsed via `boost::geometry::read_wkt`. Must be `LINESTRING(...)` — point/multi-line variants are not handled.
- **`timestamp`** (optional, comma-separated doubles): per-point timestamps. Default column name `timestamp`. If the column is missing, a single warning is logged and the field is left empty; matching proceeds.

When `timestamp` is present, the per-trajectory count of timestamps **should** equal the polyline's point count. There is no runtime check — a mismatch will surface later when downstream code indexes into the vector.

## Format 2 — CSV point

Same delimiter (`;`). One row per GPS point, sorted by `id` so consecutive rows with the same id form one trajectory. The reader buffers rows until the id changes, then emits a trajectory.

```text
id;x;y;timestamp
1;0;0;1000
1;1;0;1010
1;2;0;1020
2;5;5;
2;5;6;
```

Columns: `id` (required), `x` and `y` (required reals, default column names `x` / `y`), optional `timestamp` (single double per row).

## Format 0 — OGR trajectory (shapefile / GeoPackage)

One trajectory per feature. Geometry must be `wkbLineString`. The configured `id` field maps to the trajectory ID; `timestamp` is parsed via `string2time` if present.

CRS is *not* re-projected — must match the network's CRS.

## Field-name overrides

CLI / XML knobs (see [matcher-cli.md](matcher-cli.md)):

- `--gps_id` / `<gps><id>` → ID column (default `id`)
- `--gps_geom` / `<gps><geom>` → geometry / WKT column (default `geom`, Format 1 + 0)
- `--gps_x` / `--gps_y` / `<gps><x|y>` → coordinate columns (default `x` / `y`, Format 2 only)
- `--gps_timestamp` / `<gps><timestamp>` → timestamp column (default `timestamp`)
- `--gps_point` / `<gps><gps_point>` → switch from Format 1 to Format 2 for CSV inputs

## Validation

`GPSConfig::validate()` checks:

- The file exists.
- The extension maps to a known format.

Per-row schema validation (column presence, geometry parseability) happens at reader-construction time inside `CSVTrajectoryReader::CSVTrajectoryReader(...)`: missing `id` or `geom` columns throw a `std::runtime_error`; a missing `timestamp` column logs a warning.

## Polymatch interaction

Polymatch additionally skips trajectories with fewer than 2 GPS points (FR-019 in `specs/001-polymatch-algorithm/spec.md`). That filter lives in the matcher's `run()` loop, not the GPS reader — the reader emits whatever the file contains.
