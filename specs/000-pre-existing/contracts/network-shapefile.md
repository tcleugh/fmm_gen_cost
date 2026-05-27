# Contract: Road Network Input

The road network is consumed by every matcher (`fmm`, `stmatch`, `weightmatch`, `h3mm`, `polymatch`, and `ubodt_gen`). It carries the directed graph topology plus a polyline geometry per edge.

## Accepted formats

`Network::Network(filename, turn_ban_file, ...)` reads via GDAL/OGR (`GDALOpenEx(... GDAL_OF_VECTOR ...)`). The extension is gated to one of:

- `.shp` — ESRI Shapefile (with companion `.dbf` + `.shx`; optional `.prj`)
- `.gpkg` — GeoPackage
- `.geojson` — GeoJSON

Other formats GDAL can read in principle are rejected at the extension check; extend the allow-list in `src/network/network.cpp` to add more.

## Geometry

- Type: `wkbLineString` per feature. Mixed / multi-line / point geometries are rejected with a critical-log + thrown exception.
- Coordinate reference system: read from the file when present; defaults to EPSG:4326 with a warning when absent. **All other inputs (GPS, polygon layer, access-point layer for polymatch) MUST share this CRS — there is no on-the-fly reprojection anywhere in the pipeline.**

## Required attribute fields

| Field | Type | Default name | Description |
|---|---|---|---|
| Edge ID | Integer / Int64 | `id` | Unique identifier per edge. Persisted in matched output (`opath`/`cpath`). May be discontinuous or negative. |
| Source node | Integer / Int64 | `source` | Start node ID of the directed edge. |
| Target node | Integer / Int64 | `target` | End node ID. |

Edges are directional: `source -> target`. The graph has no automatic reverse edge; if a road is bidirectional, the input must contain two features.

## Optional attribute fields

| Field | Type | Default name | Default behavior | Description |
|---|---|---|---|---|
| Weight / cost | Real | `NO_WEIGHT` (sentinel — disabled by default) | When the field name is the sentinel `"NO_WEIGHT"`, each edge's weight is `1.0`. | Generalized cost used by `weightmatch` / `polymatch` for Dijkstra. The `Edge::cost` field shares this value. |

## Field-name overrides

Every matcher exposes the same CLI / XML knobs (see [matcher-cli.md](matcher-cli.md)):

- `--network-id` / `<network><id>` → Edge-ID field name (default `id`)
- `--source` / `<network><source>` → Source-node field name (default `source`)
- `--target` / `<network><target>` → Target-node field name (default `target`)
- `--weight` / `<network><weight>` → Weight field name (default `NO_WEIGHT` ≡ disabled)

The sentinel `"NO_WEIGHT"` is the *default* string value, not a placeholder — pass an actual column name to enable weighted edges.

## Internal indexing

- `EdgeID` (long long) — original ID from the shapefile. Carried into `opath`/`cpath` output unchanged.
- `EdgeIndex` (uint32_t) — contiguous 0-based internal index used by all `FMM::ROUTING` graphs. Derived at load.
- Same `NodeID` / `NodeIndex` distinction for nodes.

Edge IDs may collide with polygon IDs in polymatch output if a polygon shares an integer with an edge — polymatch resolves this by negating polygon IDs in `opath`/`cpath` (see `specs/001-polymatch-algorithm/contracts/polymatch-output.md`).

## Validation

`NetworkConfig::validate()` checks:

- File exists.
- Extension is in the accepted set (see top of this file).
- If `turn_ban_file != "NO_TURN_BANS"`, the turn-ban file exists and ends in `.csv`.

The actual OGR load happens lazily inside `Network::Network(...)`; field-not-found and geometry-type-mismatch errors throw `std::runtime_error` at construction.

## See also

- [turn-ban-file.md](turn-ban-file.md) — the optional CSV that accompanies a network.
- [matcher-cli.md](matcher-cli.md) — the CLI flags every matcher exposes for these fields.
