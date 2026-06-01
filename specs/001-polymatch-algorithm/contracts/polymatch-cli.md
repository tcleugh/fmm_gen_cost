# Contract: `polymatch` CLI

The `polymatch` executable mirrors `weightmatch`'s CLI shape with two additional layer configurations
(polygon, access point) and one additional algorithm parameter (through-penalty factor).

For the shared flags polymatch inherits from `NetworkConfig`, `GPSConfig`, and `ResultConfig`, see
[../../000-pre-existing/contracts/matcher-cli.md](../../000-pre-existing/contracts/matcher-cli.md).
This file documents only the polymatch-specific additions and overrides.

## Synopsis

```text
polymatch <args...>
polymatch --help
polymatch --config <file.xml>
```

## Arguments

### Network (inherited from `NetworkConfig` — same as weightmatch)

| Flag | Required | Description |
|------|----------|-------------|
| `--network`, `-n` | Y | Path to road network shapefile. |
| `--network_id`, `-i` | N | Field name for edge ID (default `id`). |
| `--source` | N | Field name for source node (default `source`). |
| `--target` | N | Field name for target node (default `target`). |
| `--weight_name` | N | Field name for edge weight (default `weight`). |
| `--cost_name` | N | Field name for edge cost (default `cost`). |
| `--mode` | N | Routing mode flag (default unchanged from existing). |

### GPS (inherited from `GPSConfig` — same as weightmatch)

| Flag | Required | Description |
|------|----------|-------------|
| `--gps`, `-g` | Y | Path to GPS CSV / shapefile. |
| `--gps_id` | N | Field name for trajectory ID. |
| `--gps_geom` | N | Field name for geometry. |
| `--gps_x`, `--gps_y` | N | X/Y column names if using point CSV. |

### Result (inherited from `ResultConfig`)

| Flag | Required | Description |
|------|----------|-------------|
| `--output`, `-o` | Y | Output CSV path. |
| `--output_fields` | N | Comma-separated list of output columns; see `polymatch-output.md`. |

### Polygon layer (NEW — `PolygonConfig`)

| Flag | Required | Description |
|------|----------|-------------|
| `--polygons` | N | Path to polygon shapefile. If omitted, polymatch falls back to weightmatch-equivalent behavior. |
| `--polygon_id_name` | N | Field name for polygon ID (default `id`). |
| `--polygon_cost_name` | N | Field name for polygon weight (default `cost`). |

### Access point layer (NEW — `AccessPointConfig`)

| Flag | Required | Description |
|------|----------|-------------|
| `--access_points` | Conditional | Path to access point shapefile. Required if `--polygons` is provided. |
| `--ap_node_id_name` | N | Field name for AP node number (default `node_id`). |
| `--ap_polygon_id_name` | N | Field name for polygon association (default `polygon_id`). |

### Algorithm parameters (NEW — `POLYMATCHConfig`; superset of `WEIGHTMATCHConfig`)

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-k`, `--candidates` | N | `8` | Candidates per GPS point. |
| `-r`, `--radius` | N | `300` | Candidate search radius. |
| `-e`, `--error` | N | `50` | GPS measurement error. |
| `--backup_k` | N | `-1` | Backup candidate count. |
| `--backup_radius` | N | `-1` | Backup search radius. |
| `-u`, `--ub_factor` | N | `10.0` | Dijkstra upper-bound multiplier. |
| `--allow_truncation` | N | `false` | Allow trimming ends if no candidates. |
| `--through_penalty_factor` | N | `1.5` | Multiplier for through-routing cost (FR-008, FR-009). |
| `--boundary_epsilon` | N | `1e-6` | Tolerance for boundary checks. |

### Runtime

| Flag | Required | Description |
|------|----------|-------------|
| `--log_level` | N | 0-trace … 6-off. |
| `--step` | N | Progress report interval. |
| `--use_omp` | N | Enable OpenMP-parallel trajectory matching (FR-017). When set, trajectories are matched concurrently using all available OpenMP threads; results are bit-identical to single-threaded runs (SC-009). |
| `--config`, `-c` | N | XML config file path (alternative to CLI flags). |

## XML config equivalent

```xml
<config>
  <input>
    <network><file>...</file><id>id</id>...</network>
    <gps>...</gps>
    <polygon>
      <file>polygons.shp</file>
      <id_name>id</id_name>
      <cost_name>cost</cost_name>
    </polygon>
    <access_point>
      <file>access_points.shp</file>
      <node_id_name>node_id</node_id_name>
      <polygon_id_name>polygon_id</polygon_id_name>
    </access_point>
  </input>
  <parameters>
    <k>8</k>
    <r>300</r>
    <gps_error>50</gps_error>
    <through_penalty_factor>2.0</through_penalty_factor>
    <boundary_epsilon>1e-6</boundary_epsilon>
  </parameters>
  <output>...</output>
</config>
```

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success (or help printed). |
| Non-zero | Invalid configuration or fatal load error (e.g., polygon shapefile missing, access-point validation failed). |

## Help output

`polymatch --help` prints the full set of supported flags in the same format as `weightmatch --help`.
