# Contract: Shared Matcher CLI Flags

`fmm`, `stmatch`, `weightmatch`, `h3mm`, `polymatch` each compose the same three sub-configs — `NetworkConfig`, `GPSConfig`, `ResultConfig` — plus their algorithm-specific options + a small set of runtime flags. This file documents the shared surface; algorithm-specific flags live in the matcher's own spec.

Every matcher accepts either:

- Direct CLI flags (`./matcher --network n.shp --gps g.csv ...`), or
- A single XML config file as the only argument (`./matcher config.xml`) — `stmatch`, `fmm`, `ubodt_gen`, and `polymatch` support this; `weightmatch` is CLI-only as of this writing.

## Shared CLI flags

### Network — `CONFIG::NetworkConfig`

| Flag | Required | Default | Description |
|---|---|---|---|
| `--network` | Y | `""` | Path to the network file. See [network-shapefile.md](network-shapefile.md) for accepted formats. |
| `--turn_ban_file` | N | `NO_TURN_BANS` | Path to optional turn-ban CSV. See [turn-ban-file.md](turn-ban-file.md). The literal `NO_TURN_BANS` disables turn-ban loading. |
| `--network_id` | N | `id` | Edge-ID field name. |
| `--source` | N | `source` | Source-node field name. |
| `--target` | N | `target` | Target-node field name. |
| `--weight` | N | `NO_WEIGHT` | Weight field name. The literal `NO_WEIGHT` disables weight loading; per-edge weight defaults to 1.0. |

### GPS — `CONFIG::GPSConfig`

| Flag | Required | Default | Description |
|---|---|---|---|
| `--gps` | Y | `""` | Path to the GPS trajectory file. See [gps-trajectory.md](gps-trajectory.md). |
| `--gps_id` | N | `id` | Trajectory-ID column. |
| `--gps_geom` | N | `geom` | WKT geometry column (CSV-trajectory + OGR formats). |
| `--gps_x` | N | `x` | X column (CSV-point format only). |
| `--gps_y` | N | `y` | Y column (CSV-point format only). |
| `--gps_timestamp` | N | `timestamp` | Optional timestamp column. |
| `--gps_point` | N | (off) | Treat the CSV input as one-row-per-point instead of one-row-per-trajectory. |

### Result — `CONFIG::ResultConfig`

| Flag | Required | Default | Description |
|---|---|---|---|
| `--output` / `-o` | Y | `""` | Output CSV path. See [match-result-csv.md](match-result-csv.md). |
| `--output_fields` | N | `""` (means defaults) | Comma-separated list selecting which optional columns to emit. Empty value keeps the writer's defaults (`opath,cpath,mgeom`). The literal value `all` enables every column. |

The accepted token list inside `--output_fields`:

```text
opath, cpath, tpath, mgeom, pgeom,
offset, error, spdist, tp, ep, length, duration, speed,
all
```

### Runtime / logging

| Flag | Required | Default | Description |
|---|---|---|---|
| `--log_level` / `-l` | N | `2` (info) | spdlog level: 0 trace, 1 debug, 2 info, 3 warn, 4 err, 5 critical, 6 off. |
| `--step` / `-s` | N | `100` | Progress-report interval in trajectories. |
| `--use_omp` | N | (off) | Match trajectories in parallel via OpenMP. The number of threads honors `OMP_NUM_THREADS`. |
| `--help` / `-h` | N | (off) | Print help text and exit. |

## XML configuration

When a single argument ending in `.xml` is passed, the matcher routes through its `load_xml` instead of `load_arg`. The XML structure is a flat `<config>` document with three input sub-blocks, one parameters block, one output block, and one other block:

```xml
<config>
  <input>
    <network>
      <file>edges.shp</file>
      <id>id</id>
      <source>source</source>
      <target>target</target>
      <weight>NO_WEIGHT</weight>
      <turn_ban_file>NO_TURN_BANS</turn_ban_file>
    </network>
    <gps>
      <file>gps.csv</file>
      <id>id</id>
      <geom>geom</geom>
      <timestamp>timestamp</timestamp>
      <x>x</x>
      <y>y</y>
      <!-- presence of <gps_point/> enables CSV-point mode -->
    </gps>
  </input>
  <parameters>
    <!-- Matcher-specific. See e.g. STMATCHConfig::load_from_xml for STMatch's keys. -->
  </parameters>
  <output>
    <file>matched.csv</file>
    <fields>
      <!-- one <name/> child per enabled column -->
      <opath/>
      <cpath/>
      <mgeom/>
    </fields>
  </output>
  <other>
    <log_level>2</log_level>
    <step>100</step>
    <use_omp/>  <!-- presence enables OpenMP -->
  </other>
</config>
```

A `<network><file>` value is required; everything else has a default matching the CLI flag table above.

## Exit codes

| Code | Meaning |
|---|---|
| `0` | Success, or `--help` printed. |
| Non-zero | Validation failed or a critical error was thrown during load / match. |

The actual exit code for failure is what `main()` returns — generally `0` even on critical errors because the binaries return early without setting a non-zero status. **Production CI should grep the log output for `[critical]` rather than rely on exit status.**

## Polymatch superset

`polymatch` accepts every flag above plus polygon-specific additions; see `specs/001-polymatch-algorithm/contracts/polymatch-cli.md`. The CLI-parity test (T089 in that feature's tasks) verifies polymatch parses the full weightmatch flag set without error.
