# Contract: Base Match Result CSV

The default output written by `IO::CSVMatchResultWriter` (`src/io/mm_writer.cpp`). All matchers except `polymatch` use this writer directly; polymatch uses `IO::PolyMMWriter`, which inherits the header layout below and appends polygon-specific columns (see `specs/001-polymatch-algorithm/contracts/polymatch-output.md`).

## Encoding

- Semicolon (`;`) is the **column** delimiter — matches the GPS reader's convention.
- Comma (`,`) is the **inner-list** delimiter (per-point values within a single column).
- Pipe (`|`) splits inter-point segments within the `tpath` column only.
- UTF-8 / ASCII text. One header line, then one row per trajectory.

## Header

The header is built by `CSVMatchResultWriter::write_header()` from the `OutputConfig` flags:

```text
id[;opath][;error][;offset][;spdist][;pgeom][;cpath][;tpath][;mgeom][;ep][;tp][;length][;duration][;speed]
```

Each optional column is emitted only when its `write_*` flag is set in `OutputConfig`. **The `id` column is always first and always emitted.** `opath`, `cpath`, `mgeom` are on by default; the rest are off by default.

## Columns

| Column | Source | Type | Description |
|---|---|---|---|
| `id` | trajectory ID | int64 | Carried from the GPS input. |
| `opath` | per-point edge ID | list&lt;int64&gt; | Comma-separated. One matched edge ID per GPS point in the trajectory, in trajectory order. `-1` marks a truncated / unmatched point. |
| `error` | per-point distance | list&lt;double&gt; | One value per GPS point: perpendicular distance from the observed point to the matched point on the chosen candidate edge. |
| `offset` | per-point arc-length | list&lt;double&gt; | Distance from each matched candidate edge's source node to the matched-point projection. |
| `spdist` | per-point shortest-path | list&lt;double&gt; | One value per GPS point starting from point 2: shortest-path distance from the previous matched candidate to this one. |
| `pgeom` | matched points | WKT LineString | The line connecting the matched-point projections, one vertex per GPS point. |
| `cpath` | complete path | list&lt;int64&gt; | Topologically connected edge ID sequence covering the whole matched route. |
| `tpath` | per-point cpath segment | list&lt;int64&gt; with `,` and `\|` | For each consecutive pair of GPS points, the slice of `cpath` traversed. Inner separator `,`; pair separator `\|`. |
| `mgeom` | matched geometry | WKT LineString | The route's actual road geometry, clipped to the matched offsets at the trip endpoints. |
| `ep` | emission probability | list&lt;double&gt; | HMM emission probability per matched candidate. |
| `tp` | transition probability | list&lt;double&gt; | HMM transition probability per consecutive matched-candidate pair. |
| `length` | per-edge length | list&lt;double&gt; | Length of each matched candidate edge. |
| `duration` | per-point dt | list&lt;double&gt; | `timestamps[i] - timestamps[i-1]` from the GPS input. Empty column if input had no timestamps. |
| `speed` | per-point speed | list&lt;double&gt; | `spdist / duration` per inter-point segment. |

The columns are written in the order shown above regardless of which subset is enabled.

## Column selection

`ResultConfig::register_arg` exposes `--output_fields` (a comma-separated list parsed via `ResultConfig::string2set`) that toggles which `write_*` flags are set. Defaults (when no `--output_fields` is passed):

- on: `opath`, `cpath`, `mgeom`
- off: everything else

When `--output_fields` is given, **only** the listed columns are emitted; the defaults are not implicitly added. To enable a single extra column on top of the defaults, list all four.

## Empty / unmatched

When a trajectory produces no match (`MatchResult` has empty `opath`/`cpath`), the writer emits the `id` followed by empty fields separated by `;`. There is no special sentinel — column count is preserved.

## Thread safety

`CSVMatchResultWriter::write_result` is not internally synchronized — callers in `WEIGHTMATCHApp` / `STMATCHApp` / `FMMApp` serialize writes implicitly because the OpenMP-parallel matching path appends to a buffer that's flushed inside a `#pragma omp critical` block.

Polymatch's `PolyMMWriter` adds its own internal `std::mutex` (FR-017) — see `specs/001-polymatch-algorithm/contracts/polymatch-output.md`.

## Polymatch delta

`PolyMMWriter` writes the same header columns above, then appends:

```text
[;polygon_ids][;entry_aps][;egress_aps][;is_through][;polygon_distances]
```

…and reinterprets the integers inside `opath` / `cpath` so that negative values are polygon IDs (negated). Full details in the polymatch-output contract.
