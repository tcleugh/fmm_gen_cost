# Contract: Real-Network Trace CSV

Output of `polymatch_traces_gen` and input to the real-network validation harness.

**Lives at**: `test/data/polymatch/real_example_area/trips.csv` (committed).

**Extends**: [`specs/000-pre-existing/contracts/gps-trajectory.md`](../../000-pre-existing/contracts/gps-trajectory.md) Format 1 (CSV trajectory, semicolon-delimited).

## Encoding

- Semicolon (`;`) is the **column** delimiter.
- Comma (`,`) is the inner separator inside WKT coordinates (e.g., `LINESTRING(0 0,1 1)`).
- UTF-8 / ASCII. LF line endings (`\n`), not CRLF.
- One header line, then one data row per trace.
- Locale-independent floating-point output: `std::locale::classic()` + `std::setprecision(9)`.

## Header (exact)

```text
id;geom;category
```

The header order is fixed. Adding columns is breaking for the byte-determinism guarantee (FR-002 / SC-002) but tolerated by polymatch's `CSVTrajectoryReader` (unknown columns are silently ignored).

## Columns

| Column | Type | Description |
|---|---|---|
| `id` | int | Stable trace identifier, unique within the batch. Assigned by category in 100-wide buckets (see ID ranges below). |
| `geom` | WKT `LINESTRING` | Trace polyline. ≥ 2 points for normal categories; 0 or 1 points for `short-trip` edge cases (rendered as a degenerate WKT — see below). Coordinates in the network's projected CRS (GDA 1994 Australia Albers, meters). |
| `category` | string | One of: `link-only`, `polygon-traversal`, `polygon-shared-ap`, `mid-polygon-start`, `mid-polygon-end`, `fully-inside`, `through-routing`, `off-network-noise`, `short-trip`, `duplicate-points`. |

### ID ranges by category

| Category | ID range |
|---|---|
| `link-only` | 1000–1099 |
| `polygon-traversal` | 1100–1199 |
| `polygon-shared-ap` | 1200–1299 |
| `mid-polygon-start` | 1300–1399 |
| `mid-polygon-end` | 1400–1499 |
| `fully-inside` | 1500–1599 |
| `through-routing` | 1600–1699 |
| `off-network-noise` | 1700–1799 |
| `short-trip` | 1800–1899 |
| `duplicate-points` | 1900–1999 |

100-wide buckets give headroom; FR-003 requires ≥ 200 traces total and FR-004 requires ≥ 20 per category, so ~20 IDs per bucket are populated in practice.

### Degenerate-geometry encoding

The matcher (and `CSVTrajectoryReader`) accept any `LINESTRING(...)`. For edge-case categories:

- **0-point trace**: rendered as `LINESTRING EMPTY`. `CSVTrajectoryReader` returns a `LineString` with zero points; the matcher's `npts < 2` skip kicks in.
- **1-point trace**: rendered as `LINESTRING(x y, x y)` (the same point repeated) so the WKT remains valid. The matcher's `npts < 2` skip still applies because the geometric segment has zero length — *but* polymatch counts WKT vertices, not unique vertices, so a duplicated single point shows as 2 points. The generator MUST instead emit `LINESTRING EMPTY` for 0-point traces and a single `LINESTRING(x y, x y)` only when the test specifically wants to exercise duplicate-handling. **For `short-trip`, the generator emits 2–3 distinct points — never 0 or 1.** 0/1-point traces are reserved for `duplicate-points` (where the duplicates are interior to a longer trace) and are not produced by this generator at all.

### Row ordering

Rows MUST be sorted by `id` ascending. This makes diffs human-readable when the fixture changes and makes the byte-deterministic-output check trivial.

## Example

```text
id;geom;category
1000;LINESTRING(2053281.123456789 -3124567.987654321,2053285.001 -3124570.5);link-only
1001;LINESTRING(2055100 -3125900,2055150 -3125850,2055200 -3125800);link-only
1100;LINESTRING(2054000 -3122000,2054100 -3122050,2054120 -3122100,2054200 -3122150);polygon-traversal
1500;LINESTRING(2054500 -3120000,2054510 -3120005,2054520 -3120010);fully-inside
1600;LINESTRING(2052500 -3119000,2056500 -3119500);through-routing
1700;LINESTRING(2053000 -3120000,2053050 -3120050,3000000 -3120000,2053100 -3120050);off-network-noise
1800;LINESTRING(2055000 -3121000,2055020 -3121000);short-trip
1900;LINESTRING(2054000 -3123000,2054001 -3123000,2054001 -3123000,2054050 -3123050);duplicate-points
```

## Determinism guarantee

With a fixed seed (default `2026`), running `polymatch_traces_gen` against unchanged input shapefiles produces a byte-identical CSV across machines (where C++ standard library shared-object behavior is consistent, which it is for `std::mt19937_64`, `std::ofstream`, and `std::locale::classic()`).

CI verification:

```sh
sha256sum test/data/polymatch/real_example_area/trips.csv
# Run the generator from scratch:
./polymatch_traces_gen --network test/data/polymatch/real_example_area/network.shp \
                       --polygons test/data/polymatch/real_example_area/polygons.shp \
                       --access_points test/data/polymatch/real_example_area/access_points.shp \
                       --seed 2026 \
                       --output /tmp/trips.csv
sha256sum /tmp/trips.csv
# Both hashes MUST match.
```

This check is not part of the default test suite (it requires shelling out to the generator binary), but a developer who modifies the generator MUST run it before committing.

## Validation when read by the harness

On load, the validation harness asserts:

- Header line equals `id;geom;category` exactly (modulo trailing whitespace).
- Every `id` is unique.
- Every `category` value is one of `kCategoryLabels` (data-model.md).
- Total row count ≥ 200 (SC-003).
- Each category present has ≥ 20 rows OR exactly 0 rows (R7 in research.md — categories the network topology can't support are exempt).

Any of these check failures fails the test loudly *before* matching starts.
