# Test Suite

## Overview

The test suite uses [Catch2 v2](https://github.com/catchorg/Catch2/tree/v2.x). Tests are excluded from the default build and must be built explicitly.

## Building the Tests

```bash
mkdir -p build && cd build
cmake ..
make tests -j$(nproc)
```

This builds all test executables into the `build/` directory:

| Executable | Tests |
|---|---|
| `algorithm_test` | Geometry and core algorithm utilities |
| `network_test` | Network loading and R-tree candidate search |
| `network_graph_test` | Boost graph routing on the node-based network |
| `fmm_test` | FMM (UBODT-based) map matching |
| `weightmatch_test` | WeightMatch (link-graph-based) map matching |

## Running Tests

Test executables use relative paths anchored two directories above their working directory (i.e., the project root). They must be run from `build/test/`:

```bash
# Run all weightmatch tests
cd build/test && ../weightmatch_test

# Run with verbose output
cd build/test && ../weightmatch_test -s

# Run only a specific tag
cd build/test && ../weightmatch_test "[weightmatch][sa4_212]"

# List all test cases without running them
cd build/test && ../weightmatch_test --list-tests
```

Other test executables follow the same pattern:

```bash
cd build/test && ../algorithm_test
cd build/test && ../network_test
cd build/test && ../fmm_test
```

## WeightMatch Tests

`weightmatch_test.cpp` contains the following test cases:

| Test Case | Tag | Description |
|---|---|---|
| `weightmatch basic matching` | `[weightmatch]` | Basic and consistency checks on the bundled `example/data/edges.shp` synthetic network |
| `weightmatch edge cases` | `[weightmatch]` | Edge-case trajectories (very short, single-point, no-candidate) on the synthetic network |
| `weightmatch SA4=212 small real network` | `[weightmatch][sa4_212]` | Basic, edge-case, turn-ban, and consistency checks on a real 3 826-edge network |
| `weightmatch SA4=17 medium real network` | `[weightmatch][sa4_17]` | Same suite on a real 65 270-edge network |
| `weightmatch SA4=212 stress (1000 trips)` | `[weightmatch][sa4_212][stress]` | 1 000 synthetic trips matched against a pre-computed golden file |
| `weightmatch SA4=17 stress (500 trips)` | `[weightmatch][sa4_17][stress]` | 500 synthetic trips matched against a pre-computed golden file |
| `weightmatch allow_truncation` | `[weightmatch][truncation]` | Verifies `allow_truncation=true/false` behaviour for all combinations of matchable/unmatched head, tail, and gap patterns |

#### allow_truncation behaviour summary

Points far from the network (`_`) produce no candidates. The sections cover:

| Pattern | allow_truncation=true | allow_truncation=false |
|---|---|---|
| `v v v v v` | match (all points) | match |
| `_ _ v v v` | match (head truncated, opath prefixed with -1s) | empty |
| `v v v _ _` | match (tail truncated, opath suffixed with -1s) | empty |
| `_ v v v _` | match (both ends truncated) | empty |
| `_ _ _ _ _` | empty (no matchable points) | empty |
| `v _ _ _ _` | empty (only 1 matchable point) | empty |
| `_ _ _ _ v` | empty (only 1 matchable point) | empty |
| `_ _ v _ _` | empty (only 1 matchable point) | empty |
| `v _ _ _ v` | empty (gap in middle — non-contiguous) | empty |
| `_ v v _ v` | empty (gap in middle — non-contiguous) | empty |
| `_ v _ v _` | empty (gap in middle — non-contiguous) | empty |

When truncation succeeds, the returned `opath` has the same length as the input trajectory, with `-1` for each unmatched GPS point at the head or tail.

### Test Data

Test data lives under `test/data/weightmatch/`:

```
test/data/weightmatch/
├── trips_basic.csv          # synthetic trips on example network
├── trips_edge_cases.csv     # edge-case trips on example network
├── expected_basic.csv       # golden opath/cpath for basic trips
├── expected_edge_cases.csv  # golden opath/cpath for edge-case trips
├── sa4_212/
│   ├── links.shp            # real 3826-edge road network
│   ├── turn_bans.csv        # turn restrictions for sa4_212
│   ├── trips_basic.csv / trips_edge_cases.csv / trips_stress.csv
│   └── expected_basic.csv / expected_basic_tb.csv / expected_edge_cases.csv / expected_stress.csv
└── sa4_17/
    ├── links.shp            # real 65270-edge road network
    ├── turn_bans.csv
    ├── trips_basic.csv / trips_edge_cases.csv / trips_stress.csv
    └── expected_basic.csv / expected_basic_tb.csv / expected_edge_cases.csv / expected_stress.csv
```

Golden files are CSV with columns `id;opath;cpath` (semicolon-separated). `opath` is a comma-separated list of matched edge IDs (one per GPS point); `cpath` is the full topologically-connected edge sequence.

### Regenerating Test Data and Golden Files

If you change the matching algorithm or add new test trips, regenerate the golden files using the Python helper:

```bash
# Requires: geopandas (pip install geopandas)
python test/generate_weightmatch_test_data.py --help
```

The script generates synthetic GPS trips from a shapefile network and writes both the trip CSV files and corresponding expected output files that the C++ tests compare against.
