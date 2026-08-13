# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FMM (Fast Map Matching) is a C++ framework that matches noisy GPS traces to a road
network using hidden Markov models (HMM) + precomputation. Upstream provides three
matchers (FMM, STMATCH, H3); **this fork adds a fourth, WeightMatch**, which is where
nearly all local development happens (see `src/mm/weightmatch/`, `src/network/link_graph_routing.*`,
turn-ban support in `src/network/`, and `test/weightmatch_test.cpp`).

Each program has both a Python binding (SWIG) and a command-line executable, and each
accepts either an XML config file (`prog config.xml`) or CLI arguments.

## Build

Dependencies: GDAL ≥2.2, Boost (serialization; +exception on Windows) ≥1.56, OpenMP,
CMake ≥3.5, SWIG (for Python bindings). See README for the Ubuntu `apt-get` line.

```bash
mkdir build && cd build   # the dir MUST be named `build`
cmake ..
make -j4
sudo make install         # installs executables to /usr/local/bin + Python bindings to site-packages
```

Build gotchas that will bite you:
- **In-source builds are disabled** (`CMAKE_DISABLE_IN_SOURCE_BUILD ON`). You must build in a subdir.
- The runtime output dir is **hardcoded** to `${SOURCE}/build` (`CMAKE_RUNTIME_OUTPUT_DIRECTORY`).
  All executables land in `build/` regardless of where you invoke cmake, and tests assume this — use `build/`.
- Build type is forced to `Release` (`-O3`, C++11, `SPDLOG_ACTIVE_LEVEL=TRACE`).
- Sources are picked up by `file(GLOB ...)` per subdirectory, then compiled into per-directory
  `OBJECT` libraries that link into one `FMMLIB` shared lib. **Adding a new `.cpp` requires re-running
  `cmake`** so the glob re-evaluates.

Executables produced: `fmm`, `stmatch`, `weightmatch`, `h3mm`, `ubodt_gen`.

## Test

Tests use **Catch2 v2** (`third_party/catch2/catch.hpp`) and are `EXCLUDE_FROM_ALL` — build them explicitly:

```bash
cd build && make tests -j$(nproc)
```

Executables (`algorithm_test`, `network_test`, `network_graph_test`, `fmm_test`, `weightmatch_test`)
build into `build/`, but **must be run from `build/test/`** because they resolve test-data paths
relative to the project root (two directories up from the CWD):

```bash
cd build/test
../weightmatch_test                          # all cases
../weightmatch_test "[weightmatch][sa4_212]" # filter by tag
../weightmatch_test -s                        # verbose
../weightmatch_test --list-tests
```

WeightMatch golden files live in `test/data/weightmatch/` as `id;opath;cpath` (semicolon-separated).
`sa4_212` (3 826 edges) and `sa4_17` (65 270 edges) are real networks with basic/edge-case/stress/turn-ban
fixtures. If you change matching behavior or add trips, **regenerate goldens** with
`python test/generate_weightmatch_test_data.py` (needs `geopandas`). Full details in `test/README.md`.

## Architecture

### The four matchers and how they differ

All matchers share the same HMM skeleton: build per-point candidate sets (KNN over an R-tree),
score **emission** probability from GPS error and **transition** probability from
shortest-path-vs-Euclidean distance in a `TransitionGraph`, then Viterbi-backtrack for the optimal
path. They differ only in how shortest paths between candidates are computed:

- **FMM** (`src/mm/fmm/`): precomputes an **UBODT** (Upper-Bounded Origin-Destination Table) with the
  separate `ubodt_gen` program, then does O(1) table lookups. Fast for small/mid networks. Node-based.
- **STMATCH** (`src/mm/stmatch/`): no precomputation; routes on the fly over the node-based graph.
  Scales to large networks.
- **WeightMatch** (`src/mm/weightmatch/`): the local addition. Routes on an **edge-based link graph**
  (see below) so it can honor **turn bans**, **turn costs**, and **generalized per-edge weights**
  (a shapefile `--weight` field, falling back to geometry length). This is the matcher to reach for
  when the network has directional/turn constraints.
- **H3** (`src/mm/h3mm/`, header-only): matches to Uber H3 hexagons instead of edges.

### Two graph representations — the key mental model

- **Node-based** (`src/network/network_graph.*`, `bidirectional_network_graph.*`): a Boost graph whose
  vertices are road nodes and edges are road edges. Used by FMM/STMATCH. Cannot express turn costs
  (a turn is a node, and node-based Dijkstra can't distinguish which edge you arrived on).
- **Edge-based link graph** (`src/network/link_graph_routing.*`): vertices **are** road edges
  (`EdgeIndex`), and an arc connects edge A→edge B iff they're topologically adjacent and not turn-banned.
  This is what lets WeightMatch model turns. It ships its own hand-rolled Dijkstra
  (`shortest_edge_to_edges`) with:
  - `IndexedMinHeap` — a decrease-key binary heap.
  - `DijkstraState` — epoch-stamped `dist`/`parent`/`seen` arrays so state is **reused across calls
    without reallocation** (`next_epoch()` is an O(1) logical clear). One `DijkstraState`+`IndexedMinHeap`
    pair is created per OpenMP thread and threaded through `match_traj` — do not allocate these per call.
  - An optional **upper-bound pruning** factor (`upper_bound_factor`): once enough goals are found,
    exploration stops at `max_found_cost * factor`. `0` disables it.

### Shared pieces

- `src/network/network.*` — loads the shapefile into `edges`/nodes, builds the R-tree, does KNN
  candidate search (`search_tr_cs_knn*`, including a fallback-radius variant), and owns the turn-ban set.
- `src/mm/transition_graph.*` — the HMM layers + Viterbi; `calc_ep`/`calc_tp` are the probability models.
- `src/mm/mm_type.hpp` — shared result types. Note `opath` (one matched edge per GPS point) vs
  `cpath` (the full topologically-connected edge sequence traversed); `indices` maps opath→cpath.
- `src/core/` geometry/GPS types, `src/io/` GPS reader + result writer, `src/config/` XML+CLI parsing,
  `src/network/type.hpp` core `Edge`/`Turn`/index typedefs.

### Per-program code layout (consistent pattern)

Every matcher is three files plus a thin `main`:
- `X_algorithm.{hpp,cpp}` — the matching logic and its `XConfig` (algorithm params).
- `X_app.{hpp,cpp}` — `XApp` that loads network/graph once and loops over trajectories
  (serial, or `#pragma omp parallel` when `--use_omp`).
- `X_app_config.{hpp,cpp}` — `XAppConfig` parsing XML or CLI into I/O + algorithm config.
- `src/app/X.cpp` — `main()` that just builds the AppConfig, constructs the App, and calls `run()`.

When adding a matcher option, it typically must be threaded through all four: the algorithm's
`Config` struct + its `register_arg`/`load_from_arg`/`register_help`, and the app config's XML path.

`src/fmm-api.hpp` is the umbrella header for the installed C++ library; `src/python/pyfmm.hpp` +
`python/fmm.i` drive the SWIG bindings.

## Repo conventions

- **Commit as you go.** After each meaningful, self-contained change, make a descriptive commit rather
  than batching many changes into one. Don't push or open PRs unless asked.
- Branches follow `tcleugh/tessa/<feature>`; changes land via PR merges to `master`.
- `brainstorms/` holds design docs (e.g. GPS-accuracy explorations); `.claude/skills/` provides the
  `brainstorm` and `scopestorm` workflows used to produce them.
- CI (`.travis.yml`) builds all executables, installs, smoke-runs each binary with no args, and runs
  `example/python/fmm_test.py`.
