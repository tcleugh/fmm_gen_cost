# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a fork of [FMM (Fast Map Matching)](https://github.com/cyang-kth/fmm) — a C++11 map matching framework that matches noisy GPS trajectories to road networks using Hidden Markov Models. This fork adds a **weightmatch** algorithm that routes on a link graph (edge-based graph with turn costs) and supports generalized cost functions via the `Edge::cost` field.

## Build Commands

```bash
# Full build (from project root)
mkdir -p build && cd build && cmake .. && make -j$(nproc)

# Rebuild after changes (from build/)
make -j$(nproc)

# Build tests (tests are excluded from default build)
cd build && cmake .. && make tests

# Run individual test executables
./build/test/algorithm_test
./build/test/fmm_test
./build/test/network_test
./build/test/network_graph_test

# Polymatch real-network validation suite (default-on; runs as part of polymatch_test).
# Use the tag to run only the real-network subset during iteration.
./build/polymatch_test '[real_network]'
```

Executables are output to `build/`: `fmm`, `stmatch`, `weightmatch`, `h3mm`, `ubodt_gen`.

## Dependencies

C++11, CMake >=3.5, GDAL >=2.2, Boost (Graph, Geometry, Serialization) >=1.56, OpenMP, SWIG (for Python bindings).

## Architecture

### Namespaces and Key Directories

- **`FMM::NETWORK`** (`src/network/`) — Road network loading (GDAL/shapefile), R-tree spatial index for candidate search, boost graph wrapper (`NetworkGraph`) with Dijkstra/A* routing.
- **`FMM::ROUTING`** (`src/network/link_graph_routing.hpp/.cpp`) — **Fork addition.** Edge-based `LinkGraph` where vertices are `EdgeIndex` values, used by `weightmatch`. Includes reusable `DijkstraState` and `IndexedMinHeap` to avoid per-query allocation.
- **`FMM::MM`** (`src/mm/`) — Map matching algorithms sharing a common HMM transition graph (`TransitionGraph`) and result types (`MatchResult`).
- **`FMM::PYTHON`** (`src/python/`) — POD types (`PyMatchResult`, `PyCandidate`) for SWIG Python bindings.
- **`FMM::CONFIG`** (`src/config/`) — XML and CLI argument parsing for network, GPS, and result output configuration.

### Map Matching Algorithms

All algorithms follow the same pattern: build candidate points via R-tree KNN search, construct a `TransitionGraph` (HMM), update transition/emission probabilities between layers, run Viterbi backtracking, then build a complete path (`C_Path`).

| Algorithm | Class | Routing | Use Case |
|-----------|-------|---------|----------|
| **FMM** | `FastMapMatch` | Precomputed UBODT lookup | Small/medium networks; requires `ubodt_gen` precomputation step |
| **STMatch** | `STMATCH` | Online Dijkstra on `CompositeGraph` (network + dummy candidate nodes) | Large networks; no precomputation |
| **WeightMatch** | `WEIGHTMATCH` | Online Dijkstra on `LinkGraph` (edge-based, uses `Edge::cost`) | Networks with turn costs or generalized edge costs |

### Core Data Flow

1. `Network` loads shapefile → builds edges with `id`, `source`, `target`, `weight`, `cost`, `geom`
2. `Network::search_tr_cs_knn()` finds candidate edges for each GPS point via R-tree
3. `TransitionGraph` constructs HMM layers from candidates, computes emission probabilities
4. Algorithm-specific routing fills transition probabilities between candidate layers
5. Viterbi `backtrack()` yields `TGOpath` → converted to `O_Path` (per-point edge) and `C_Path` (complete edge sequence)
6. `Network::complete_path_to_geometry()` builds matched geometry from `C_Path`

### Key Types (`src/network/type.hpp`, `src/mm/mm_type.hpp`)

- `NodeID`/`EdgeID` (long long) — original IDs from shapefile, can be discontinuous/negative
- `NodeIndex`/`EdgeIndex` (unsigned int) — contiguous 0-based indices used internally
- `Edge` — has both `weight` (for standard routing) and `cost` (for generalized cost routing)
- `Candidate` — a GPS point snapped to an edge, with offset and distance
- `O_Path` — optimal path (one EdgeID per GPS point)
- `C_Path` — complete path (topologically connected edge sequence)

### Python Bindings

SWIG-based bindings defined in `python/fmm.i`, built as `fmm.py` + `_fmm.so`. The Python API uses POD types from `src/python/pyfmm.hpp` to avoid exposing internal C++ objects.

### Path convention

Never hardcode an absolute path like `/workspace/...` into source, scripts, or
tests — the project root is mounted on the sandbox at `/workspace`, but it sits
at a different location outside the sandbox. Use repo-relative paths
(preferred) or have CMake inject the absolute path via
`target_compile_definitions` for binaries that need to address fixtures
regardless of cwd (`polymatch_test`'s `POLYMATCH_FIXTURE_DIR` is the canonical
example).

<!-- SPECKIT START -->
For additional context about technologies to be used, project structure,
shell commands, and other important information, read the current plan at
`specs/002-real-network-validation/plan.md`.
<!-- SPECKIT END -->
