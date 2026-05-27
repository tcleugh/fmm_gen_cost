# Solving the build toolchain in the sandbox

Notes from the 2026-05-27 session getting polymatch to compile inside a docker_gvisor sandbox.

## The problem

The default `docker_gvisor` sandbox image has Python 3.13 but no `cmake`, no `g++`, no GDAL/Boost dev headers, no OpenMP. Running as `devuser` (uid 1010), no sudo, no apt access. So we can't build the C++ project as-is.

## The fix — sandbox profile patch

Added `dockerfile_patches.user` to `/workspace/.vlc-claudit/profiles/docker_gvisor.toml`:

```toml
[impl.docker_gvisor.dockerfile_patches]
user = '''
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake g++ build-essential libgdal-dev libboost-all-dev libomp-dev swig pkg-config \
 && rm -rf /var/lib/apt/lists/*
USER devuser
'''
```

**Confirmed working after sandbox restart on 2026-05-27.** Toolchain present:

```
/usr/bin/cmake
/usr/bin/g++
/usr/bin/pkg-config
libgdal.so.36 (Debian bookworm package)
libboost_serialization.so.1.83.0
libgomp.so.1
swig 4.3.0
python 3.13.5
```

The `USER root` / `USER devuser` switches inside the patch worked — no need for the conda fallback.

## Source-side fixes also required

The toolchain patch was necessary but not sufficient. The upstream FMM code was written against older Boost (~1.56) and Python (with `distutils`), and doesn't compile against bookworm's Boost 1.83 + Python 3.13 out of the box. Two source-level fixes:

### 1. Python 3.13 dropped `distutils`

`python/CMakeLists.txt:20-22` called:

```cmake
execute_process(
  COMMAND python -c "from distutils.sysconfig import get_python_lib;print(get_python_lib())"
  ...)
```

`distutils` was removed in Python 3.12+. The call now errors silently, leaving `PYTHON_SITE_PACKAGES` empty → `install TARGETS given no LIBRARY DESTINATION`.

**Durable fix:** use `sysconfig` (stdlib since Python 2.7, replaces distutils for this purpose) and `${PYTHON_EXECUTABLE}` (already resolved by `find_package(PythonInterp)` above):

```cmake
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
  ...)
```

### 2. Boost 1.83 geometry utilities need ≥ C++14

`boost::geometry::util::remove_cptrref` (used by `boost::geometry::tag<>`, called from every `bg::wkt()` / `bg::distance()` / etc.) is defined with `std::remove_const_t` and friends in `/usr/include/boost/geometry/util/type_traits_std.hpp:113`. Those `_t` alias templates are C++14. Under `-std=c++11`, the struct has no nested `type` and you get cascading errors like:

```
no type named 'type' in 'struct boost::geometry::util::remove_cptrref<const ...>'
```

The project was set to C++11 in `CMakeLists.txt:30`:

```cmake
set(CMAKE_CXX_STANDARD 11)
```

**Note:** a `-DCMAKE_CXX_STANDARD=14` cmd-line override does NOT work — the `set()` in CMakeLists.txt is unconditional and overwrites the cache value.

**Durable fix:** bump the project to C++14 (no source changes needed — the FMM code is already C++14-compatible):

```cmake
set(CMAKE_CXX_STANDARD 14)
```

C++14 is the minimum to compile Boost 1.83 geometry. C++17 also works.

### 3. GDAL ≥ 3 returns `const OGRSpatialReference*`

`src/network/network.cpp:134` had:

```cpp
OGRSpatialReference *ogrsr = ogrFDefn->GetGeomFieldDefn(0)->GetSpatialRef();
```

GDAL 3.x changed `GetSpatialRef()` to return `const OGRSpatialReference*`. Fix is to qualify the local with `const` — we only call `GetEPSGGeogCS()` on it, which is const-safe.

### 4. Missing include in new `poly_mm_writer.cpp`

`src/io/poly_mm_writer.cpp` streams `O_Path` / `C_Path` (both `std::vector<long long>`) via `operator<<`, defined in `src/util/util.hpp`. The file was missing `#include "util/util.hpp"`. `mm_writer.cpp` had it; the new file didn't.

## Build verification (post-fix)

```bash
rm -rf /workspace/build  # old build dir referenced /home/tessa/git/fmm_gen_cost from a different host
mkdir -p /workspace/build && cd /workspace/build && cmake .. && make -j$(nproc)
```

**Confirmed clean build on 2026-05-27.** Outputs in `/workspace/build/`:

- `fmm`, `stmatch`, `weightmatch`, `polymatch`, `h3mm`, `ubodt_gen` (all six executables)
- `_fmm.so` + `python/fmm.py` (SWIG Python bindings)

`./polymatch --help` runs and prints the expected argument list.

## Test suite status

`make tests` builds **algorithm_test, polymatch_test, weightmatch_test** cleanly after one more change:

- `test/CMakeLists.txt`: `polymatch_test` was missing `$<TARGET_OBJECTS:WEIGHTMATCH_OBJ>` in its link list. POLYMATCH's algorithm calls into WEIGHTMATCH so the object library must be linked.
- `third_party/catch2/catch.hpp`: glibc ≥ 2.34 made `MINSIGSTKSZ` a runtime `sysconf()` call, breaking Catch2 v2.11.1's `static constexpr std::size_t sigStackSize = 32768 >= MINSIGSTKSZ ? ...;`. Replaced with a constant 32KB (matches upstream Catch2 fix in v2.13.7). For a more durable fix, bump `third_party/catch2/catch.hpp` to v2.13.10 — leaving for another day.

**Three upstream test files still fail to compile** — `network_test.cpp`, `fmm_test.cpp`, `network_graph_test.cpp` — but this is **pre-existing breakage on the polymatch branch**, not a toolchain issue. The polymatch branch added a required `turn_ban_filename` parameter to `Network::Network(...)` (network.hpp:62) but those upstream tests still call the one-arg form. Unrelated to the sandbox restoration; out of scope for this doc.

## What didn't work (early-session fallback attempts)

Tried Miniforge + mamba in `~/miniforge3` (rootless install path, matches the upstream example pattern). Two attempts:

1. **First attempt**: `mamba create -n fmm -c conda-forge cmake make gxx_linux-64 gcc_linux-64 gdal libgdal boost-cpp boost llvm-openmp swig` — libmamba's SAT solver spun on CPU for ~17 minutes resolving conda-forge. Took ~810MB of package cache. Eventually completed env resolution (371.7s logged), then started a second resolve cycle. Killed before the env directory was created.
2. **Second attempt** with warm cache, dropped `swig`: also got stuck in resolve, killed by user request.

Miniforge is a viable rootless path in principle, but the cold solve is too slow to be practical for an interactive session. The apt-based profile patch is much better.

## Branches and state

- Branch `001-polymatch-algorithm` has the polymatch skeleton (Phase 1 + Phase 2 from `specs/001-polymatch-algorithm/tasks.md` — T001, T002, T004–T019 marked `[X]`; T003 Python fixture generator still pending).
- With the two source-side fixes above, the skeleton compiles cleanly.
