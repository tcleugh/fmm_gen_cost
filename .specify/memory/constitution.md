<!--
## Sync Impact Report

**Version change**: [NEW — 1.0.0 initial ratification]

**Principles established**:
- I. Performance-First (new)
- II. Algorithmic Correctness (new)
- III. Test-Driven Development (new)
- IV. Dual-API Stability (new)
- V. Separation of Concerns (new)

**Added sections**:
- Core Principles (5 principles)
- Technical Constraints
- Development Workflow
- Governance

**Removed sections**: N/A (initial version)

**Templates reviewed**:
- ✅ `.specify/templates/plan-template.md` — Constitution Check gate present; compatible with principles
- ✅ `.specify/templates/spec-template.md` — User story structure compatible; no conflicts
- ✅ `.specify/templates/tasks-template.md` — Task phases compatible; no conflicts
- ✅ `.specify/templates/checklist-template.md` — No conflicts identified

**Deferred TODOs**: None
-->

# FMM Fork (WeightMatch) Constitution

## Core Principles

### I. Performance-First

C++ performance constraints are non-negotiable in hot paths. Specifically:

- Routing queries MUST NOT allocate heap memory per query; reuse `DijkstraState` and `IndexedMinHeap` across calls.
- R-tree candidate search MUST remain the sole spatial index for KNN lookups — do not add secondary indexes without
  benchmarked justification.
- OpenMP parallelism MUST be preserved for trajectory-level parallel matching; no synchronization primitives that
  serialize the outer loop.
- Algorithmic complexity of any new routing or matching step MUST be documented with its Big-O bound in the PR
  description.

**Rationale**: Map matching at scale (millions of GPS points, million-edge networks) is the primary value proposition.
Regressions in throughput are not acceptable without explicit user opt-in.

### II. Algorithmic Correctness

All map matching algorithms MUST produce results that satisfy HMM and topological invariants:

- Every `C_Path` MUST be topologically connected — consecutive edges MUST share a node.
- Transition probabilities MUST be non-negative and finite; infinite/NaN values in the HMM graph are a crash-class bug.
- Viterbi backtracking MUST terminate; cycles in the `TransitionGraph` are forbidden.
- Duplicate GPS points (zero-distance observations) MUST be handled gracefully without producing degenerate HMM layers.

**Rationale**: Incorrect matched paths silently corrupt downstream analyses. Correctness failures are worse than
performance failures because they may go undetected.

### III. Test-Driven Development (NON-NEGOTIABLE)

All changes to algorithmic or routing code MUST be accompanied by tests before merging:

- New algorithm variants MUST have unit tests covering: normal trajectories, edge trajectories (single point, duplicate
  points, out-of-network points), and at least one integration test against a known shapefile fixture.
- Bug fixes MUST include a regression test that reproduces the original failure.
- Tests are built separately (`make tests`) and all four suites (`algorithm_test`, `fmm_test`, `network_test`,
  `network_graph_test`) MUST pass clean before any PR is merged.
- Test code follows the same C++11 standard as production code.

**Rationale**: The fork introduced the `weightmatch` algorithm without full test parity at the time. This principle
closes that gap and prevents future regressions.

### IV. Dual-API Stability

The C++ API and the Python SWIG bindings MUST remain consistent:

- Public C++ types exposed through SWIG (`src/python/pyfmm.hpp`) MUST be POD structs — no internal C++ objects
  (smart pointers, STL containers with non-trivial destructors) may cross the language boundary.
- Removing or renaming a public C++ method or type referenced in `python/fmm.i` requires a corresponding update to
  the `.i` file and a MAJOR version bump of any released Python package.
- The Python API MUST be validated end-to-end (import + match a single trajectory) before any release.

**Rationale**: Downstream users interact with this project primarily through the Python API. Breaking it silently is
unacceptable.

### V. Separation of Concerns

Namespace boundaries MUST be respected and not blurred:

- `FMM::NETWORK` owns spatial loading, R-tree, and graph construction. It MUST NOT contain HMM logic.
- `FMM::ROUTING` owns routing algorithms (`LinkGraph`, `DijkstraState`, `IndexedMinHeap`). It MUST NOT depend on
  `FMM::MM` types.
- `FMM::MM` owns HMM construction, emission/transition probability computation, Viterbi decoding, and result assembly.
  It MAY depend on `FMM::NETWORK` and `FMM::ROUTING` but not vice versa.
- `FMM::CONFIG` and `FMM::PYTHON` are leaf namespaces — they may depend on core namespaces but MUST NOT be depended
  on by them.

**Rationale**: Clean namespace layering enables isolated testing and prevents the codebase from collapsing into a ball
of mud as new algorithms are added.

## Technical Constraints

- **Language**: C++11 strictly. No C++14/17/20 features without explicit team agreement and compiler-support
  verification.
- **Build system**: CMake >=3.5. Build artifacts go to `build/`; no in-source builds.
- **Required dependencies**: GDAL >=2.2, Boost (Graph, Geometry, Serialization) >=1.56, OpenMP. These MUST NOT be
  made optional without a feature-flag mechanism.
- **Optional dependencies**: SWIG (Python bindings only). Builds without SWIG MUST still produce all C++ executables.
- **Executables**: `fmm`, `stmatch`, `weightmatch`, `h3mm`, `ubodt_gen`. All five MUST build cleanly on Linux and macOS.
- **No dynamic polymorphism in hot paths**: Virtual dispatch in the HMM inner loop is prohibited. Use templates or
  static dispatch instead.

## Development Workflow

- All feature work MUST start from a `specs/` specification (via `/speckit-specify`) before implementation.
- Implementation plans (`plan.md`) MUST pass the Constitution Check gate before Phase 0 research proceeds.
- PRs targeting `weightmatch` or `LinkGraph` code MUST include benchmark numbers if any routing path is modified.
- All four test suites MUST pass before a PR may be merged: `algorithm_test`, `fmm_test`, `network_test`,
  `network_graph_test`.
- Commit messages MUST reference the relevant algorithm or namespace (e.g., `weightmatch:`, `network:`, `mm:`).
- `CLAUDE.md` is the authoritative runtime development guidance document for AI-assisted work in this repository.

## Governance

This constitution supersedes all other informal practices and conventions. Amendments require:

1. A written rationale explaining what principle is being added, modified, or removed.
2. A version bump following semantic versioning:
   - MAJOR — principle removal or backward-incompatible redefinition.
   - MINOR — new principle or materially expanded guidance.
   - PATCH — clarifications, wording fixes, non-semantic refinements.
3. The `Last Amended` date MUST be updated to the ISO date of the amendment.
4. All dependent templates (`.specify/templates/`) MUST be reviewed for consistency after any MINOR or MAJOR amendment.

All PRs and code reviews MUST verify compliance with the principles above. Violations must be called out in review and
resolved before merge — not deferred to follow-up issues.

**Version**: 1.0.0 | **Ratified**: 2026-05-26 | **Last Amended**: 2026-05-26
