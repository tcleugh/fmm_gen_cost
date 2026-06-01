# Phase 1 Data Model: Real-Network Validation

This feature does **not** add new entities to the matcher's data model. It does introduce three small data structures (one persisted to disk, two transient in-test) that the generator and harness use.

---

## Entity: `TraceCategory`

A label-only enum / string set capturing the ten categories from spec FR-004.

| Label | When it applies |
|---|---|
| `link-only` | Trajectory has no polygon within `radius` of any GPS point; all candidates are link candidates. |
| `polygon-traversal` | Trajectory enters at least one polygon via an AP and exits via a (possibly different) AP; ≥ 1 GPS observation inside. |
| `polygon-shared-ap` | Trajectory crosses two polygons that share an AP; matched cpath includes both polygon IDs with the shared AP serving as egress of one + entry of the next. |
| `mid-polygon-start` | First GPS observation is inside a polygon; matched first segment is that polygon with `entry_ap == kNoAccessPoint`. |
| `mid-polygon-end` | Symmetric — last GPS inside a polygon; last segment has `egress_ap == kNoAccessPoint`. |
| `fully-inside` | Every GPS observation inside one polygon; matched result has a single polygon segment with both APs absent. |
| `through-routing` | The optimal `PolyLinkGraph` route between two distant link endpoints crosses one polygon as a shortcut, but no GPS observation falls inside that polygon. |
| `off-network-noise` | At least one interior GPS point lies >> `radius` from any network edge. |
| `short-trip` | Trajectory has 2 or 3 GPS points. |
| `duplicate-points` | At least one pair of consecutive GPS points has identical coordinates. |

Concrete representation:

```cpp
constexpr const char* kCategoryLabels[] = {
  "link-only", "polygon-traversal", "polygon-shared-ap",
  "mid-polygon-start", "mid-polygon-end", "fully-inside",
  "through-routing", "off-network-noise", "short-trip", "duplicate-points"
};
```

Stored as a `std::string` in the harness's parallel-vector lookup; written as-is to the CSV's `category` column.

---

## Entity: `GeneratedTrace`

Transient — produced by the generator, immediately serialized to CSV, not kept in memory by the harness.

| Field | Type | Description |
|---|---|---|
| `id` | `int` | Stable trace identifier; unique within the batch (FR-005). Allocated by category in 100-wide buckets (link-only 1000–1099, polygon-traversal 1100–1199, …). |
| `geom` | `CORE::LineString` | The GPS polyline. ≥ 2 points (or 0/1 for short/edge-case categories). |
| `category` | `std::string` | One of `kCategoryLabels`. |

Persisted as one CSV row: `<id>;LINESTRING(x1 y1, x2 y2, ...);<category>`.

---

## Entity: `Invariant`

Each is a free function in the harness anonymous namespace:

```cpp
using InvariantFn = std::function<std::optional<std::string>(
    const PolyMatchResult& pm,
    const std::optional<MatchResult>& wm_link_only,  // populated only for link-only category
    const Trajectory& traj,
    const std::string& category,
    const Network& net,
    const PolygonLayer& poly,
    const AccessPointLayer& aps)>;
```

Returns `std::nullopt` on pass; a one-line human-readable failure message on fail (e.g., `"cpath[3]=42 -> cpath[4]=-7 but polygon 7 has no AP at edge 42's endpoint"`).

Four invariants for v1, one per FR-011..FR-014:

| ID | Spec FR | Applies to | Pass condition |
|---|---|---|---|
| `cpath-topology` | FR-011 | every trace | Every consecutive link-link pair in `cpath` shares a network node; every link↔polygon transition uses one of that polygon's APs. |
| `is-through-has-aps` | FR-012 | every trace | Every `PolygonSegment` with `is_through == true` has `entry_ap != kNoAccessPoint && egress_ap != kNoAccessPoint`. |
| `link-only-eq-weightmatch` | FR-013 | traces with `category == "link-only"` | `pm.base.opath == wm.opath && pm.base.cpath == wm.cpath`. |
| `distance-inside-finite` | FR-014 | every trace | For every `PolygonSegment`, `distance_inside` is finite and `>= 0.0`. |

---

## Entity: `ViolationLedger`

Transient harness helper that aggregates per-invariant failures.

```cpp
struct ViolationLedger {
  struct Entry { int trace_id; std::string reason; };
  std::map<std::string /*invariant_id*/, std::vector<Entry>> per_invariant;
  std::map<std::string /*invariant_id*/, size_t> pass_count;

  void record_pass(const std::string& invariant_id);
  void record_fail(const std::string& invariant_id, int trace_id, const std::string& reason);
  void print_summary(std::ostream& os) const;   // one line per invariant + up to 10 trace IDs
  bool any_failures() const;
};
```

At end of TEST_CASE the harness calls `ledger.print_summary(std::cerr)` for visibility, then `CHECK(!ledger.any_failures())` to gate the test.

---

## Relationships

```text
TraceCategory (label)
       ▲
       │ stamped
       │
GeneratedTrace ──serialized──► trips.csv ──read──► Trajectory + category-string
                                                          │
                                                          ▼
                                                  POLYMATCH::match_traj  ──► PolyMatchResult
                                                          │                          │
                                                          ▼                          │
                                                  WEIGHTMATCH::match_traj  ─► MatchResult (link-only)
                                                          │                          │
                                                          └────────► Invariant ◄─────┘
                                                                         │
                                                                         ▼
                                                                   ViolationLedger
                                                                         │
                                                                         ▼
                                                                  Catch2 CHECK
```

---

## Invariants the data model maintains

| Invariant | Where enforced |
|---|---|
| Every trace ID is unique within the batch (FR-005). | Generator's per-category id-range allocation; cross-checked by harness load (asserts no duplicates). |
| Every category present in the CSV is one of the ten in `kCategoryLabels`. | Harness sanity check on CSV load. |
| Categories with zero traces are skipped, not failed (R7). | `ViolationLedger::print_summary` emits "skipped — no traces produced" for an empty per-invariant tally on the affected category-filter invariants. |
| The CSV is byte-deterministic for a fixed seed (FR-002 / SC-002). | Generator's `std::mt19937_64` + `setprecision(9)` + `locale::classic()`. Verifiable by `sha256sum` across two runs. |
