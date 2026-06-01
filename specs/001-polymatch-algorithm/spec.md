# Feature Specification: PolyMatch — Link-Polygon Map Matching Algorithm

**Feature Branch**: `001-polymatch-algorithm`

**Created**: 2026-05-26

**Status**: Draft

**Input**: User description: "build a new matching algorithm (based on weightmatch) called polymatch. It will also read in a polygon layer treated as first class objects for routing. Trips can be routed through a combination of links and polygons"

## Clarifications

### Session 2026-05-26

- Q: How should the system handle polygon ID 0 at load? → A: Reject polygon ID 0 at load with a descriptive error (ID 0 is reserved by the output format).
- Q: How should the loader handle two polygons sharing the same ID? → A: Halt with descriptive error listing the offending ID.
- Q: For polygons with holes (inner rings), how is "inside" defined? → A: Holes excluded — standard Boost.Geometry semantics. A point inside a hole is treated as outside the polygon. For multipolygons, inside any sub-polygon counts as inside.
- Q: Should the loader detect and reject CRS mismatches between the network and the polygon / access point shapefiles? → A: No — CRS is not checked. The system trusts the user to provide consistent CRSes.
- Q: How should polymatch handle a trajectory with 0 or 1 GPS points? → A: Warn and skip the trajectory; emit no matched-path output row; continue processing remaining trajectories in the batch.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Match GPS Trip Through Mixed Link-Polygon Network (Priority: P1)

A transport analyst has a road network shapefile and a separate polygon shapefile representing traversable areas such as parking lots, ferry terminals, or pedestrian plazas. They run polymatch with both layers as input and receive a matched path that seamlessly routes through road links and polygon areas in a single continuous result.

**Why this priority**: This is the core capability. Without it the feature delivers no value. Every other story depends on it.

**Independent Test**: Given a GPS trajectory that enters a polygon area, crosses it, and exits onto a different road link, running polymatch produces a matched path that correctly includes the polygon traversal with no topological gap between the entry link, the polygon, and the exit link.

**Acceptance Scenarios**:

1. **Given** a road network and polygon layer are configured, **When** a GPS trajectory is submitted that passes through a polygon area, **Then** the matched output includes the polygon as a named segment with its entry and egress access points recorded alongside road links.
2. **Given** a GPS trajectory that stays entirely on road links, **When** polymatch processes it, **Then** the result is identical to a weightmatch result on the same network (no regression).
3. **Given** a GPS trajectory that stays entirely within a single polygon, **When** polymatch processes it, **Then** the matched output identifies that polygon as the sole matched segment; entry and egress access points are absent because the trip has no link transitions.
4. **Given** a trajectory that crosses two adjacent polygons without returning to a link between them, **When** polymatch processes it, **Then** the matched output includes both polygons as consecutive segments connected via their shared access point.
5. **Given** a GPS trajectory where the optimal path routes through a polygon but no GPS observations fall inside it, **When** polymatch processes it, **Then** the matched output includes the polygon as a through-routing segment with both entry and egress access points recorded and the through-routing cost applied.

---

### User Story 2 - Configure Polygon Layer as Routing Input (Priority: P2)

A data engineer configures polymatch by pointing it at a polygon shapefile and an access point shapefile alongside the road network. The polygon shapefile carries an ID field and an optional generalized cost field. The access point shapefile contains one feature per polygon-access-point relationship: each feature is a point on a polygon boundary with a node number and a polygon ID. An access point shared between multiple polygons appears as multiple features with the same node number, one per polygon. Polymatch loads and validates all three inputs before any matching begins.

**Why this priority**: Without polygon layer loading and access point configuration, no mixed-mode matching can occur. This is independently testable by verifying load success/failure before running any trajectories.

**Independent Test**: Running polymatch with a valid polygon shapefile and access point shapefile produces a confirmation that N polygons and M access points were loaded; running with an invalid path or invalid shapefile produces a clear, descriptive error before any matching is attempted.

**Acceptance Scenarios**:

1. **Given** a valid polygon shapefile and access point shapefile are specified, **When** polymatch starts, **Then** all polygons and their access points are loaded and available for routing.
2. **Given** a polygon shapefile with a cost field is provided, **When** polymatch routes through polygons, **Then** the cost field value influences the match result in the same way edge cost does for links.
3. **Given** an invalid or missing polygon shapefile or access point shapefile is specified, **When** polymatch starts, **Then** a descriptive error is produced and matching does not proceed.
4. **Given** no polygon layer is specified, **When** polymatch starts, **Then** it falls back to link-only matching equivalent to weightmatch behavior.
5. **Given** a polygon exists in the polygon shapefile but has no features in the access point shapefile referencing its ID, **When** polymatch starts, **Then** that polygon is treated as unreachable and a warning is emitted.
6. **Given** a THROUGH_PENALTY_FACTOR value is supplied via configuration, **When** polymatch routes a trajectory through a polygon with no GPS points inside it, **Then** the through-routing cost reflects the supplied factor; changing the factor to a higher value makes through-routing less preferred relative to alternative paths.

---

### User Story 3 - Inspect Matched Path Distinguishing Links from Polygons (Priority: P3)

A researcher receives a polymatch result and needs to identify which segments of the matched path are road links and which are polygon traversals, so they can apply different analysis logic to each type.

**Why this priority**: Consumers of the output need to distinguish path types for downstream analysis. Independently testable against a known fixture with a ground-truth polygon traversal.

**Independent Test**: A matched path for a trajectory known to cross one polygon produces output where the polygon segment is identifiable as distinct from adjacent link segments without additional parsing.

**Acceptance Scenarios**:

1. **Given** a matched path that includes polygon traversals with GPS points inside the polygon, **When** the output is inspected, **Then** each polygon segment is labeled with its polygon ID, the access point used to enter it, and the access point used to exit it; for segments where the trip starts or ends mid-polygon, the absent access point is represented as empty.
2. **Given** matched path output in the standard structured format, **When** parsed, **Then** link segments and polygon segments are distinguishable, and every polygon segment carries its entry and egress access point identifiers (or empty where absent).
3. **Given** a matched path that includes a through-routing polygon segment (no GPS observations inside the polygon), **When** the output is inspected, **Then** the polygon segment is present with both entry and egress access points recorded and is distinguishable from segments where GPS observations fell inside the polygon.
4. **Given** any polygon segment in the matched output, **When** the output is inspected, **Then** the segment carries a numeric distance value representing the Euclidean path length spent inside the polygon (entry-AP-to-first-inside-GPS + between-inside-GPS + last-inside-GPS-to-egress-AP for traversals, or entry-AP-to-egress-AP for through-routing).

---

### Edge Cases

- A GPS point that is simultaneously inside a polygon and within the search radius of a road link must be a candidate for both; the HMM must choose the best assignment. The polygon candidate's emission distance is zero; the link candidate's emission distance is computed normally.
- A GPS point on the exact boundary of a polygon must be treated as inside (emission distance zero). The boundary includes both outer rings and inner hole rings.
- A GPS point inside a polygon's hole (inner ring) must be treated as outside the polygon; its emission distance to the polygon equals the distance to the nearest ring (outer or hole).
- Two polygons that share an access point must be routable consecutively without requiring an intermediate link.
- A polygon with no adjacent road links but with access points shared with other polygons must still be a valid matched segment for GPS points contained within it.
- A GPS trajectory that begins inside a polygon must produce a valid matched start, recording the egress access point even when no entry access point exists (the trip starts mid-polygon).
- A GPS trajectory that ends inside a polygon must produce a valid matched end, recording the access point used to enter the polygon even when no egress access point is used.
- An access point feature in the shapefile whose geometry does not lie on the boundary of its associated polygon must be rejected with a descriptive error.
- A polygon with no features referencing it in the access point shapefile must be treated as unreachable and excluded from candidate generation, with a warning.
- An empty polygon layer (zero polygons) must not crash the algorithm and must fall back to link-only matching.
- A trajectory with zero or one GPS points must be skipped with a per-trajectory warning identifying the trajectory ID; remaining trajectories in the batch must continue to be matched.
- Polygons with invalid geometry (self-intersecting, zero-area) must be skipped with a warning, not silently accepted or cause a crash.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST accept a polygon layer (shapefile) and an access point shapefile as named input parameters alongside the road network, configurable via command-line arguments and XML configuration.
- **FR-002**: The system MUST load all valid polygons from the polygon layer and make them available as routing primitives before matching begins.
- **FR-003**: The system MUST treat each polygon as a first-class routing object with a unique identifier and an optional generalized cost value.
- **FR-004**: The system MUST load the access point shapefile, where each feature is a point with a node number and a polygon ID; features sharing the same node number represent a single access point shared between multiple polygons. The system MUST associate each access point with all polygons it references, and — when its node number matches a source or target node ID of any road link in the network — with all such links, before matching begins. Access points whose node number does not match any network node ID are polygon-shared only (no link attachment).
- **FR-005**: The system MUST validate the access point shapefile before matching and halt with a descriptive error if any of the following conditions are found:
  - An access point feature's geometry does not lie on the boundary of its declared polygon.
  - An access point feature references a polygon ID that does not exist in the polygon shapefile.
  - Features sharing the same node number have differing geometries (contradictory shared access point definitions).
- **FR-006**: The system MUST generate polygon candidates for GPS points that fall inside a polygon or within the global search radius of any of that polygon's access points. A GPS point MAY be a candidate for both a link and a polygon simultaneously; the HMM MUST evaluate both in the same candidate layer and select the best assignment. When computing the emission probability distance for a polygon candidate:
  - A GPS point **inside** the polygon, or on its **boundary**, MUST use distance **zero**.
  - A GPS point **outside** the polygon MUST use the **minimum straight-line distance from the GPS point to the polygon boundary**.
  - **Inside / outside / boundary** are defined per standard GIS semantics: a polygon's boundary includes its outer ring and any inner hole rings; a GPS point lying inside a hole is **outside** the polygon. For multipolygons (multiple disjoint outer rings), a point is inside if it lies inside any sub-polygon (and not in any of that sub-polygon's holes).
- **FR-007**: The system MUST allow routing transitions between all combinations: link→link, link→polygon, polygon→link, and polygon→polygon; transitions to or from a polygon MUST pass through one of that polygon's defined access points. The system MUST additionally support trajectories that begin or end mid-polygon: the first GPS point MAY be matched to a polygon with no entry access point, and the last GPS point MAY be matched to a polygon with no egress access point.
- **FR-008**: The system MUST apply polygon weight values when computing polygon transition costs according to the following cost model:
  - **Entry (outside → inside)**: cost = `polygon_weight × straight-line distance(entry access point, matched point inside polygon)`
  - **Egress (inside → outside)**: cost = `polygon_weight × straight-line distance(matched point inside polygon, egress access point)`
  - **Within-polygon (inside → inside, same polygon)**: the transition path distance MUST be set to the straight-line (Euclidean) distance between the two GPS points, mirroring the existing same-link override in WeightMatch (see `src/mm/weightmatch/weightmatch_algorithm.cpp:323-324`). This yields the most favorable transition probability via `calc_tp(eu_dist, eu_dist)` and effectively ignores routing cost within the polygon.
  - **Through-routing (outside → polygon → outside)**: cost = `polygon_weight × straight-line distance(entry access point, egress access point) × THROUGH_PENALTY_FACTOR`
- **FR-009**: THROUGH_PENALTY_FACTOR MUST be a configurable numeric parameter specifiable via command-line argument and XML configuration; its default value is `1.5`.
- **FR-010**: The system MUST produce a matched path output that identifies each segment as either a road link or a polygon; each polygon segment MUST record the access point used to enter and the access point used to exit (egress) it. Where a trip starts or ends mid-polygon (no link transition on that side), the absent access point MUST be represented as empty. Through-routing segments (no GPS observations inside the polygon) MUST be explicitly marked as through-routing so they are distinguishable from segments where GPS observations fell inside the polygon.
- **FR-016**: For each polygon segment in the matched path, the output MUST include a numeric **distance spent inside the polygon**, defined as the Euclidean path length contributed by that segment:
  - For a traversal segment with at least one GPS observation inside the polygon: the sum of `distance(entry_AP, first inside GPS)` + `Σ distance(consecutive inside GPS pairs)` + `distance(last inside GPS, egress_AP)`. Terms involving an absent access point (mid-polygon start or end) are omitted.
  - For a through-routing segment (no GPS observations inside): the straight-line distance from entry access point to egress access point.
  - For a trip entirely inside a single polygon: the sum of distances between consecutive GPS points only (no entry or egress contribution).
- **FR-017**: The system MUST support OpenMP-parallel matching across trajectories. The implementation MUST:
  - Treat all shared inputs (road network, polygon layer, access point layer, polygon-aware routing graph including the through-routing cost table) as read-only and safe to share across threads after construction.
  - Maintain per-thread instances of mutable routing state (Dijkstra state, indexed heap, HMM transition graph) so threads do not contend on these structures.
  - Serialize result output such that two trajectories matched on different threads never produce interleaved or corrupted output rows.
  - Produce results that are bit-for-bit identical regardless of thread count for the same input — parallelism MUST NOT alter matching outcomes.
- **FR-011**: The system MUST expose polymatch through the same command-line interface as weightmatch and stmatch.
- **FR-012**: The system MUST fall back to link-only matching when no polygon layer is provided or when the polygon layer contains zero valid polygons, with no change to output format.
- **FR-013**: The system MUST skip invalid polygon geometries (self-intersecting, zero-area) with a descriptive warning and continue processing valid polygons.
- **FR-014**: The system MUST exclude any polygon with no features referencing it in the access point shapefile from candidate generation, emitting a warning per excluded polygon.
- **FR-015**: The system MUST produce topologically valid matched paths: consecutive link segments MUST share a node; a link-to-polygon or polygon-to-link transition MUST occur through a shared access point attached to both; a polygon-to-polygon transition MUST occur through an access point node number that appears under both polygons in the access point shapefile.
- **FR-018**: The system MUST validate the polygon shapefile at load time and halt with a descriptive error if any of the following conditions are found:
  - A polygon feature has `id == 0` (ID 0 is reserved by the output format and cannot appear in input data).
  - Two or more polygon features share the same `id` (the error message MUST list the offending ID).
- **FR-019**: The system MUST handle trajectories with zero or one GPS points gracefully: emit a per-trajectory warning identifying the trajectory ID, skip the trajectory, and continue processing the remaining trajectories in the batch. No matched-path output row is produced for skipped trajectories.

### Key Entities

- **Polygon**: An area feature with a unique ID, a geometric boundary, and an optional generalized cost value; represents a traversable zone (e.g., parking lot, ferry terminal, pedestrian plaza).
- **Access Point**: A point on the boundary of one or more polygons that defines where routing may enter or exit a polygon. Represented in the access point shapefile as one feature per polygon relationship; features sharing the same node number identify a single physical access point that is shared across multiple polygons. An access point is attached to either one or more road links (enabling link↔polygon transitions) or to multiple polygons (enabling polygon↔polygon transitions).
- **PolyMatch Candidate**: A GPS point associated with either a road link or a polygon, carrying proximity score and snapping information used by the HMM.
- **Hybrid Matched Path**: An ordered sequence of matched segments, each identified as either a road link or a polygon traversal. Each polygon segment carries the entry and egress access point IDs (empty where the trip starts or ends mid-polygon), a flag indicating whether it is a through-routing segment, the numeric distance spent inside the polygon (FR-016), and forms a topologically connected route from the first to the last GPS observation.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A GPS trajectory that crosses a polygon area is matched with the polygon segment correctly included, verified against a ground-truth fixture with a known polygon traversal.
- **SC-002**: Link-only GPS trajectories produce results identical to weightmatch output (same matched edge sequences on the same network) — zero regressions.
- **SC-003**: Polymatch processes a 1,000-point GPS trajectory on a network with 100 polygons in under 10 seconds on a single core; with N OpenMP threads, throughput on a batch of independent trajectories scales near-linearly (≥ 0.7 × N speedup up to the available core count) compared to the single-core baseline.
- **SC-004**: All four existing test suites (algorithm_test, fmm_test, network_test, network_graph_test) pass without modification after polymatch is introduced.
- **SC-005**: Polygon segments in the matched output are identifiable without additional parsing in the command-line output format; each polygon segment includes entry and egress access point identifiers, represented as empty where the trip starts or ends mid-polygon.
- **SC-006**: The four polygon transition cost types (entry, egress, within-polygon, through-routing) each produce verifiable costs against known geometry fixtures — within-polygon transitions use the GPS Euclidean distance as path distance (matching the same-link override); entry/egress costs match the formula; through-routing cost changes proportionally when THROUGH_PENALTY_FACTOR is varied.
- **SC-007**: Emission probability distances for polygon candidates are verifiably zero for GPS points inside the polygon and equal to the minimum boundary distance for points outside, confirmed against fixtures with known geometries.
- **SC-008**: For each polygon segment in the matched output, the reported distance inside the polygon matches the formula in FR-016 within floating-point tolerance, verified against fixtures with hand-computed expected values for traversal, through-routing, mid-polygon start, mid-polygon end, and fully-inside cases.
- **SC-009**: Matching a batch of trajectories with N OpenMP threads produces bit-for-bit identical output rows to a single-threaded run (after sorting by trajectory ID), verified across N ∈ {1, 2, 4, 8}.

## Assumptions

- The polygon layer is provided as an ESRI shapefile with an ID field; other formats (GeoJSON, PostGIS) are out of scope for v1.
- Polygons and road links share the same coordinate reference system; no on-the-fly reprojection is performed. The system does not validate CRS consistency between input shapefiles — the user is responsible for ensuring matching CRSes.
- A GPS point inside a polygon is always a polygon candidate; a GPS point outside but within the global search radius of any of that polygon's access points is also a polygon candidate.
- Polygon-to-polygon routing is allowed only when the two polygons share an access point; geometric adjacency alone does not imply routability.
- Polymatch is an online algorithm requiring no precomputed lookup table, consistent with weightmatch and stmatch.
- No Python API bindings are provided for polymatch; it is a command-line-only algorithm.
- The default value of THROUGH_PENALTY_FACTOR is `1.5`, subject to refinement through testing.
- The polygon layer optional cost field uses the same column name convention as the road network cost field (`cost`), unless overridden in configuration.
