// TraceGenerator implementation. Produces the deterministic real-network
// trace batch for specs/002-real-network-validation. See research.md R3-R4
// for the determinism strategy and per-category constructors.

#include "network/trace_generator.hpp"

#include <algorithm>
#include <boost/geometry.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <locale>
#include <set>
#include <sstream>
#include <unordered_set>

#include "util/debug.hpp"

namespace FMM {
namespace NETWORK {

namespace bg = boost::geometry;

// ---------------- Construction ----------------

TraceGenerator::TraceGenerator(const Network& network,
                               const PolygonLayer& polygons,
                               const AccessPointLayer& access_points,
                               const ROUTING::LinkGraph& link_graph,
                               const ROUTING::PolyLinkGraph& poly_graph,
                               uint64_t seed)
    : network_(network),
      polygons_(polygons),
      access_points_(access_points),
      link_graph_(link_graph),
      poly_graph_(poly_graph),
      rng_(seed) {}

// ---------------- Helpers (file-local) ----------------

namespace {

double point_dist(const FMM::CORE::Point& a, const FMM::CORE::Point& b) {
  double dx = bg::get<0>(a) - bg::get<0>(b);
  double dy = bg::get<1>(a) - bg::get<1>(b);
  return std::sqrt(dx * dx + dy * dy);
}

FMM::CORE::Point edge_midpoint(const Edge& e) {
  const auto& g = e.geom;
  int n = g.get_num_points();
  if (n == 0) return FMM::CORE::Point(0, 0);
  int mid = n / 2;
  return FMM::CORE::Point(g.get_x(mid), g.get_y(mid));
}

void append_polyline_with_noise(FMM::CORE::LineString& out, const Edge& e,
                                 double noise_sigma, std::mt19937_64& rng,
                                 bool skip_first) {
  std::normal_distribution<double> nd(0.0, noise_sigma);
  const auto& g = e.geom;
  int n = g.get_num_points();
  for (int i = (skip_first ? 1 : 0); i < n; ++i) {
    double x = g.get_x(i) + nd(rng);
    double y = g.get_y(i) + nd(rng);
    out.add_point(x, y);
  }
}

// Build a polyline by walking `hops` link arcs starting from `start_edge`.
// Returns the EdgeIndex sequence the walk took.
std::vector<EdgeIndex> walk(const ROUTING::LinkGraph& g, EdgeIndex start,
                            int hops, std::mt19937_64& rng) {
  std::vector<EdgeIndex> path{start};
  for (int i = 0; i < hops; ++i) {
    const auto& nbrs = g.neighbors(path.back());
    if (nbrs.empty()) break;
    std::uniform_int_distribution<size_t> pick(0, nbrs.size() - 1);
    path.push_back(nbrs[pick(rng)].to);
  }
  return path;
}

// Build a LineString from a sequence of edges, sampling `n_points` evenly
// along the polyline with Gaussian noise. Returns geom with at least
// min(n_points, 2) points.
FMM::CORE::LineString polyline_from_edge_path(
    const std::vector<EdgeIndex>& path, const std::vector<Edge>& edges,
    int n_points, double noise_sigma, std::mt19937_64& rng) {
  // Flatten path into a single point sequence.
  std::vector<FMM::CORE::Point> pts;
  for (size_t k = 0; k < path.size(); ++k) {
    const Edge& e = edges[path[k]];
    int n = e.geom.get_num_points();
    int start = (k == 0) ? 0 : 1;
    for (int i = start; i < n; ++i) {
      pts.emplace_back(e.geom.get_x(i), e.geom.get_y(i));
    }
  }
  // Sample n_points evenly indexed.
  FMM::CORE::LineString out;
  if (pts.size() < 2) {
    for (auto& p : pts) {
      out.add_point(bg::get<0>(p), bg::get<1>(p));
    }
    return out;
  }
  std::normal_distribution<double> nd(0.0, noise_sigma);
  for (int i = 0; i < n_points; ++i) {
    double t = (n_points == 1) ? 0.5 : double(i) / double(n_points - 1);
    size_t idx = std::min<size_t>(pts.size() - 1,
                                  size_t(std::round(t * (pts.size() - 1))));
    double x = bg::get<0>(pts[idx]) + nd(rng);
    double y = bg::get<1>(pts[idx]) + nd(rng);
    out.add_point(x, y);
  }
  return out;
}

// Format a single GeneratedTrace as a CSV row using classic locale +
// setprecision(9). One row, no trailing newline.
std::string format_row(const GeneratedTrace& t) {
  std::ostringstream os;
  os.imbue(std::locale::classic());
  os << std::setprecision(9);
  os << t.id << ";";
  // WKT LINESTRING
  int n = t.geom.get_num_points();
  if (n == 0) {
    os << "LINESTRING EMPTY";
  } else {
    os << "LINESTRING(";
    for (int i = 0; i < n; ++i) {
      if (i > 0) os << ",";
      os << t.geom.get_x(i) << " " << t.geom.get_y(i);
    }
    os << ")";
  }
  os << ";" << to_label(t.category);
  return os.str();
}

}  // namespace

// ---------------- Random-walk helper ----------------

FMM::CORE::LineString TraceGenerator::random_walk_trace(int hops, int n_points,
                                                        double noise_sigma,
                                                        double search_radius) {
  const auto& edges = network_.get_edges();
  if (edges.empty()) return FMM::CORE::LineString{};

  // Try up to 20 random start edges to find one whose midpoint is >
  // search_radius from any polygon.
  std::uniform_int_distribution<size_t> pick_edge(0, edges.size() - 1);
  for (int attempt = 0; attempt < 20; ++attempt) {
    EdgeIndex start = static_cast<EdgeIndex>(pick_edge(rng_));
    auto mid = edge_midpoint(edges[start]);
    auto nearby = polygons_.polygons_within_radius(mid, search_radius);
    if (nearby.empty()) {
      auto path = walk(link_graph_, start, hops, rng_);
      auto geom = polyline_from_edge_path(path, edges, n_points, noise_sigma,
                                          rng_);
      if (geom.get_num_points() >= 2) return geom;
    }
  }
  // Fall back: any start, no polygon constraint.
  EdgeIndex start = static_cast<EdgeIndex>(pick_edge(rng_));
  auto path = walk(link_graph_, start, hops, rng_);
  return polyline_from_edge_path(path, edges, n_points, noise_sigma, rng_);
}

// ---------------- Category generators ----------------

std::vector<GeneratedTrace> TraceGenerator::generate_link_only(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  out.reserve(n_per_category);
  std::uniform_int_distribution<int> hop_pick(3, 8);
  std::uniform_int_distribution<int> pt_pick(5, 12);  // ≥ 5 to feed duplicate-points
  // Matcher uses radius=300m for candidates. To ensure no polygon is a
  // candidate at any GPS layer, every emitted point must be > matcher_radius
  // from any polygon. Use 350m to give a comfortable margin.
  constexpr double kKeepoutRadius = 350.0;
  int attempts = 0;
  const int kMaxAttempts = n_per_category * 30;
  while ((int)out.size() < n_per_category && attempts++ < kMaxAttempts) {
    int hops = hop_pick(rng_);
    int pts = pt_pick(rng_);
    auto geom = random_walk_trace(hops, pts, /*noise_sigma=*/15.0,
                                  /*search_radius=*/kKeepoutRadius);
    if (geom.get_num_points() < 2) continue;
    // Reject if any GPS point has a polygon within the keepout radius.
    bool ok = true;
    for (int i = 0; i < geom.get_num_points(); ++i) {
      FMM::CORE::Point p(geom.get_x(i), geom.get_y(i));
      if (!polygons_.polygons_within_radius(p, kKeepoutRadius).empty()) {
        ok = false;
        break;
      }
    }
    if (!ok) continue;
    out.push_back({1000 + (int)out.size(), geom, TraceCategory::LinkOnly});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_polygon_traversal(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  const auto& aps = access_points_.access_points();
  const auto& edges = network_.get_edges();

  // Find polygons with ≥ 2 link-attached APs.
  struct Cand { PolygonIndex p; std::vector<AccessPointIndex> link_aps; };
  std::vector<Cand> candidates;
  for (PolygonIndex p = 0; p < polygons_.size(); ++p) {
    Cand c{p, {}};
    for (auto idx : access_points_.aps_for_polygon(p)) {
      if (aps[idx].attached_node.has_value() &&
          !aps[idx].attached_edges.empty()) {
        c.link_aps.push_back(idx);
      }
    }
    if (c.link_aps.size() >= 2) candidates.push_back(std::move(c));
  }
  if (candidates.empty()) return out;

  std::uniform_int_distribution<size_t> pick_c(0, candidates.size() - 1);
  std::normal_distribution<double> nd(0.0, 10.0);
  for (int i = 0; i < n_per_category; ++i) {
    const Cand& c = candidates[pick_c(rng_)];
    std::uniform_int_distribution<size_t> pick_ap(0, c.link_aps.size() - 1);
    AccessPointIndex ap_in = c.link_aps[pick_ap(rng_)];
    AccessPointIndex ap_out;
    do { ap_out = c.link_aps[pick_ap(rng_)]; } while (ap_out == ap_in);
    const auto& a = aps[ap_in];
    const auto& b = aps[ap_out];
    const auto& poly = polygons_.polygons()[c.p];

    // Compute a polygon centroid as a midpoint inside.
    auto centroid_pt = bg::return_centroid<FMM::CORE::Point>(poly.geom);

    FMM::CORE::LineString geom;
    // Pre-AP point: pick the source edge endpoint that's NOT the AP node,
    // so the trace starts outside the polygon.
    EdgeIndex e_in = a.attached_edges[0];
    const Edge& ein = edges[e_in];
    NodeIndex other_node = (ein.source == *a.attached_node) ? ein.target : ein.source;
    auto other_pt = network_.get_node_geom_from_idx(other_node);
    geom.add_point(bg::get<0>(other_pt) + nd(rng_),
                   bg::get<1>(other_pt) + nd(rng_));
    geom.add_point(bg::get<0>(a.point) + nd(rng_),
                   bg::get<1>(a.point) + nd(rng_));
    geom.add_point(bg::get<0>(centroid_pt) + nd(rng_),
                   bg::get<1>(centroid_pt) + nd(rng_));
    geom.add_point(bg::get<0>(b.point) + nd(rng_),
                   bg::get<1>(b.point) + nd(rng_));
    EdgeIndex e_out = b.attached_edges[0];
    const Edge& eout = edges[e_out];
    NodeIndex other_node_out = (eout.source == *b.attached_node) ? eout.target : eout.source;
    auto other_pt_out = network_.get_node_geom_from_idx(other_node_out);
    geom.add_point(bg::get<0>(other_pt_out) + nd(rng_),
                   bg::get<1>(other_pt_out) + nd(rng_));

    out.push_back({1100 + i, geom, TraceCategory::PolygonTraversal});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_polygon_shared_ap(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  const auto& aps = access_points_.access_points();
  // Find APs that touch ≥ 2 polygons.
  std::vector<AccessPointIndex> shared;
  for (size_t i = 0; i < aps.size(); ++i) {
    if (aps[i].polygons.size() >= 2) shared.push_back(static_cast<AccessPointIndex>(i));
  }
  if (shared.empty()) return out;  // category not supported by this fixture

  std::uniform_int_distribution<size_t> pick(0, shared.size() - 1);
  std::normal_distribution<double> nd(0.0, 10.0);
  for (int i = 0; i < n_per_category; ++i) {
    const auto& ap = aps[shared[pick(rng_)]];
    PolygonIndex pa = ap.polygons[0];
    PolygonIndex pb = ap.polygons[1];
    auto ca = bg::return_centroid<FMM::CORE::Point>(polygons_.polygons()[pa].geom);
    auto cb = bg::return_centroid<FMM::CORE::Point>(polygons_.polygons()[pb].geom);
    FMM::CORE::LineString geom;
    geom.add_point(bg::get<0>(ca) + nd(rng_), bg::get<1>(ca) + nd(rng_));
    geom.add_point(bg::get<0>(ap.point) + nd(rng_),
                   bg::get<1>(ap.point) + nd(rng_));
    geom.add_point(bg::get<0>(cb) + nd(rng_), bg::get<1>(cb) + nd(rng_));
    out.push_back({1200 + i, geom, TraceCategory::PolygonSharedAp});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_mid_polygon_start(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  const auto& aps = access_points_.access_points();
  // Need polygons with at least one link-attached AP to route OUT.
  std::vector<PolygonIndex> candidates;
  for (PolygonIndex p = 0; p < polygons_.size(); ++p) {
    for (auto idx : access_points_.aps_for_polygon(p)) {
      if (aps[idx].attached_node.has_value() &&
          !aps[idx].attached_edges.empty()) {
        candidates.push_back(p);
        break;
      }
    }
  }
  if (candidates.empty()) return out;

  std::uniform_int_distribution<size_t> pick(0, candidates.size() - 1);
  std::normal_distribution<double> nd(0.0, 8.0);
  for (int i = 0; i < n_per_category; ++i) {
    PolygonIndex p = candidates[pick(rng_)];
    auto centroid = bg::return_centroid<FMM::CORE::Point>(polygons_.polygons()[p].geom);
    // Find one link-attached AP.
    AccessPointIndex ap_idx = 0;
    for (auto idx : access_points_.aps_for_polygon(p)) {
      if (aps[idx].attached_node.has_value() &&
          !aps[idx].attached_edges.empty()) {
        ap_idx = idx;
        break;
      }
    }
    const auto& ap = aps[ap_idx];
    FMM::CORE::LineString geom;
    geom.add_point(bg::get<0>(centroid) + nd(rng_),
                   bg::get<1>(centroid) + nd(rng_));
    geom.add_point(bg::get<0>(ap.point) + nd(rng_),
                   bg::get<1>(ap.point) + nd(rng_));
    // exit toward an adjacent network node
    EdgeIndex e_idx = ap.attached_edges[0];
    const Edge& e = network_.get_edges()[e_idx];
    NodeIndex other = (e.source == *ap.attached_node) ? e.target : e.source;
    auto op = network_.get_node_geom_from_idx(other);
    geom.add_point(bg::get<0>(op) + nd(rng_), bg::get<1>(op) + nd(rng_));
    out.push_back({1300 + i, geom, TraceCategory::MidPolygonStart});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_mid_polygon_end(
    int n_per_category) {
  // Symmetric of mid-polygon-start.
  auto starts = generate_mid_polygon_start(n_per_category);
  std::vector<GeneratedTrace> out;
  for (size_t i = 0; i < starts.size(); ++i) {
    GeneratedTrace t = std::move(starts[i]);
    // Reverse the polyline.
    FMM::CORE::LineString rev;
    int n = t.geom.get_num_points();
    for (int j = n - 1; j >= 0; --j) {
      rev.add_point(t.geom.get_x(j), t.geom.get_y(j));
    }
    t.geom = rev;
    t.id = 1400 + static_cast<int>(i);
    t.category = TraceCategory::MidPolygonEnd;
    out.push_back(std::move(t));
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_fully_inside(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  if (polygons_.size() == 0) return out;
  std::uniform_int_distribution<PolygonIndex> pick(0, polygons_.size() - 1);
  for (int i = 0; i < n_per_category; ++i) {
    PolygonIndex p = pick(rng_);
    const auto& poly_geom = polygons_.polygons()[p].geom;
    auto bbox = bg::return_envelope<bg::model::box<FMM::CORE::Point>>(poly_geom);
    double xmin = bg::get<0>(bbox.min_corner());
    double ymin = bg::get<1>(bbox.min_corner());
    double xmax = bg::get<0>(bbox.max_corner());
    double ymax = bg::get<1>(bbox.max_corner());
    std::uniform_real_distribution<double> ux(xmin, xmax);
    std::uniform_real_distribution<double> uy(ymin, ymax);
    FMM::CORE::LineString geom;
    int tries = 0;
    while (geom.get_num_points() < 4 && tries < 200) {
      FMM::CORE::Point cand(ux(rng_), uy(rng_));
      if (bg::covered_by(cand, poly_geom)) {
        geom.add_point(bg::get<0>(cand), bg::get<1>(cand));
      }
      ++tries;
    }
    if (geom.get_num_points() < 2) continue;
    out.push_back({1500 + i, geom, TraceCategory::FullyInside});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_through_routing(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  const auto& aps = access_points_.access_points();
  // Polygons with ≥ 2 link-attached APs (same set as polygon-traversal but
  // we sample GPS points OUTSIDE the polygon to make it a pure shortcut).
  struct Cand { PolygonIndex p; std::vector<AccessPointIndex> link_aps; };
  std::vector<Cand> candidates;
  for (PolygonIndex p = 0; p < polygons_.size(); ++p) {
    Cand c{p, {}};
    for (auto idx : access_points_.aps_for_polygon(p)) {
      if (aps[idx].attached_node.has_value() &&
          !aps[idx].attached_edges.empty()) {
        c.link_aps.push_back(idx);
      }
    }
    if (c.link_aps.size() >= 2) candidates.push_back(std::move(c));
  }
  if (candidates.empty()) return out;

  std::uniform_int_distribution<size_t> pick(0, candidates.size() - 1);
  std::normal_distribution<double> nd(0.0, 5.0);
  for (int i = 0; i < n_per_category; ++i) {
    const Cand& c = candidates[pick(rng_)];
    std::uniform_int_distribution<size_t> pick_ap(0, c.link_aps.size() - 1);
    AccessPointIndex ap_in = c.link_aps[pick_ap(rng_)];
    AccessPointIndex ap_out;
    do { ap_out = c.link_aps[pick_ap(rng_)]; } while (ap_out == ap_in);
    const auto& a = aps[ap_in];
    const auto& b = aps[ap_out];
    const auto& poly_geom = polygons_.polygons()[c.p].geom;
    auto bbox = bg::return_envelope<bg::model::box<FMM::CORE::Point>>(poly_geom);
    double xmin = bg::get<0>(bbox.min_corner());
    double ymin = bg::get<1>(bbox.min_corner());
    double xmax = bg::get<0>(bbox.max_corner());
    double ymax = bg::get<1>(bbox.max_corner());
    double pad = std::max(xmax - xmin, ymax - ymin);  // ~polygon diameter
    EdgeIndex e_in = a.attached_edges[0];
    const Edge& ein = network_.get_edges()[e_in];
    NodeIndex on_a = (ein.source == *a.attached_node) ? ein.target : ein.source;
    auto pa = network_.get_node_geom_from_idx(on_a);
    EdgeIndex e_out = b.attached_edges[0];
    const Edge& eout = network_.get_edges()[e_out];
    NodeIndex on_b = (eout.source == *b.attached_node) ? eout.target : eout.source;
    auto pb = network_.get_node_geom_from_idx(on_b);
    FMM::CORE::LineString geom;
    geom.add_point(bg::get<0>(pa) + nd(rng_), bg::get<1>(pa) + nd(rng_));
    geom.add_point(bg::get<0>(pb) + nd(rng_), bg::get<1>(pb) + nd(rng_));
    (void)pad;  // documented intent — endpoint nodes are outside the polygon
    out.push_back({1600 + i, geom, TraceCategory::ThroughRouting});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_off_network_noise(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  out.reserve(n_per_category);
  std::normal_distribution<double> nd(0.0, 1500.0);  // ~ 5 × radius (300m)
  // Generate enough base traces that we can produce n_per_category off-network
  // ones (each needs ≥ 3 points). Request 2x to be safe.
  auto base = generate_link_only(n_per_category * 2);
  for (size_t i = 0; i < base.size() && (int)out.size() < n_per_category; ++i) {
    GeneratedTrace t = std::move(base[i]);
    if (t.geom.get_num_points() < 3) continue;
    int mid = t.geom.get_num_points() / 2;
    FMM::CORE::LineString patched;
    for (int j = 0; j < t.geom.get_num_points(); ++j) {
      patched.add_point(t.geom.get_x(j), t.geom.get_y(j));
      if (j == mid) {
        patched.add_point(t.geom.get_x(j) + nd(rng_),
                          t.geom.get_y(j) + nd(rng_));
      }
    }
    t.geom = patched;
    t.id = 1700 + (int)out.size();
    t.category = TraceCategory::OffNetworkNoise;
    out.push_back(std::move(t));
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_short_trip(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  out.reserve(n_per_category);
  std::uniform_int_distribution<int> pt_pick(2, 3);
  int attempts = 0;
  const int kMaxAttempts = n_per_category * 6;
  while ((int)out.size() < n_per_category && attempts++ < kMaxAttempts) {
    int pts = pt_pick(rng_);
    auto geom = random_walk_trace(/*hops=*/1, pts, /*noise=*/5.0,
                                  /*search_radius=*/100.0);
    if (geom.get_num_points() < 2) continue;
    out.push_back({1800 + (int)out.size(), geom, TraceCategory::ShortTrip});
  }
  return out;
}

std::vector<GeneratedTrace> TraceGenerator::generate_duplicate_points(
    int n_per_category) {
  std::vector<GeneratedTrace> out;
  out.reserve(n_per_category);
  // Request 2x base traces so the ≥ 5 point filter doesn't undershoot.
  auto base = generate_link_only(n_per_category * 2);
  for (size_t i = 0; i < base.size() && (int)out.size() < n_per_category; ++i) {
    GeneratedTrace t = std::move(base[i]);
    if (t.geom.get_num_points() < 5) continue;
    int mid = t.geom.get_num_points() / 2;
    FMM::CORE::LineString patched;
    for (int j = 0; j < t.geom.get_num_points(); ++j) {
      if (j == mid) {
        patched.add_point(t.geom.get_x(j - 1), t.geom.get_y(j - 1));
      } else {
        patched.add_point(t.geom.get_x(j), t.geom.get_y(j));
      }
    }
    t.geom = patched;
    t.id = 1900 + (int)out.size();
    t.category = TraceCategory::DuplicatePoints;
    out.push_back(std::move(t));
  }
  return out;
}

// ---------------- CSV writer ----------------

bool TraceGenerator::write_csv(const std::string& out_path,
                               const std::vector<GeneratedTrace>& traces_in) const {
  // Sort by id (FR-005 + contract requires id-ascending order).
  std::vector<GeneratedTrace> traces = traces_in;
  std::sort(traces.begin(), traces.end(),
            [](const GeneratedTrace& a, const GeneratedTrace& b) {
              return a.id < b.id;
            });

  std::ofstream out(out_path, std::ios::binary);
  if (!out.is_open()) return false;
  out.imbue(std::locale::classic());
  out << "id;geom;category\n";
  for (const auto& t : traces) {
    out << format_row(t) << "\n";
  }
  return out.good();
}

} // namespace NETWORK
} // namespace FMM
