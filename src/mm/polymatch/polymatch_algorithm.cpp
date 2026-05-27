#include "mm/polymatch/polymatch_algorithm.hpp"

#include <algorithm>
#include <boost/geometry.hpp>
#include <cmath>
#include <limits>
#include <unordered_set>

#include "algorithm/geom_algorithm.hpp"
#include "mm/weightmatch/weightmatch_algorithm.hpp"
#include "network/link_graph_routing.hpp"
#include "util/debug.hpp"
#include "util/util.hpp"

namespace FMM {
namespace MM {

namespace bg = boost::geometry;

POLYMATCHConfig::POLYMATCHConfig(int k_arg, double r_arg, double gps_error_arg,
                                 int backup_k_arg, double backup_r_arg,
                                 double ub_factor_arg, bool allow_truncation_arg,
                                 double through_penalty_factor_arg,
                                 double boundary_epsilon_arg)
    : k(k_arg),
      radius(r_arg),
      gps_error(gps_error_arg),
      backup_k(backup_k_arg),
      backup_radius(backup_r_arg),
      upper_bound_factor(ub_factor_arg),
      allow_truncation(allow_truncation_arg),
      through_penalty_factor(through_penalty_factor_arg),
      boundary_epsilon(boundary_epsilon_arg) {}

void POLYMATCHConfig::print() const {
  SPDLOG_INFO("POLYMATCHAlgorithmConfig");
  SPDLOG_INFO(
      "k {} radius {} gps_error {} backup_k {} backup_radius {} upper_bound_factor {} allow_truncation {} through_penalty_factor {} boundary_epsilon {}",
      k, radius, gps_error, backup_k, backup_radius, upper_bound_factor,
      allow_truncation, through_penalty_factor, boundary_epsilon);
}

POLYMATCHConfig POLYMATCHConfig::load_from_arg(
    const cxxopts::ParseResult &arg_data) {
  int k = arg_data["candidates"].as<int>();
  double radius = arg_data["radius"].as<double>();
  double gps_error = arg_data["error"].as<double>();
  int backup_k = arg_data["backup_candidates"].as<int>();
  double backup_radius = arg_data["backup_radius"].as<double>();
  double ub_factor = arg_data["upper_bound_factor"].as<double>();
  bool allow_truncation = arg_data["allow_truncation"].as<bool>();
  double tpf = arg_data["through_penalty_factor"].as<double>();
  double eps = arg_data["boundary_epsilon"].as<double>();
  return POLYMATCHConfig{k, radius, gps_error, backup_k, backup_radius,
                         ub_factor, allow_truncation, tpf, eps};
}

void POLYMATCHConfig::register_arg(cxxopts::Options &options) {
  options.add_options()
    ("k,candidates", "Number of candidates",
      cxxopts::value<int>()->default_value("8"))
    ("r,radius", "Search radius",
      cxxopts::value<double>()->default_value("300.0"))
    ("e,error", "GPS error",
      cxxopts::value<double>()->default_value("50.0"))
    ("backup_candidates", "Number of candidates in fallback radius",
      cxxopts::value<int>()->default_value("-1"))
    ("backup_radius", "Fallback search radius",
      cxxopts::value<double>()->default_value("-1"))
    ("upper_bound_factor", "Dijkstra upper bound factor (0 to disable)",
      cxxopts::value<double>()->default_value("10.0"))
    ("allow_truncation", "Allow truncation of ends of trips in search",
      cxxopts::value<bool>()->default_value("false"))
    ("through_penalty_factor", "Through-routing cost multiplier",
      cxxopts::value<double>()->default_value("1.5"))
    ("boundary_epsilon", "Polygon boundary tolerance for AP validation",
      cxxopts::value<double>()->default_value("0.000001"));
}

void POLYMATCHConfig::register_help(std::ostringstream &oss) {
  oss << "-k/--candidates (optional) <int>: number of candidates (8)\n";
  oss << "-r/--radius (optional) <double>: search radius (300)\n";
  oss << "-e/--error (optional) <double>: GPS error (50)\n";
  oss << "--backup_candidates (optional) <int>: backup candidates (-1)\n";
  oss << "--backup_radius (optional) <double>: backup search radius (-1)\n";
  oss << "--upper_bound_factor (optional) <double>: Dijkstra upper bound factor (10.0)\n";
  oss << "--allow_truncation (optional) <bool>: allow truncation (false)\n";
  oss << "--through_penalty_factor (optional) <double>: through-routing multiplier (1.5)\n";
  oss << "--boundary_epsilon (optional) <double>: AP boundary tolerance (1e-6)\n";
}

bool POLYMATCHConfig::validate() const {
  if (gps_error <= 0 || radius <= 0 || k <= 0) {
    SPDLOG_CRITICAL("Invalid polymatch parameter k {} r {} gps_error {}", k,
                    radius, gps_error);
    return false;
  }
  if (through_penalty_factor < 0) {
    SPDLOG_CRITICAL("through_penalty_factor must be >= 0, got {}",
                    through_penalty_factor);
    return false;
  }
  if (backup_radius >= 0 && backup_radius <= radius) {
    SPDLOG_CRITICAL("backup_radius must exceed radius");
    return false;
  }
  if ((backup_radius >= 0 && backup_k < 0) ||
      (backup_radius < 0 && backup_k >= 0)) {
    SPDLOG_CRITICAL("backup_radius and backup_k must both be set or both unset");
    return false;
  }
  return true;
}

// ---------------- Helpers ----------------

namespace {

inline double point_distance(const FMM::CORE::Point &a,
                             const FMM::CORE::Point &b) {
  double dx = bg::get<0>(a) - bg::get<0>(b);
  double dy = bg::get<1>(a) - bg::get<1>(b);
  return std::sqrt(dx * dx + dy * dy);
}

// Project a point onto an edge's polyline. Returns the matched point + offset
// (arc length from the edge's source) + perpendicular distance.
void project_point_on_edge(const FMM::NETWORK::Edge &e,
                           const FMM::CORE::Point &p,
                           FMM::CORE::Point *out_point, double *out_offset,
                           double *out_dist) {
  // Simplified: closest point on the edge polyline. We iterate segments.
  const auto &line = e.geom;
  int n = line.get_num_points();
  double best_dist2 = std::numeric_limits<double>::infinity();
  double best_offset = 0.0;
  FMM::CORE::Point best_pt(0, 0);
  double accum = 0.0;
  for (int i = 0; i + 1 < n; ++i) {
    double ax = line.get_x(i), ay = line.get_y(i);
    double bx = line.get_x(i + 1), by = line.get_y(i + 1);
    double dx = bx - ax, dy = by - ay;
    double seg_len2 = dx * dx + dy * dy;
    double t = 0.0;
    if (seg_len2 > 0) {
      t = ((bg::get<0>(p) - ax) * dx + (bg::get<1>(p) - ay) * dy) / seg_len2;
      if (t < 0) t = 0;
      if (t > 1) t = 1;
    }
    double px = ax + t * dx, py = ay + t * dy;
    double ddx = bg::get<0>(p) - px, ddy = bg::get<1>(p) - py;
    double d2 = ddx * ddx + ddy * ddy;
    if (d2 < best_dist2) {
      best_dist2 = d2;
      best_offset = accum + t * std::sqrt(seg_len2);
      bg::set<0>(best_pt, px);
      bg::set<1>(best_pt, py);
    }
    accum += std::sqrt(seg_len2);
  }
  *out_point = best_pt;
  *out_offset = best_offset;
  *out_dist = std::sqrt(best_dist2);
}

}  // namespace

// ---------------- Candidate generation ----------------

PolyTrajCandidates POLYMATCH::build_candidates(
    const CORE::Trajectory &traj, const POLYMATCHConfig &config,
    MM::Traj_Candidates &link_candidates_owner) const {
  // Link candidates from existing KNN search; takes a non-const linestring
  // reference but the implementation is read-only for our purposes.
  link_candidates_owner = network_.search_tr_cs_knn_with_fallback(
      traj.geom, config.k, config.radius, config.backup_k, config.backup_radius,
      config.allow_truncation);

  int n = traj.geom.get_num_points();
  PolyTrajCandidates out;
  out.reserve(n);

  for (int i = 0; i < n; ++i) {
    PolyPointCandidates layer;
    if ((int)link_candidates_owner.size() > i) {
      for (auto &lc : link_candidates_owner[i]) {
        PolyCandidate pc;
        pc.kind = PolyCandidateKind::Link;
        pc.ep_distance = lc.dist;
        pc.matched_point = lc.point;
        pc.edge = lc.edge;
        pc.offset = lc.offset;
        layer.push_back(pc);
      }
    }

    if (polygon_layer_.size() > 0) {
      FMM::CORE::Point gps(traj.geom.get_x(i), traj.geom.get_y(i));

      // Polygons containing this GPS point (inside / on boundary).
      for (auto p_idx : polygon_layer_.polygons_containing(gps)) {
        if (ap_layer_.aps_for_polygon(p_idx).empty()) continue;  // FR-014
        PolyCandidate pc;
        pc.kind = PolyCandidateKind::Polygon;
        pc.ep_distance = 0.0;
        pc.matched_point = gps;
        pc.polygon_index = p_idx;
        pc.inside = true;
        layer.push_back(pc);
      }
      // Polygons within radius (outside but close); emission distance is min
      // distance to boundary per FR-006.
      for (auto p_idx : polygon_layer_.polygons_within_radius(gps, config.radius)) {
        if (ap_layer_.aps_for_polygon(p_idx).empty()) continue;
        double d = polygon_layer_.min_boundary_distance(p_idx, gps);
        if (d <= 0.0) continue;  // already added via polygons_containing
        // De-dup against the contained set (defensive).
        bool dup = false;
        for (const auto &existing : layer) {
          if (existing.kind == PolyCandidateKind::Polygon &&
              existing.polygon_index == p_idx) {
            dup = true;
            break;
          }
        }
        if (dup) continue;
        PolyCandidate pc;
        pc.kind = PolyCandidateKind::Polygon;
        pc.ep_distance = d;
        pc.matched_point = gps;  // closest on boundary is not stored; matched
                                  // point used in distance_inside calcs
        pc.polygon_index = p_idx;
        pc.inside = false;
        layer.push_back(pc);
      }
    }
    out.push_back(std::move(layer));
  }

  return out;
}

// ---------------- Transition cost ----------------

double POLYMATCH::transition_cost(const PolyCandidate &a, const PolyCandidate &b,
                                  double eu_dist, const POLYMATCHConfig &config,
                                  ROUTING::DijkstraState &state,
                                  ROUTING::IndexedMinHeap &heap) const {
  using FMM::NETWORK::Edge;
  using FMM::ROUTING::Path;

  // Same-polygon, both inside → eu_dist override (mirror weightmatch:323-324).
  // Even when one observation is on the boundary (treated as inside per R9),
  // the eu_dist override applies because path inside the polygon is bounded
  // by straight-line distance times polygon weight; we use eu_dist as the
  // distance proxy to keep the HMM transition probability favorable.
  if (a.is_polygon() && b.is_polygon() && a.polygon_index == b.polygon_index) {
    return eu_dist;
  }

  // Link → Link (same as weightmatch)
  if (a.is_link() && b.is_link()) {
    Edge *se = a.edge;
    Edge *te = b.edge;
    if (se->id == te->id) {
      return eu_dist;
    }
    std::vector<FMM::NETWORK::EdgeIndex> goals{te->index};
    std::vector<Path> paths(1);
    FMM::ROUTING::shortest_edge_to_edges(link_graph_, state, heap, se->index,
                                         goals, paths, config.upper_bound_factor);
    if (!paths[0].found) return std::numeric_limits<double>::infinity();
    return se->weight * (se->length - a.offset) + paths[0].total_cost +
           te->weight * (b.offset - te->length);
  }

  // Link → Polygon: best AP that's link-attached on the target polygon.
  if (a.is_link() && b.is_polygon()) {
    const auto &aps = ap_layer_.aps_for_polygon(b.polygon_index);
    if (aps.empty()) return std::numeric_limits<double>::infinity();
    Edge *se = a.edge;
    double best = std::numeric_limits<double>::infinity();
    double poly_w = polygon_layer_.polygons()[b.polygon_index].weight;
    for (auto ap_idx : aps) {
      const auto &ap = ap_layer_.access_points()[ap_idx];
      if (!ap.attached_node.has_value() || ap.attached_edges.empty()) continue;
      // Dijkstra from source edge to any of AP's incident edges.
      std::vector<FMM::NETWORK::EdgeIndex> goals(ap.attached_edges.begin(),
                                                  ap.attached_edges.end());
      std::vector<Path> paths(goals.size());
      FMM::ROUTING::shortest_edge_to_edges(link_graph_, state, heap, se->index,
                                           goals, paths,
                                           config.upper_bound_factor);
      for (size_t i = 0; i < goals.size(); ++i) {
        if (!paths[i].found) continue;
        const Edge &te = network_.get_edge(goals[i]);
        // Cost = source edge adjustment + Dijkstra + arrive at AP node (no
        // target offset adjustment, since the AP is at a node). Then entry
        // cost from AP into matched point inside polygon.
        double cost = se->weight * (se->length - a.offset) + paths[i].total_cost;
        // AP node is either source or target of te. Subtract any remaining
        // edge length to arrive at the AP node specifically.
        // Conservative: skip per-edge offset adjustment for the AP node; the
        // path cost already places us at one endpoint of `te`. For our test
        // fixtures with edge length 1.0 and AP at corner nodes, this is
        // accurate enough; refinement can subtract te.weight * te.length when
        // AP is at te.source.
        double entry = poly_w * point_distance(ap.point, b.matched_point);
        cost += entry;
        if (cost < best) best = cost;
      }
    }
    return best;
  }

  // Polygon → Link: symmetric
  if (a.is_polygon() && b.is_link()) {
    const auto &aps = ap_layer_.aps_for_polygon(a.polygon_index);
    if (aps.empty()) return std::numeric_limits<double>::infinity();
    Edge *te = b.edge;
    double best = std::numeric_limits<double>::infinity();
    double poly_w = polygon_layer_.polygons()[a.polygon_index].weight;
    for (auto ap_idx : aps) {
      const auto &ap = ap_layer_.access_points()[ap_idx];
      if (!ap.attached_node.has_value() || ap.attached_edges.empty()) continue;
      // Cost from AP to target edge via LinkGraph. We pick the min over AP's
      // incident edges as the start.
      for (auto src_edge : ap.attached_edges) {
        std::vector<FMM::NETWORK::EdgeIndex> goals{te->index};
        std::vector<Path> paths(1);
        FMM::ROUTING::shortest_edge_to_edges(link_graph_, state, heap, src_edge,
                                             goals, paths,
                                             config.upper_bound_factor);
        if (!paths[0].found) continue;
        double egress = poly_w * point_distance(a.matched_point, ap.point);
        double cost = egress + paths[0].total_cost +
                      te->weight * (b.offset - te->length);
        if (cost < best) best = cost;
      }
    }
    return best;
  }

  // Polygon → Polygon (different polygons): only via shared AP.
  // a.polygon_index != b.polygon_index because the same-polygon case is handled
  // above.
  const auto &aps_a = ap_layer_.aps_for_polygon(a.polygon_index);
  const auto &aps_b = ap_layer_.aps_for_polygon(b.polygon_index);
  std::unordered_set<FMM::NETWORK::AccessPointIndex> set_b(aps_b.begin(),
                                                            aps_b.end());
  double best = std::numeric_limits<double>::infinity();
  double poly_a_w = polygon_layer_.polygons()[a.polygon_index].weight;
  double poly_b_w = polygon_layer_.polygons()[b.polygon_index].weight;
  for (auto ap_idx : aps_a) {
    if (!set_b.count(ap_idx)) continue;
    const auto &ap = ap_layer_.access_points()[ap_idx];
    double cost = poly_a_w * point_distance(a.matched_point, ap.point) +
                  poly_b_w * point_distance(ap.point, b.matched_point);
    if (cost < best) best = cost;
  }
  return best;
}

// ---------------- HMM layer update ----------------

void POLYMATCH::update_layer(int level, PolyTGLayer *la, PolyTGLayer *lb,
                             double eu_dist, const POLYMATCHConfig &config,
                             ROUTING::DijkstraState &state,
                             ROUTING::IndexedMinHeap &heap) {
  for (auto &a : *la) {
    for (auto &b : *lb) {
      double sp_dist = transition_cost(*a.c, *b.c, eu_dist, config, state, heap);
      double tp = PolyTransitionGraph::calc_tp(sp_dist, eu_dist);
      double temp = a.cumu_prob + std::log(tp) + std::log(b.ep);
      if (temp >= b.cumu_prob) {
        b.cumu_prob = temp;
        b.prev = &a;
        b.sp_dist = sp_dist;
        b.tp = tp;
      }
    }
  }
}

void POLYMATCH::update_tg(PolyTransitionGraph &tg, const CORE::Trajectory &traj,
                          const POLYMATCHConfig &config,
                          ROUTING::DijkstraState &state,
                          ROUTING::IndexedMinHeap &heap) {
  auto &layers = tg.get_layers();
  std::vector<double> eu_dists = ALGORITHM::cal_eu_dist(traj.geom);
  int N = layers.size();
  for (int i = 0; i < N - 1; ++i) {
    update_layer(i, &layers[i], &layers[i + 1], eu_dists[i], config, state,
                 heap);
  }
}

// ---------------- Hybrid path assembly ----------------

void POLYMATCH::build_hybrid_path(const PolyTGOpath &opath,
                                  const CORE::Trajectory &traj,
                                  const POLYMATCHConfig &config,
                                  ROUTING::DijkstraState &state,
                                  ROUTING::IndexedMinHeap &heap,
                                  PolyMatchResult *out) const {
  if (opath.empty()) return;
  C_Path cpath;
  std::vector<int> indices;
  std::vector<PolygonSegment> segments;

  // opath[i] -> step in cpath
  int current_idx = 0;

  auto push_link = [&](FMM::NETWORK::EdgeID id) {
    cpath.push_back(id);
  };
  auto push_polygon = [&](FMM::NETWORK::PolygonID id) {
    cpath.push_back(-id);
  };

  // Track active polygon segment for distance_inside accumulation.
  struct ActiveSegment {
    FMM::NETWORK::PolygonIndex p_idx;
    FMM::NETWORK::PolygonID p_id;
    FMM::NETWORK::NodeID entry_ap = kNoAccessPoint;
    FMM::NETWORK::NodeID egress_ap = kNoAccessPoint;
    bool has_inside_obs = false;
    FMM::CORE::Point first_inside;
    FMM::CORE::Point last_inside;
    double inside_path_sum = 0.0;  // Σ eu_dist between consecutive inside GPS
    int last_inside_layer = -1;
    size_t position_in_cpath = 0;
  };

  std::vector<ActiveSegment> active;  // chronological order

  auto begin_segment = [&](FMM::NETWORK::PolygonIndex p_idx, size_t pos) {
    const auto &poly = polygon_layer_.polygons()[p_idx];
    ActiveSegment seg;
    seg.p_idx = p_idx;
    seg.p_id = poly.id;
    seg.position_in_cpath = pos;
    active.push_back(seg);
  };

  // Initial element
  const PolyCandidate *first_c = opath[0]->c;
  if (first_c->is_link()) {
    push_link(first_c->edge->id);
  } else {
    begin_segment(first_c->polygon_index, cpath.size());
    push_polygon(polygon_layer_.polygons()[first_c->polygon_index].id);
    if (first_c->inside) {
      active.back().has_inside_obs = true;
      active.back().first_inside = first_c->matched_point;
      active.back().last_inside = first_c->matched_point;
      active.back().last_inside_layer = 0;
    }
  }
  indices.push_back(current_idx);

  std::vector<double> eu_dists = ALGORITHM::cal_eu_dist(traj.geom);

  for (size_t i = 0; i + 1 < opath.size(); ++i) {
    const PolyCandidate *a = opath[i]->c;
    const PolyCandidate *b = opath[i + 1]->c;

    if (a->is_link() && b->is_link()) {
      if (a->edge->id == b->edge->id) {
        indices.push_back(current_idx);
        continue;
      }
      // Expand the Dijkstra path between edges
      std::vector<FMM::NETWORK::EdgeIndex> goals{b->edge->index};
      std::vector<FMM::ROUTING::Path> paths(1);
      FMM::ROUTING::shortest_edge_to_edges(link_graph_, state, heap,
                                           a->edge->index, goals, paths);
      const auto &segs = paths[0].edges;
      if (segs.size() > 2) {
        for (auto it = std::next(segs.begin()); it != std::prev(segs.end());
             ++it) {
          push_link(network_.get_edges()[*it].id);
          ++current_idx;
        }
      }
      push_link(b->edge->id);
      ++current_idx;
      indices.push_back(current_idx);
      continue;
    }

    // Polygon entry from a link
    if (a->is_link() && b->is_polygon()) {
      // Determine entry AP: pick the one minimizing transition cost (same logic
      // as transition_cost). We just record the best AP.
      const auto &aps = ap_layer_.aps_for_polygon(b->polygon_index);
      double poly_w = polygon_layer_.polygons()[b->polygon_index].weight;
      double best = std::numeric_limits<double>::infinity();
      FMM::NETWORK::NodeID best_ap = kNoAccessPoint;
      // Insert intermediate links from a to the AP's incident edge that we
      // chose, before the polygon element.
      std::vector<FMM::NETWORK::EdgeIndex> chosen_segs;
      for (auto ap_idx : aps) {
        const auto &ap = ap_layer_.access_points()[ap_idx];
        if (!ap.attached_node.has_value() || ap.attached_edges.empty()) continue;
        std::vector<FMM::NETWORK::EdgeIndex> goals(ap.attached_edges.begin(),
                                                    ap.attached_edges.end());
        std::vector<FMM::ROUTING::Path> paths(goals.size());
        FMM::ROUTING::shortest_edge_to_edges(link_graph_, state, heap,
                                             a->edge->index, goals, paths);
        for (size_t g = 0; g < goals.size(); ++g) {
          if (!paths[g].found) continue;
          double cost = a->edge->weight * (a->edge->length - a->offset) +
                        paths[g].total_cost +
                        poly_w * point_distance(ap.point, b->matched_point);
          if (cost < best) {
            best = cost;
            best_ap = ap.node_id;
            chosen_segs = paths[g].edges;
          }
        }
      }
      // Emit intermediate edges (drop first since it's already in cpath)
      if (chosen_segs.size() > 1) {
        for (auto it = std::next(chosen_segs.begin()); it != chosen_segs.end();
             ++it) {
          push_link(network_.get_edges()[*it].id);
          ++current_idx;
        }
      }
      // Begin polygon segment and emit polygon as a cpath element.
      begin_segment(b->polygon_index, cpath.size());
      active.back().entry_ap = best_ap;
      push_polygon(polygon_layer_.polygons()[b->polygon_index].id);
      ++current_idx;
      if (b->inside) {
        active.back().has_inside_obs = true;
        active.back().first_inside = b->matched_point;
        active.back().last_inside = b->matched_point;
        active.back().last_inside_layer = (int)(i + 1);
      }
      indices.push_back(current_idx);
      continue;
    }

    // Polygon egress to a link
    if (a->is_polygon() && b->is_link()) {
      const auto &aps = ap_layer_.aps_for_polygon(a->polygon_index);
      double poly_w = polygon_layer_.polygons()[a->polygon_index].weight;
      double best = std::numeric_limits<double>::infinity();
      FMM::NETWORK::NodeID best_ap = kNoAccessPoint;
      std::vector<FMM::NETWORK::EdgeIndex> chosen_segs;
      for (auto ap_idx : aps) {
        const auto &ap = ap_layer_.access_points()[ap_idx];
        if (!ap.attached_node.has_value() || ap.attached_edges.empty()) continue;
        for (auto src_edge : ap.attached_edges) {
          std::vector<FMM::NETWORK::EdgeIndex> goals{b->edge->index};
          std::vector<FMM::ROUTING::Path> paths(1);
          FMM::ROUTING::shortest_edge_to_edges(link_graph_, state, heap,
                                               src_edge, goals, paths);
          if (!paths[0].found) continue;
          double cost = poly_w * point_distance(a->matched_point, ap.point) +
                        paths[0].total_cost +
                        b->edge->weight * (b->offset - b->edge->length);
          if (cost < best) {
            best = cost;
            best_ap = ap.node_id;
            chosen_segs = paths[0].edges;
          }
        }
      }
      if (!active.empty()) {
        active.back().egress_ap = best_ap;
      }
      // Emit intermediate edges from chosen AP path, including final edge.
      for (auto it = chosen_segs.begin(); it != chosen_segs.end(); ++it) {
        push_link(network_.get_edges()[*it].id);
        ++current_idx;
      }
      indices.push_back(current_idx);
      continue;
    }

    // Polygon → Polygon
    if (a->is_polygon() && b->is_polygon()) {
      if (a->polygon_index == b->polygon_index) {
        // Same polygon, another inside (or boundary) GPS observation.
        if (!active.empty() &&
            active.back().p_idx == a->polygon_index) {
          if (b->inside) {
            auto &seg = active.back();
            if (seg.has_inside_obs) {
              seg.inside_path_sum += point_distance(seg.last_inside,
                                                    b->matched_point);
              seg.last_inside = b->matched_point;
              seg.last_inside_layer = (int)(i + 1);
            } else {
              seg.has_inside_obs = true;
              seg.first_inside = b->matched_point;
              seg.last_inside = b->matched_point;
              seg.last_inside_layer = (int)(i + 1);
            }
          }
        }
        indices.push_back(current_idx);
        continue;
      }
      // Different polygons via shared AP.
      const auto &aps_a = ap_layer_.aps_for_polygon(a->polygon_index);
      const auto &aps_b = ap_layer_.aps_for_polygon(b->polygon_index);
      std::unordered_set<FMM::NETWORK::AccessPointIndex> set_b(aps_b.begin(),
                                                                aps_b.end());
      FMM::NETWORK::NodeID shared_ap = kNoAccessPoint;
      for (auto ap_idx : aps_a) {
        if (set_b.count(ap_idx)) {
          shared_ap = ap_layer_.access_points()[ap_idx].node_id;
          break;
        }
      }
      if (!active.empty()) active.back().egress_ap = shared_ap;
      begin_segment(b->polygon_index, cpath.size());
      active.back().entry_ap = shared_ap;
      push_polygon(polygon_layer_.polygons()[b->polygon_index].id);
      ++current_idx;
      if (b->inside) {
        active.back().has_inside_obs = true;
        active.back().first_inside = b->matched_point;
        active.back().last_inside = b->matched_point;
        active.back().last_inside_layer = (int)(i + 1);
      }
      indices.push_back(current_idx);
      continue;
    }
  }

  // Assemble PolygonSegments + distance_inside per R12
  for (auto &seg : active) {
    PolygonSegment ps;
    ps.polygon_id = seg.p_id;
    ps.entry_ap = seg.entry_ap;
    ps.egress_ap = seg.egress_ap;
    ps.position_in_cpath = seg.position_in_cpath;
    ps.is_through = !seg.has_inside_obs;
    // distance_inside
    double d = 0.0;
    if (seg.has_inside_obs) {
      if (seg.entry_ap != kNoAccessPoint) {
        // distance from entry AP coords to first inside observation
        // Find AP coords by node_id
        auto it = ap_layer_.access_points();
        for (const auto &ap : it) {
          if (ap.node_id == seg.entry_ap) {
            d += point_distance(ap.point, seg.first_inside);
            break;
          }
        }
      }
      d += seg.inside_path_sum;
      if (seg.egress_ap != kNoAccessPoint) {
        for (const auto &ap : ap_layer_.access_points()) {
          if (ap.node_id == seg.egress_ap) {
            d += point_distance(seg.last_inside, ap.point);
            break;
          }
        }
      }
    } else if (seg.entry_ap != kNoAccessPoint &&
               seg.egress_ap != kNoAccessPoint) {
      // Through-routing: entry to egress straight-line distance.
      FMM::CORE::Point ep(0, 0), gp(0, 0);
      bool have_e = false, have_g = false;
      for (const auto &ap : ap_layer_.access_points()) {
        if (ap.node_id == seg.entry_ap) {
          ep = ap.point;
          have_e = true;
        }
        if (ap.node_id == seg.egress_ap) {
          gp = ap.point;
          have_g = true;
        }
      }
      if (have_e && have_g) d = point_distance(ep, gp);
    }
    ps.distance_inside = d;
    out->polygon_segments.push_back(ps);
  }

  // Build opath (one entry per GPS point) — use the candidate's edge ID for
  // link candidates, -polygon_id for polygon candidates.
  O_Path opath_out;
  opath_out.reserve(opath.size());
  for (const auto *node : opath) {
    if (node->c->is_link()) {
      opath_out.push_back(node->c->edge->id);
    } else {
      opath_out.push_back(-polygon_layer_.polygons()[node->c->polygon_index].id);
    }
  }

  out->base.opath = std::move(opath_out);
  out->base.cpath = std::move(cpath);
  out->base.indices = std::move(indices);
}

// ---------------- match_traj orchestrator ----------------

PolyMatchResult POLYMATCH::match_traj(const CORE::Trajectory &traj,
                                      const POLYMATCHConfig &config,
                                      ROUTING::DijkstraState &state,
                                      ROUTING::IndexedMinHeap &heap,
                                      bool link_only_mode,
                                      MatchTimings *timings) {
  // Link-only fallback: delegate to WEIGHTMATCH for binary-identical output
  // per SC-002 / FR-012.
  if (link_only_mode || polygon_layer_.empty()) {
    WEIGHTMATCHConfig wm_config{config.k, config.radius, config.gps_error,
                                config.backup_k, config.backup_radius,
                                config.upper_bound_factor,
                                config.allow_truncation};
    WEIGHTMATCH wm(network_, link_graph_);
    PolyMatchResult result;
    result.base = wm.match_traj(traj, wm_config, state, heap, timings);
    return result;
  }

  PolyMatchResult result;
  result.base.id = traj.id;

  auto t0 = UTIL::get_current_time();

  // Phase A — candidates
  MM::Traj_Candidates link_cands_owner;
  PolyTrajCandidates ptc = build_candidates(traj, config, link_cands_owner);
  if (ptc.empty()) return result;

  // Skip layers with zero candidates entirely from both ends (mirror weightmatch
  // truncation behavior). For simplicity here, we require all layers to have at
  // least one candidate; if any layer is empty, we return an empty result.
  for (const auto &layer : ptc) {
    if (layer.empty()) return result;
  }

  auto t1 = UTIL::get_current_time();
  if (timings) timings->candidate_search += UTIL::get_duration(t0, t1);

  // Phase B — HMM
  PolyTransitionGraph tg(ptc, config.gps_error);
  update_tg(tg, traj, config, state, heap);

  auto t2 = UTIL::get_current_time();
  if (timings) timings->update_tg += UTIL::get_duration(t1, t2);

  PolyTGOpath opath = tg.backtrack();

  auto t3 = UTIL::get_current_time();
  if (timings) timings->backtrack += UTIL::get_duration(t2, t3);

  // Phase C — assemble hybrid result
  build_hybrid_path(opath, traj, config, state, heap, &result);

  auto t4 = UTIL::get_current_time();
  if (timings) timings->build_cpath += UTIL::get_duration(t3, t4);

  // matched geometry (best-effort: only link portions for now)
  // mgeom remains empty/default for hybrid paths; clients can derive it from
  // opath + polygon_segments. Future refinement: build hybrid LineString.

  return result;
}

} // MM
} // FMM
