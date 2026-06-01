#include "mm/polymatch/poly_transition_graph.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "util/debug.hpp"

namespace FMM {
namespace MM {

PolyTransitionGraph::PolyTransitionGraph(const PolyTrajCandidates &tc,
                                         double gps_error) {
  layers_.reserve(tc.size());
  for (const auto &point_cands : tc) {
    PolyTGLayer layer;
    layer.reserve(point_cands.size());
    for (const auto &c : point_cands) {
      double ep = calc_ep(c.ep_distance, gps_error);
      layer.push_back(
          PolyTGNode{&c, nullptr, ep, 0.0,
                     -std::numeric_limits<double>::infinity(), 0.0});
    }
    layers_.push_back(std::move(layer));
  }
  if (!layers_.empty()) {
    // Initialize layer-0 cumu_prob to log(ep) and clear back-pointers.
    for (auto &node : layers_.front()) {
      node.cumu_prob = std::log(node.ep);
      node.prev = nullptr;
    }
  }
}

double PolyTransitionGraph::calc_tp(double sp_dist, double eu_dist) {
  if (eu_dist == 0.0) eu_dist = 1e-6;
  return eu_dist >= sp_dist ? 1.0 : eu_dist / sp_dist;
}

double PolyTransitionGraph::calc_ep(double dist, double error) {
  double a = dist / error;
  return std::exp(-0.5 * a * a);
}

std::vector<PolyTGLayer> &PolyTransitionGraph::get_layers() { return layers_; }

PolyTGOpath PolyTransitionGraph::backtrack() {
  PolyTGOpath opath;
  if (layers_.empty()) return opath;

  PolyTGNode *track = nullptr;
  double final_prob = -std::numeric_limits<double>::infinity();
  auto &last_layer = layers_.back();
  for (auto &node : last_layer) {
    if (node.cumu_prob > final_prob) {
      final_prob = node.cumu_prob;
      track = &node;
    }
  }
  if (final_prob == -std::numeric_limits<double>::infinity() || !track) {
    return opath;
  }
  while (track != nullptr) {
    opath.push_back(track);
    track = track->prev;
  }
  std::reverse(opath.begin(), opath.end());
  return opath;
}

} // MM
} // FMM
