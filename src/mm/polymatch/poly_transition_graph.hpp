#ifndef FMM_MM_POLYMATCH_POLY_TRANSITION_GRAPH_HPP_
#define FMM_MM_POLYMATCH_POLY_TRANSITION_GRAPH_HPP_

#include <vector>

#include "mm/polymatch/poly_candidate.hpp"

namespace FMM {
namespace MM {

// Polymatch-internal Viterbi node, parallel to FMM::MM::TGNode but referencing
// PolyCandidate (which may be a link OR polygon snap).
struct PolyTGNode {
  const PolyCandidate *c = nullptr;
  PolyTGNode *prev = nullptr;
  double ep = 0.0;              // emission probability
  double tp = 0.0;              // transition probability from prev
  double cumu_prob = 0.0;       // cumulative log-probability
  double sp_dist = 0.0;         // shortest-path cost from prev (for reporting)
};

using PolyTGLayer = std::vector<PolyTGNode>;
using PolyTGOpath = std::vector<const PolyTGNode *>;

// Parallel to FMM::MM::TransitionGraph. Constructs layers from
// PolyTrajCandidates, computes emission probability per layer-0 candidate, and
// supports Viterbi backtracking.
class PolyTransitionGraph {
 public:
  PolyTransitionGraph(const PolyTrajCandidates &tc, double gps_error);

  static double calc_tp(double sp_dist, double eu_dist);
  static double calc_ep(double dist, double error);

  std::vector<PolyTGLayer> &get_layers();
  PolyTGOpath backtrack();

 private:
  std::vector<PolyTGLayer> layers_;
};

} // MM
} // FMM

#endif // FMM_MM_POLYMATCH_POLY_TRANSITION_GRAPH_HPP_
