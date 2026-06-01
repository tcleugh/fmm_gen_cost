#ifndef FMM_NETWORK_TRACE_CATEGORY_HPP_
#define FMM_NETWORK_TRACE_CATEGORY_HPP_

#include <array>
#include <string>

namespace FMM {
namespace NETWORK {

// The ten labels from specs/002-real-network-validation FR-004.
// Stored in the trips.csv `category` column. Single-sourced here so the
// generator (src/app/polymatch_traces_gen.cpp + src/network/trace_generator.cpp)
// and the validation harness (test/polymatch_test.cpp) agree on spelling.
enum class TraceCategory : uint8_t {
  LinkOnly = 0,
  PolygonTraversal,
  PolygonSharedAp,
  MidPolygonStart,
  MidPolygonEnd,
  FullyInside,
  ThroughRouting,
  OffNetworkNoise,
  ShortTrip,
  DuplicatePoints,
};

constexpr std::array<const char*, 10> kCategoryLabels = {
  "link-only",
  "polygon-traversal",
  "polygon-shared-ap",
  "mid-polygon-start",
  "mid-polygon-end",
  "fully-inside",
  "through-routing",
  "off-network-noise",
  "short-trip",
  "duplicate-points",
};

inline const char* to_label(TraceCategory c) {
  return kCategoryLabels[static_cast<size_t>(c)];
}

// Convenience for harness CSV-parse paths that come in as raw strings.
inline bool is_valid_label(const std::string& s) {
  for (const char* lbl : kCategoryLabels) {
    if (s == lbl) return true;
  }
  return false;
}

} // namespace NETWORK
} // namespace FMM

#endif // FMM_NETWORK_TRACE_CATEGORY_HPP_
