#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace fem {

// Indices are 0-based.
struct Mesh {
  std::vector<std::array<double, 2>> points;          // N x 2 coordinates
  std::vector<std::array<std::size_t, 3>> triangles;  // M x 3 corner point indices
  std::vector<bool> on_boundary;                      // size N, true for boundary points
};

}  // namespace fem
