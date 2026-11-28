#pragma once
#include <array>
#include <cstddef>
#include <vector>

// This header contains the definition of the data structure for the mesh. It includes three lists:
// 1. points, which stores the coordinates of each vertex
// 2. triangles, which contains the global vertex indices of each triangle
// 3. on_boundary, which tells for each vertex whether it lies on the boundary

namespace fem {

// Note: Indices are 0-based.
struct Mesh {
  std::vector<std::array<double, 2>> points;          // 1. N x 2 coordinates
  std::vector<std::array<std::size_t, 3>> triangles;  // 2. M x 3 corner point indices
  std::vector<bool> on_boundary;                      // 3. size N, true for boundary points
};

}
