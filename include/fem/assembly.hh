#pragma once
#include <cmath>

#include "fem/mesh.hh"
#include "hdnum.hh"

namespace fem {

struct LinearSystem {
  hdnum::SparseMatrix<double> A;
  hdnum::Vector<double> b;
};

// Assembles the P1 stiffness matrix and load vector for -Laplace(u) = f
// and applies Dirichlet conditions u = g on boundary points (rows/cols eliminated, symmetric).
template <class Problem>
LinearSystem assemble_poisson(const Mesh& mesh, const Problem& problem) {
  const std::size_t n = mesh.points.size();
  hdnum::SparseMatrix<double>::builder builder(n, n);
  hdnum::Vector<double> b(n, 0.0);

  for (const auto& tri : mesh.triangles) {
    const auto& p0 = mesh.points[tri[0]];
    const auto& p1 = mesh.points[tri[1]];
    const auto& p2 = mesh.points[tri[2]];

    // Gradients of the barycentric basis functions are constant: grad_i = (bi, ci) / (2*area)
    const double bs[3] = {p1[1] - p2[1], p2[1] - p0[1], p0[1] - p1[1]};
    const double cs[3] = {p2[0] - p1[0], p0[0] - p2[0], p1[0] - p0[0]};
    const double det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
    const double area = 0.5 * std::abs(det);

    // 3-point edge-midpoint rule (exact for quadratics) for the load vector
    const std::array<double, 2> mid[3] = {
        {0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])},
        {0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1])},
        {0.5 * (p2[0] + p0[0]), 0.5 * (p2[1] + p0[1])}};
    const double fsum = problem.f(mid[0][0], mid[0][1]) + problem.f(mid[1][0], mid[1][1]) +
                        problem.f(mid[2][0], mid[2][1]);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        builder[tri[i]][tri[j]] += (bs[i] * bs[j] + cs[i] * cs[j]) / (4.0 * area);
      // Basis function i has value 1/2 at the two adjacent midpoints and 0 at the opposite one
      const double fi = 0.5 * (fsum - problem.f(mid[(i + 1) % 3][0], mid[(i + 1) % 3][1]));
      b[tri[i]] += area / 3.0 * fi;
    }
  }

  // Dirichlet: move known values to the right-hand side, then set identity rows
  std::vector<double> g(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    if (mesh.on_boundary[i]) g[i] = problem.g(mesh.points[i][0], mesh.points[i][1]);

  for (std::size_t i = 0; i < n; ++i) {
    if (mesh.on_boundary[i]) continue;
    for (auto& [j, v] : builder[i])
      if (mesh.on_boundary[j]) b[i] -= v * g[j];
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!mesh.on_boundary[i]) {
      for (auto it = builder[i].begin(); it != builder[i].end();) {
        if (mesh.on_boundary[it->first]) it = builder[i].erase(it);
        else ++it;
      }
    } else {
      builder[i].clear();
      builder[i][i] = 1.0;
      b[i] = g[i];
    }
  }

  return {builder.build(), b};
}

}  // namespace fem
