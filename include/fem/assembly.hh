#pragma once
#include <cmath>

#include "fem/mesh.hh"
#include "hdnum.hh"

namespace fem {

struct LinearSystem {
  hdnum::SparseMatrix<double> A;
  hdnum::Vector<double> b;
};

// To Implement: assemble the P1 stiffness matrix and load vector for -Laplace(u) = f
// and apply Dirichlet conditions u = g on boundary points.
template <class Problem>
LinearSystem assemble_poisson(const Mesh& mesh, const Problem& problem) {
  const std::size_t n = mesh.points.size();
  hdnum::SparseMatrix<double>::builder builder(n, n);
  hdnum::Vector<double> b(n, 0.0);
  (void)problem;
  return {builder.build(), b};
}

} 
