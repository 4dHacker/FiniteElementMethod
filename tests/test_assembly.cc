#include <cmath>
#include <iostream>

#include "fem/assembly.hh"

// Unit square, 4 triangles around the center point (index 4, the only interior point).
// Exact solution u = x(1-x) + y(1-y) has -Laplace(u) = 4; g = u on the boundary.
// P1 on this coarse mesh is not exact, so we compare against the hand-computed discrete value.
struct Exact {
  static double f(double, double) { return 4.0; }
  static double g(double x, double y) { return x * (1 - x) + y * (1 - y); }
};

int main() {
  fem::Mesh m;
  m.points = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}};
  m.triangles = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
  m.on_boundary = {true, true, true, true, false};

  auto sys = fem::assemble_poisson(m, Exact{});
  auto LU = hdnum::sparseLUDecomposition(sys.A);
  auto u = hdnum::sparseLUSolver(LU, sys.b);

  // Hand computation: K_cc = 4, corner values g = 0, load = 4 * (1/3 total area) -> u_c = 1/3
  const double err = std::abs(u[4] - 1.0 / 3.0);
  std::cout << "u(center) = " << u[4] << ", error = " << err << "\n";
  return err < 1e-12 ? 0 : 1;
}
