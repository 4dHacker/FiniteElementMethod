#pragma once
#include <cmath>

namespace fem {

// Model problem: -Laplace(u) = f in Omega, u = g on the boundary.
// Replace with your own differential equation.
struct PoissonProblem {
  static double f(double /*x*/, double /*y*/) { return 1.0; }
  static double g(double /*x*/, double /*y*/) { return 0.0; }
};

}  // namespace fem
