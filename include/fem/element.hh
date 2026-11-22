#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fem {

using Point2 = std::array<double, 2>;
using Matrix2 = std::array<std::array<double, 2>, 2>;
using Matrix3 = std::array<std::array<double, 3>, 3>;

// Reference points and weights from Bastian's 4 November 2026 handout.
// Weights include the reference area 1/2; they do not sum to 1.
inline constexpr std::array<Point2, 3> p1_quadrature_points = {{
    {1.0 / 6, 1.0 / 6}, {2.0 / 3, 1.0 / 6}, {1.0 / 6, 2.0 / 3}}};
inline constexpr std::array<double, 3> p1_quadrature_weights = {1.0 / 6, 1.0 / 6, 1.0 / 6};
inline constexpr std::array<Point2, 3> p1_reference_gradients = {{{-1, -1}, {1, 0}, {0, 1}}};

inline std::array<double, 3> p1_shape_values(const Point2& reference_point) {
  return {1 - reference_point[0] - reference_point[1], reference_point[0], reference_point[1]};
}

// Small, fixed-size local data; HDNUM is used later for global sparse algebra.
// Matrices are indexed [row][column]. Corner order is preserved, including CW.
struct P1Element {
  Matrix2 jacobian{};
  double signed_determinant = 0;
  double area = 0;
  Matrix2 inverse_transpose{};
  std::array<Point2, 3> gradients{};
  std::array<Point2, 3> quadrature_points{};  // Physical, not reference, coordinates.
  std::array<double, 3> source_values{};
  Matrix3 stiffness{};
  std::array<double, 3> load{};
};

namespace element_detail {
inline void finite(double value) {
  if (!std::isfinite(value))
    throw std::runtime_error("P1 element: non-finite input or arithmetic result");
}
}  // namespace element_detail

// Local Poisson element for -Laplace(u)=f, with scalar diffusion K=1.
// Source must be callable as source(x,y). No boundary conditions are applied.
template <class Source>
P1Element compute_p1_element(const std::array<Point2, 3>& vertices, const Source& source) {
  using element_detail::finite;
  for (const auto& vertex : vertices)
    for (double coordinate : vertex) finite(coordinate);

  P1Element element;
  const double a = vertices[1][0] - vertices[0][0];
  const double b = vertices[2][0] - vertices[0][0];
  const double c = vertices[1][1] - vertices[0][1];
  const double d = vertices[2][1] - vertices[0][1];
  for (double value : {a, b, c, d}) finite(value);
  element.jacobian = {{{a, b}, {c, d}}};

  // Relative shape check: uniformly small/large triangles are not automatically
  // treated as degenerate. Extremely flat cells are rejected at double precision.
  const double scale = std::max({std::abs(a), std::abs(b), std::abs(c), std::abs(d)});
  if (scale == 0) throw std::runtime_error("P1 element: degenerate triangle");
  const double scaled_det = (a / scale) * (d / scale) - (b / scale) * (c / scale);
  if (std::abs(scaled_det) <= 64 * std::numeric_limits<double>::epsilon())
    throw std::runtime_error("P1 element: degenerate or numerically flat triangle");

  const double determinant = a * d - b * c;
  finite(determinant);
  element.signed_determinant = determinant;
  element.area = 0.5 * std::abs(determinant);
  if (element.area == 0) throw std::runtime_error("P1 element: triangle area underflow");
  element.inverse_transpose = {{{d / determinant, -c / determinant},
                                {-b / determinant, a / determinant}}};
  for (const auto& row : element.inverse_transpose)
    for (double value : row) finite(value);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t axis = 0; axis < 2; ++axis) {
      element.gradients[i][axis] = element.inverse_transpose[axis][0] * p1_reference_gradients[i][0]
                                + element.inverse_transpose[axis][1] * p1_reference_gradients[i][1];
      finite(element.gradients[i][axis]);
    }
  }
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      element.stiffness[i][j] = element.area
          * (element.gradients[i][0] * element.gradients[j][0]
             + element.gradients[i][1] * element.gradients[j][1]);
      finite(element.stiffness[i][j]);
    }
  }

  for (std::size_t q = 0; q < 3; ++q) {
    const auto& reference = p1_quadrature_points[q];
    const Point2 physical = {vertices[0][0] + a * reference[0] + b * reference[1],
                             vertices[0][1] + c * reference[0] + d * reference[1]};
    for (double coordinate : physical) finite(coordinate);
    element.quadrature_points[q] = physical;
    const double value = source(physical[0], physical[1]);
    finite(value);
    element.source_values[q] = value;
    const auto phi = p1_shape_values(reference);
    for (std::size_t i = 0; i < 3; ++i) {
      element.load[i] += std::abs(determinant) * p1_quadrature_weights[q] * value * phi[i];
      finite(element.load[i]);
    }
  }
  return element;
}
}  // namespace fem
