#include <iostream>
#include <limits>
#include <string>

#include "fem/element.hh"

namespace {
int checks = 0;
void check(bool condition, const std::string& message) {
  ++checks;
  if (!condition) throw std::runtime_error("Test failed: " + message);
}
void near(double actual, double expected, const std::string& message) {
  check(std::isfinite(actual) && std::abs(actual - expected) <= 1e-12 * std::max(1.0, std::abs(expected)),
        message);
}
template <class Function>
void reject(const Function& operation, const std::string& diagnostic) {
  try { operation(); }
  catch (const std::runtime_error& error) {
    check(std::string(error.what()).find(diagnostic) != std::string::npos, diagnostic);
    return;
  }
  throw std::runtime_error("Expected rejection: " + diagnostic);
}
const auto one = [](double, double) { return 1.0; };
const std::array<fem::Point2, 3> reference = {{{0, 0}, {1, 0}, {0, 1}}};
const std::array<fem::Point2, 3> skew = {{{1, 1}, {3, 1}, {2, 2}}};
}

int main() {
  try {
    const auto e = fem::compute_p1_element(reference, one);
    near(e.area, 0.5, "reference area");
    const fem::Matrix3 expected = {{{1, -0.5, -0.5}, {-0.5, 0.5, 0}, {-0.5, 0, 0.5}}};
    for (std::size_t i = 0; i < 3; ++i) {
      near(e.load[i], 1.0 / 6, "constant load");
      for (std::size_t j = 0; j < 3; ++j) near(e.stiffness[i][j], expected[i][j], "reference stiffness");
      for (std::size_t axis = 0; axis < 2; ++axis)
        near(e.gradients[i][axis], fem::p1_reference_gradients[i][axis], "reference gradient");
    }

    // Independent closed-form polynomial moments on the reference triangle.
    std::array<double, 6> moments{};
    for (std::size_t q = 0; q < 3; ++q) {
      const auto [x, y] = fem::p1_quadrature_points[q];
      const std::array<double, 6> monomials = {1, x, y, x*x, x*y, y*y};
      for (std::size_t k = 0; k < 6; ++k) moments[k] += fem::p1_quadrature_weights[q] * monomials[k];
    }
    const std::array<double, 6> integrals = {0.5, 1.0/6, 1.0/6, 1.0/12, 1.0/24, 1.0/12};
    for (std::size_t k = 0; k < 6; ++k) near(moments[k], integrals[k], "quadrature moment");
    for (std::size_t i = 0; i < 3; ++i) {
      const auto phi = fem::p1_shape_values(reference[i]);
      for (std::size_t j = 0; j < 3; ++j) near(phi[j], i == j ? 1 : 0, "nodal basis property");
    }

    const auto affine = [](double x, double y) { return 1 + x + 2*y; };
    const auto s = fem::compute_p1_element(skew, affine);
    near(s.area, 1, "skew area");
    const fem::Matrix2 inverse = {{{0.5, 0}, {-0.5, 1}}};
    const std::array<fem::Point2, 3> gradients = {{{-0.5, -0.5}, {0.5, -0.5}, {0, 1}}};
    const fem::Matrix3 skew_matrix = {{{0.5, 0, -0.5}, {0, 0.5, -0.5}, {-0.5, -0.5, 1}}};
    const std::array<double, 3> skew_load = {7.0/4, 23.0/12, 2};
    const std::array<fem::Point2, 3> physical_points = {{{1.5, 7.0/6}, {2.5, 7.0/6}, {2, 5.0/3}}};
    const std::array<double, 3> source_values = {29.0/6, 35.0/6, 19.0/3};
    for (std::size_t i = 0; i < 2; ++i)
      for (std::size_t j = 0; j < 2; ++j) near(s.inverse_transpose[i][j], inverse[i][j], "inverse transpose");
    for (std::size_t i = 0; i < 3; ++i) {
      near(s.load[i], skew_load[i], "affine load from exact triangle moments");
      near(s.source_values[i], source_values[i], "source at physical point");
      for (std::size_t axis = 0; axis < 2; ++axis) {
        near(s.gradients[i][axis], gradients[i][axis], "skew gradients");
        near(s.quadrature_points[i][axis], physical_points[i][axis], "physical quadrature point");
      }
      double row_sum = 0;
      for (std::size_t j = 0; j < 3; ++j) {
        near(s.stiffness[i][j], skew_matrix[i][j], "skew stiffness");
        near(s.stiffness[i][j], s.stiffness[j][i], "symmetry");
        row_sum += s.stiffness[i][j];
      }
      near(row_sum, 0, "constant function has zero energy");
    }

    // Reordering must only permute local rows/columns, not change physics.
    std::array<std::size_t, 3> permutation = {0, 1, 2};
    do {
      const std::array<fem::Point2, 3> reordered = {skew[permutation[0]], skew[permutation[1]], skew[permutation[2]]};
      const auto p = fem::compute_p1_element(reordered, affine);
      near(p.area, s.area, "area under permutation");
      for (std::size_t i = 0; i < 3; ++i) {
        near(p.load[i], s.load[permutation[i]], "permuted load");
        for (std::size_t j = 0; j < 3; ++j)
          near(p.stiffness[i][j], s.stiffness[permutation[i]][permutation[j]], "permuted stiffness");
      }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    const auto clockwise = fem::compute_p1_element(std::array<fem::Point2, 3>{skew[0], skew[2], skew[1]}, one);
    check(clockwise.signed_determinant < 0 && clockwise.area > 0, "signed determinant and positive area");

    for (double scale : {1e-6, 2.0, 1e6}) {
      auto transformed = reference;
      for (auto& vertex : transformed) {
        const auto original = vertex;
        vertex = {-scale * original[1], scale * original[0]};  // 90-degree rotation and scaling.
      }
      const auto t = fem::compute_p1_element(transformed, one);
      near(t.area, scale*scale/2, "scaled area");
      for (std::size_t i = 0; i < 3; ++i) {
        near(t.load[i] / (scale*scale), e.load[i], "scaled load");
        for (std::size_t j = 0; j < 3; ++j) near(t.stiffness[i][j], e.stiffness[i][j], "scale/rotation invariance in 2D");
      }
    }
    auto shifted = reference;
    for (auto& vertex : shifted) { vertex[0] += 5; vertex[1] -= 7; }
    const auto translated = fem::compute_p1_element(shifted, one);
    for (std::size_t i = 0; i < 3; ++i) {
      near(translated.load[i], e.load[i], "translation-invariant constant load");
      for (std::size_t j = 0; j < 3; ++j) near(translated.stiffness[i][j], e.stiffness[i][j], "translation invariance");
    }
    // Energy of the exactly represented affine field u=2x-3y+4 is area*13.
    double energy = 0;
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j) {
        const double ui = 2*skew[i][0] - 3*skew[i][1] + 4;
        const double uj = 2*skew[j][0] - 3*skew[j][1] + 4;
        energy += ui * s.stiffness[i][j] * uj;
      }
    near(energy, 13 * s.area, "affine-field energy");
    const auto zero = fem::compute_p1_element(reference, [](double, double) { return 0.0; });
    for (double value : zero.load) near(value, 0, "zero source");

    reject([] { fem::compute_p1_element(std::array<fem::Point2, 3>{{{0,0},{1,0},{2,0}}}, one); }, "degenerate");
    reject([] { fem::compute_p1_element(std::array<fem::Point2, 3>{{{0,0},{0,0},{1,1}}}, one); }, "degenerate");
    reject([] { fem::compute_p1_element(std::array<fem::Point2, 3>{{{0,0},{1,0},{1,1e-16}}}, one); }, "flat");
    auto invalid = reference;
    invalid[0][0] = std::numeric_limits<double>::quiet_NaN();
    reject([&] { fem::compute_p1_element(invalid, one); }, "non-finite");
    invalid[0][0] = std::numeric_limits<double>::infinity();
    reject([&] { fem::compute_p1_element(invalid, one); }, "non-finite");
    reject([] { fem::compute_p1_element(reference, [](double, double) { return std::numeric_limits<double>::infinity(); }); }, "non-finite");
    reject([] { fem::compute_p1_element(reference, [](double, double) { return std::numeric_limits<double>::quiet_NaN(); }); }, "non-finite");
    reject([] { fem::compute_p1_element(std::array<fem::Point2, 3>{{{0,0},{1e-200,0},{0,1e-200}}}, one); }, "underflow");
    reject([] { fem::compute_p1_element(std::array<fem::Point2, 3>{{{0,0},{1e308,0},{0,1e308}}}, one); }, "non-finite");
    std::cout << checks << " element checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
