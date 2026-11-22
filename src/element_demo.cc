#include <iomanip>
#include <iostream>
#include <string>

#include "fem/element.hh"

int main(int argc, char** argv) {
  if (argc != 2 || (std::string(argv[1]) != "reference" && std::string(argv[1]) != "skew")) {
    std::cerr << "Usage: element_demo <reference|skew>\n";
    return 1;
  }
  const bool skew = std::string(argv[1]) == "skew";
  const std::array<fem::Point2, 3> vertices = skew
      ? std::array<fem::Point2, 3>{{{1, 1}, {3, 1}, {2, 2}}}
      : std::array<fem::Point2, 3>{{{0, 0}, {1, 0}, {0, 1}}};
  try {
    const auto source = [skew](double x, double y) { return skew ? 1 + x + 2 * y : 1.0; };
    const auto e = fem::compute_p1_element(vertices, source);
    std::cout << std::setprecision(17) << "Case: " << argv[1]
              << "\nSource: " << (skew ? "f(x,y)=1+x+2y" : "f(x,y)=1") << '\n';
    for (std::size_t i = 0; i < 3; ++i)
      std::cout << "P[" << i << "] = " << vertices[i][0] << ' ' << vertices[i][1] << '\n';
    std::cout << "B (rows):\n";
    for (const auto& row : e.jacobian) std::cout << row[0] << ' ' << row[1] << '\n';
    std::cout << "signed det(B) = " << e.signed_determinant << "\narea = " << e.area
              << "\nB^(-T) (rows):\n";
    for (const auto& row : e.inverse_transpose) std::cout << row[0] << ' ' << row[1] << '\n';
    for (std::size_t i = 0; i < 3; ++i)
      std::cout << "G[" << i << "] = " << e.gradients[i][0] << ' ' << e.gradients[i][1] << '\n';
    for (std::size_t q = 0; q < 3; ++q) {
      const auto phi = fem::p1_shape_values(fem::p1_quadrature_points[q]);
      std::cout << "Q[" << q << "] = " << e.quadrature_points[q][0] << ' '
                << e.quadrature_points[q][1] << "; f(Q) = " << e.source_values[q]
                << "; reference weight = " << fem::p1_quadrature_weights[q]
                << "; phi = " << phi[0] << ' ' << phi[1] << ' ' << phi[2] << '\n';
    }
    std::cout << "E (rows):\n";
    for (const auto& row : e.stiffness)
      std::cout << row[0] << ' ' << row[1] << ' ' << row[2] << '\n';
    std::cout << "r = " << e.load[0] << ' ' << e.load[1] << ' ' << e.load[2] << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
