#include <iostream>

#include "fem/mesh_reader.hh"
#include "hdnum.hh"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " mesh.msh\n";
    return 1;
  }
  //reading mesh from file
  const fem::Mesh mesh = fem::read_gmsh_msh2(argv[1]);
  const auto& points = mesh.points;
  const auto& triangles = mesh.triangles;

  //right-hand side f(x, y) (example)
  auto f = [](double x, double y) { return 2.0 * x * (1.0 - x) + 2.0 * y * (1.0 - y); };

  std::cout << points.size() << " points, " << triangles.size() << " triangles\n";

  //initialising sparse matrix builder and vectors
  const std::size_t N = points.size();
  hdnum::SparseMatrix<double>::builder builder(N, N);
  hdnum::Vector<double> B(N, 0.0);
  hdnum::Vector<double> z(N, 0.0);

  //generating Matrix
  for (const auto& T : triangles) {
    //getting triangle vertices
    auto [i0, i1, i2] = T;
    const auto& [x0, y0] = points[i0];
    const auto& [x1, y1] = points[i1];
    const auto& [x2, y2] = points[i2];
    
    //area of the triangle
    double area = 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
    
    //Gradients of Hutfunction
    double beta[3] = {y1-y2, y2-y0, y0-y1};
    double gamma[3] = {x2-x1, x0-x2, x1-x0};
    
    //local stiffness matrix for the triangle
    double A_loc[3][3];
    double c = f((x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0)*area/3;
    
    for (int i = 0; i < 3; i++) {
      auto a = T[i];
      B[a] += c;
      for (int j = 0; j < 3; j++) {
        auto b = T[j];
        A_loc[i][j] = (beta[i] * beta[j] + gamma[i] * gamma[j]) / (4.0 * area);
        builder[a][b] += A_loc[i][j];
      }
    }
  }
  auto A = builder.build();
}
