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

  std::cout << points.size() << " points, " << triangles.size() << " triangles\n";

  //initialising sparse matrix builder and vectors
  const std::size_t N = points.size();
  hdnum::SparseMatrix<double>::builder builder(N, N);
  hdnum::Vector<double> b(N, 0.0);
  hdnum::Vector<double> z(N, 0.0);

  //initialising matrix form
  for (const auto& T : triangles) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double a = T[i];
        double b = T[j];
        builder.addEntry(a, b, 1.0);
      }
    } 
  }
  auto A = builder.build();

}
