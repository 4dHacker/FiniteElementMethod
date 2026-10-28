#include <iostream>

#include "fem/assembly.hh"
#include "fem/mesh_reader.hh"
#include "fem/problem.hh"
#include "fem/vtk_writer.hh"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " mesh.msh [output.vtk]\n";
    return 1;
  }
  const std::string out = argc > 2 ? argv[2] : "results/solution.vtk";

  fem::Mesh mesh = fem::read_gmsh_msh2(argv[1]);
  std::cout << mesh.points.size() << " points, " << mesh.triangles.size() << " triangles\n";

  auto sys = fem::assemble_poisson(mesh, fem::PoissonProblem{});
  auto LU = hdnum::sparseLUDecomposition(sys.A);
  hdnum::Vector<double> u = hdnum::sparseLUSolver(LU, sys.b);

  fem::write_vtk(out, mesh, u);
  std::cout << "Wrote " << out << "\n";
}
