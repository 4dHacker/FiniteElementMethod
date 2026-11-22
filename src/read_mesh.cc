#include <algorithm>
#include <exception>
#include <iostream>

#include "fem/mesh_reader.hh"

// Reader-only demo: deliberately does not call the unfinished assembly/solver.
int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: read_mesh <mesh.msh>\n";
    return 1;
  }
  try {
    const auto mesh = fem::read_gmsh_msh2(std::string(argv[1]));
    std::cout << "Points: " << mesh.points.size()
              << "\nTriangles: " << mesh.triangles.size()
              << "\nBoundary points: "
              << std::count(mesh.on_boundary.begin(), mesh.on_boundary.end(), true) << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
