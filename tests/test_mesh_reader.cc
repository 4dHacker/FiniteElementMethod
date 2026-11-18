#include <algorithm>
#include <iostream>

#include "fem/mesh_reader.hh"

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Expected paths to an MSH2 and an MSH4 mesh\n";
    return 1;
  }

  const auto msh2 = fem::read_gmsh_msh(argv[1]);
  if (msh2.points.size() != 142 || msh2.triangles.empty() ||
      std::none_of(msh2.on_boundary.begin(), msh2.on_boundary.end(),
                   [](bool boundary) { return boundary; })) {
    std::cerr << "MSH2 mesh was not read correctly\n";
    return 1;
  }

  const auto msh4 = fem::read_gmsh_msh(argv[2]);
  if (msh4.points.size() != 3 || msh4.triangles.size() != 1 ||
      std::count(msh4.on_boundary.begin(), msh4.on_boundary.end(), true) != 3 ||
      msh4.triangles[0] != std::array<std::size_t, 3>{0, 1, 2}) {
    std::cerr << "MSH4 mesh was not read correctly\n";
    return 1;
  }

  return 0;
}
