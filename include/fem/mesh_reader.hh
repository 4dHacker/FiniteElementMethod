#pragma once
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include "fem/mesh.hh"

namespace fem {

// Reads an ASCII Gmsh file in format 2.2 (generate with: gmsh -2 -format msh2 x.geo).
// Triangles (type 2) become cells; line elements (type 1) mark boundary points.
inline Mesh read_gmsh_msh2(const std::string& filename) {
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("Cannot open mesh file: " + filename);

  Mesh mesh;
  std::map<std::size_t, std::size_t> id_to_index;  // Gmsh ids are 1-based, maybe sparse
  std::vector<std::array<std::size_t, 2>> lines;
  std::string tok;

  while (in >> tok) {
    if (tok == "$MeshFormat") {
      double version;
      int file_type;
      in >> version >> file_type;
      if (version < 2.0 || version >= 3.0 || file_type != 0)
        throw std::runtime_error("Need ASCII msh format 2 (gmsh -format msh2)");
    } else if (tok == "$Nodes") {
      std::size_t n;
      in >> n;
      for (std::size_t i = 0; i < n; ++i) {
        std::size_t id;
        double x, y, z;
        in >> id >> x >> y >> z;
        id_to_index[id] = mesh.points.size();
        mesh.points.push_back({x, y});
      }
    } else if (tok == "$Elements") {
      std::size_t n;
      in >> n;
      for (std::size_t i = 0; i < n; ++i) {
        std::size_t id, type, ntags, tag;
        in >> id >> type >> ntags;
        for (std::size_t t = 0; t < ntags; ++t) in >> tag;
        if (type == 2) {
          std::size_t a, b, c;
          in >> a >> b >> c;
          mesh.triangles.push_back({id_to_index.at(a), id_to_index.at(b), id_to_index.at(c)});
        } else if (type == 1) {
          std::size_t a, b;
          in >> a >> b;
          lines.push_back({id_to_index.at(a), id_to_index.at(b)});
        } else {
          std::string rest;
          std::getline(in, rest);  // skip points and unsupported element types
        }
      }
    }
  }

  mesh.on_boundary.assign(mesh.points.size(), false);
  for (auto& l : lines) mesh.on_boundary[l[0]] = mesh.on_boundary[l[1]] = true;
  return mesh;
}

}  // namespace fem
