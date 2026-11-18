#pragma once
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fem/mesh.hh"

namespace fem {
namespace detail {

inline void read_msh2_nodes(std::istream& in, Mesh& mesh,
                           std::map<std::size_t, std::size_t>& id_to_index) {
  std::size_t n;
  in >> n;
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t id;
    double x, y, z;
    in >> id >> x >> y >> z;
    id_to_index[id] = mesh.points.size();
    mesh.points.push_back({x, y});
  }
}

inline void read_msh4_nodes(std::istream& in, Mesh& mesh,
                            std::map<std::size_t, std::size_t>& id_to_index) {
  std::size_t num_blocks, num_nodes, min_tag, max_tag;
  in >> num_blocks >> num_nodes >> min_tag >> max_tag;
  for (std::size_t block = 0; block < num_blocks; ++block) {
    int entity_dim, entity_tag, parametric;
    std::size_t block_size;
    in >> entity_dim >> entity_tag >> parametric >> block_size;

    std::vector<std::size_t> node_tags(block_size);
    for (auto& tag : node_tags) in >> tag;

    for (const auto tag : node_tags) {
      double x, y, z;
      in >> x >> y >> z;
      id_to_index[tag] = mesh.points.size();
      mesh.points.push_back({x, y});
      if (parametric) {
        double parameter;
        for (int i = 0; i < entity_dim; ++i) in >> parameter;
      }
    }
  }
}

inline void read_msh2_elements(std::istream& in, Mesh& mesh,
                               const std::map<std::size_t, std::size_t>& id_to_index,
                               std::vector<std::array<std::size_t, 2>>& lines) {
  std::size_t n;
  in >> n;
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t id, type, num_tags, tag;
    in >> id >> type >> num_tags;
    for (std::size_t t = 0; t < num_tags; ++t) in >> tag;
    if (type == 2) {
      std::size_t a, b, c;
      in >> a >> b >> c;
      mesh.triangles.push_back(
          {id_to_index.at(a), id_to_index.at(b), id_to_index.at(c)});
    } else if (type == 1) {
      std::size_t a, b;
      in >> a >> b;
      lines.push_back({id_to_index.at(a), id_to_index.at(b)});
    } else {
      std::string rest;
      std::getline(in, rest);
    }
  }
}

inline void read_msh4_elements(std::istream& in, Mesh& mesh,
                               const std::map<std::size_t, std::size_t>& id_to_index,
                               std::vector<std::array<std::size_t, 2>>& lines) {
  std::size_t num_blocks, num_elements, min_tag, max_tag;
  in >> num_blocks >> num_elements >> min_tag >> max_tag;
  std::string row;
  std::getline(in, row);

  for (std::size_t block = 0; block < num_blocks; ++block) {
    int entity_dim, entity_tag, element_type;
    std::size_t block_size;
    in >> entity_dim >> entity_tag >> element_type >> block_size;
    std::getline(in, row);

    for (std::size_t i = 0; i < block_size; ++i) {
      if (!std::getline(in, row))
        throw std::runtime_error("Malformed Gmsh 4.x element section");
      std::istringstream element(row);
      std::size_t id;
      element >> id;
      if (element_type == 2) {
        std::size_t a, b, c;
        if (!(element >> a >> b >> c))
          throw std::runtime_error("Malformed Gmsh 4.x triangle element");
        mesh.triangles.push_back(
            {id_to_index.at(a), id_to_index.at(b), id_to_index.at(c)});
      } else if (element_type == 1) {
        std::size_t a, b;
        if (!(element >> a >> b))
          throw std::runtime_error("Malformed Gmsh 4.x line element");
        lines.push_back({id_to_index.at(a), id_to_index.at(b)});
      }
    }
  }
}

}  // namespace detail

inline Mesh read_gmsh_msh(const std::string& filename) {
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("Cannot open mesh file: " + filename);

  Mesh mesh;
  std::map<std::size_t, std::size_t> id_to_index;
  std::vector<std::array<std::size_t, 2>> lines;
  std::string tok;
  double version = 0.0;

  while (in >> tok) {
    if (tok == "$MeshFormat") {
      int file_type;
      int data_size;
      in >> version >> file_type >> data_size;
      const bool supported_version =
          (version >= 2.0 && version < 3.0) || (version >= 4.0 && version < 5.0);
      if (!supported_version || file_type != 0)
        throw std::runtime_error("Need ASCII Gmsh msh format 2.x or 4.x");
    } else if (tok == "$Nodes") {
      if (version < 4.0)
        detail::read_msh2_nodes(in, mesh, id_to_index);
      else
        detail::read_msh4_nodes(in, mesh, id_to_index);
    } else if (tok == "$Elements") {
      if (version < 4.0)
        detail::read_msh2_elements(in, mesh, id_to_index, lines);
      else
        detail::read_msh4_elements(in, mesh, id_to_index, lines);
    }
  }

  mesh.on_boundary.assign(mesh.points.size(), false);
  for (const auto& line : lines)
    mesh.on_boundary[line[0]] = mesh.on_boundary[line[1]] = true;
  return mesh;
}

inline Mesh read_gmsh_msh2(const std::string& filename) {
  return read_gmsh_msh(filename);
}

}  // namespace fem
