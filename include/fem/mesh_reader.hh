#pragma once
#include "fem/mesh.hh"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <istream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {

namespace gmsh_detail {
inline void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error("Gmsh reader: " + message);
}

inline std::string line(std::istream& input) {
  std::string value;
  while (std::getline(input, value)) {
    const auto begin = value.find_first_not_of(" \t\r");
    if (begin != std::string::npos)
      return value.substr(begin, value.find_last_not_of(" \t\r") - begin + 1);
  }
  throw std::runtime_error("Gmsh reader: unexpected end of file");
}

template<class... Values>
inline void parse(const std::string& text, Values&... values) {
  std::istringstream row(text);
  require(static_cast<bool>((row >> ... >> values)), "malformed record: " + text);
  row >> std::ws;
  require(row.eof(), "unexpected fields: " + text);
}

inline std::int64_t count(std::istream& input) {
  std::int64_t n;
  parse(line(input), n);
  require(n >= 0, "negative record count");
  return n;
}

inline void end(std::istream& input, const std::string& section) {
  require(line(input) == "$End" + section, "missing $End" + section);
}
} // namespace gmsh_detail

/** Read Gmsh MSH 2.2 ASCII, supporting element types 1 (line), 2 (triangle),
 * and 15 (point). Point elements are validated but not retained. Other element
 * types are rejected to avoid silently reading a partial or higher-order mesh.
 * Nodes must lie in z=0. Nondegenerate triangle ordering is preserved.
 * Type-1 line endpoints mark on_boundary, following the team's Mesh contract.
 * Physical names/tags are validated but are not exposed by that contract.
 * Throws std::runtime_error on invalid input or unsupported formats.
 */
inline Mesh read_gmsh_msh2(std::istream& input) {
  using namespace gmsh_detail;
  using Id = std::int64_t;
  require(line(input) == "$MeshFormat", "expected $MeshFormat at start");
  std::string version;
  int binary, data_size;
  parse(line(input), version, binary, data_size);
  require(version == "2.2" && binary == 0 && data_size == 8,
          "only MSH 2.2 ASCII is supported; export with gmsh -format msh2");
  end(input, "MeshFormat");

  Mesh mesh;
  std::map<std::pair<int, Id>, std::string> physical_names;
  std::map<Id, std::size_t> node_index;
  struct Element { Id id; int type; Id physical; Id entity; std::vector<Id> nodes; };
  std::vector<Element> elements;
  std::set<std::string> seen;
  while (true) {
    input >> std::ws;
    if (input.eof()) break;
    require(static_cast<bool>(input), "input stream failure");
    const auto section = line(input);
    require(section.size() > 1 && section[0] == '$' && section.rfind("$End", 0) != 0,
            "expected section header, got " + section);
    const bool first = seen.insert(section).second;
    if (section == "$Nodes" || section == "$Elements" || section == "$PhysicalNames")
      require(first, "duplicate section " + section);
    require(section != "$MeshFormat", "duplicate $MeshFormat");
    if (section == "$Nodes") {
      const auto n = count(input);
      for (std::int64_t i = 0; i < n; ++i) {
        Id id;
        double x, y, z;
        parse(line(input), id, x, y, z);
        require(id > 0, "node ID must be positive");
        require(std::isfinite(x) && std::isfinite(y) && std::isfinite(z),
                "non-finite node coordinate");
        require(z == 0.0, "only meshes in the xy plane (z=0) are supported");
        require(node_index.emplace(id, mesh.points.size()).second, "duplicate node ID");
        mesh.points.push_back({x, y});
      }
      end(input, "Nodes");
    } else if (section == "$Elements") {
      const auto n = count(input);
      std::set<Id> ids;
      for (std::int64_t i = 0; i < n; ++i) {
        const auto text = line(input);
        std::istringstream row(text);
        Element element{};
        std::int64_t ntags;
        require(static_cast<bool>(row >> element.id >> element.type >> ntags),
                "malformed element header");
        require(element.id > 0 && ids.insert(element.id).second,
                "invalid or duplicate element ID");
        require(ntags >= 0, "negative element tag count");
        require(element.type == 1 || element.type == 2 || element.type == 15,
                "unsupported element type " + std::to_string(element.type)
                + "; use first-order triangles and lines");
        for (std::int64_t tag = 0; tag < ntags; ++tag) {
          Id value;
          require(static_cast<bool>(row >> value), "missing element tag");
          if (tag == 0) element.physical = value;
          if (tag == 1) element.entity = value;
        }
        const int nnodes = element.type == 2 ? 3 : (element.type == 1 ? 2 : 1);
        for (int j = 0; j < nnodes; ++j) {
          Id node;
          require(static_cast<bool>(row >> node), "missing element node");
          element.nodes.push_back(node);
        }
        row >> std::ws;
        require(row.eof(), "unexpected fields in element record");
        elements.push_back(element);
      }
      end(input, "Elements");
    } else if (section == "$PhysicalNames") {
      const auto n = count(input);
      for (std::int64_t i = 0; i < n; ++i) {
        std::istringstream row(line(input));
        int dimension;
        Id tag;
        std::string name;
        require(static_cast<bool>(row >> dimension >> tag), "malformed physical group");
        row >> std::ws;
        require(row.peek() == '"', "physical group name must be quoted");
        require(static_cast<bool>(row >> std::quoted(name)), "invalid physical group name");
        row >> std::ws;
        require(row.eof() && dimension >= 0 && dimension <= 2 && tag > 0,
                "invalid physical group record");
        require(physical_names.emplace(std::make_pair(dimension, tag), name).second,
                "duplicate physical group");
      }
      end(input, "PhysicalNames");
    } else {
      // Additional sections (e.g. NodeData) are not needed for mesh geometry.
      const auto closing = "$End" + section.substr(1);
      while (line(input) != closing) {}
    }
  }
  require(seen.count("$Nodes") && seen.count("$Elements"), "missing Nodes or Elements section");
  require(!mesh.points.empty(), "mesh has no vertices");
  mesh.on_boundary.assign(mesh.points.size(), false);
  for (const auto& element : elements) {
    std::vector<std::size_t> indices;
    for (const auto node : element.nodes) {
      const auto it = node_index.find(node);
      require(it != node_index.end(), "element references unknown node " + std::to_string(node));
      indices.push_back(it->second);
    }
    if (element.type == 2) {
      const auto& a = mesh.points[indices[0]];
      const auto& b = mesh.points[indices[1]];
      const auto& c = mesh.points[indices[2]];
      const long double bx = static_cast<long double>(b[0]) - a[0];
      const long double by = static_cast<long double>(b[1]) - a[1];
      const long double cx = static_cast<long double>(c[0]) - a[0];
      const long double cy = static_cast<long double>(c[1]) - a[1];
      const auto determinant = bx * cy - by * cx;
      require(std::isfinite(determinant) && determinant != 0, "degenerate triangle");
      mesh.triangles.push_back({indices[0], indices[1], indices[2]});
    } else if (element.type == 1) {
      const auto& a = mesh.points[indices[0]];
      const auto& b = mesh.points[indices[1]];
      require(a != b, "zero-length line element");
      mesh.on_boundary[indices[0]] = mesh.on_boundary[indices[1]] = true;
    }
  }
  require(!mesh.triangles.empty(), "mesh has no linear triangles");
  return mesh;
}

inline Mesh read_gmsh_msh2(const std::string& filename) {
  std::ifstream input(filename);
  gmsh_detail::require(input.is_open(), "cannot open file " + filename);
  return read_gmsh_msh2(input);
}
} // namespace fem
