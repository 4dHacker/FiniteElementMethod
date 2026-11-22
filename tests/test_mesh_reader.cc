#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "fem/mesh_reader.hh"

namespace {
int checks = 0;
void check(bool condition, const std::string& message) {
  ++checks;
  if (!condition) throw std::runtime_error("Test failed: " + message);
}
std::string change(std::string text, const std::string& from, const std::string& to) {
  const auto position = text.find(from);
  if (position == std::string::npos) throw std::runtime_error("Invalid test fixture replacement");
  return text.replace(position, from.size(), to);
}
fem::Mesh read(const std::string& text) {
  std::istringstream input(text);
  return fem::read_gmsh_msh2(input);
}
void reject(const std::string& text, const std::string& diagnostic) {
  try { read(text); }
  catch (const std::runtime_error& error) {
    check(std::string(error.what()).find(diagnostic) != std::string::npos,
          "expected " + diagnostic + ", got " + error.what());
    return;
  }
  throw std::runtime_error("Invalid mesh accepted: " + diagnostic);
}

// Four triangles around an interior vertex. IDs are sparse and deliberately
// unsorted, so boundary and connectivity cannot be inferred from ID ordering.
const std::string square = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
1 11 "outer boundary"
2 21 "domain"
$EndPhysicalNames
$Nodes
5
30 0 0 0
10 1 0 0
80 1 1 0
20 0 1 0
90 0.5 0.5 0
$EndNodes
$Elements
8
1 1 2 11 1 30 10
2 1 2 11 2 10 80
3 1 2 11 3 80 20
4 1 2 11 4 20 30
5 2 2 21 1 30 10 90
6 2 2 21 1 10 80 90
7 2 2 21 1 80 20 90
8 2 2 21 1 20 30 90
$EndElements
)";
}  // namespace

int main(int argc, char** argv) {
  try {
    check(argc == 2, "repository mesh path provided");
    const auto mesh = read(square);
    check(mesh.points.size() == 5 && mesh.triangles.size() == 4, "square sizes");
    check(mesh.points[4] == std::array<double, 2>{0.5, 0.5}, "coordinates preserved");
    check(mesh.triangles[1] == std::array<std::size_t, 3>{1, 2, 4}, "sparse unsorted IDs");
    check(mesh.on_boundary == std::vector<bool>({true, true, true, true, false}),
          "line endpoints are boundary; center is interior");
    double area = 0;
    for (const auto& triangle : mesh.triangles) {
      const auto& a = mesh.points[triangle[0]];
      const auto& b = mesh.points[triangle[1]];
      const auto& c = mesh.points[triangle[2]];
      area += 0.5 * std::abs((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]));
    }
    check(area == 1.0, "total area");
    check(read(change(square, "5 2 2 21 1 30 10 90", "5 2 0 30 10 90"))
          .triangles == mesh.triangles, "no element tags");
    check(read(change(square, "5 2 2 21 1 30 10 90", "5 2 4 21 1 0 99 30 10 90"))
          .triangles == mesh.triangles, "extra element tags");
    check(read(change(square, "$Elements\n8", "$Elements\n9\n9 15 0 90"))
          .on_boundary == mesh.on_boundary, "point element does not mark boundary");
    check(read(square + "$NodeData\na\n$EndNodeData\n$NodeData\nb\n$EndNodeData\n")
          .triangles == mesh.triangles, "skip repeated data sections");
    check(read(change(square, "5 2 2 21 1 30 10 90", "5 2 2 21 1 30 90 10"))
          .triangles[0] == std::array<std::size_t, 3>{0, 4, 1}, "clockwise order preserved");
    std::string crlf;
    for (const auto character : square) {
      if (character == '\n') crlf += '\r';
      crlf += character;
    }
    check(read(crlf).triangles == mesh.triangles, "CRLF file");
    const std::string minimal = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$Elements\n1\n1 2 0 30 10 20\n$EndElements\n"
        "$Nodes\n3\n30 0 0 0\n10 1 0 0\n20 0 1 0\n$EndNodes\n";
    check(read(minimal).triangles.size() == 1, "elements before nodes; no names or lines");
    check(read(minimal).on_boundary == std::vector<bool>({false, false, false}),
          "no implicit boundary inference");
    reject(change(square, "2.2 0 8", "4.1 0 8"), "only MSH 2.2 ASCII");
    reject(change(square, "2.2 0 8", "2.2 1 8"), "only MSH 2.2 ASCII");
    reject(change(square, "10 1 0 0", "30 1 0 0"), "duplicate node ID");
    reject(change(square, "30 0 0 0", "0 0 0 0"), "node ID must be positive");
    reject(change(square, "30 0 0 0", "30 0 0 1"), "xy plane");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 2 2 21 1 30 10 999"), "unknown node");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 2 2 21 1 30 10 10"), "degenerate triangle");
    reject(change(square, "90 0.5 0.5 0", "90 0.5 0 0"), "degenerate triangle");
    reject(change(square, "1 1 2 11 1 30 10", "1 1 2 11 1 30 30"), "zero-length line");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 3 2 21 1 30 10 90 20"), "unsupported element type");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 9 2 21 1 30 10 90"), "unsupported element type");
    reject(change(square, "8 2 2 21 1 20 30 90", "5 2 2 21 1 20 30 90"), "duplicate element ID");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 2 -1 30 10 90"), "negative element tag count");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 2 2 21 1 30 10"), "missing element node");
    reject(change(square, "5 2 2 21 1 30 10 90", "5 2 2 21 1 30 10 90 20"), "unexpected fields");
    reject(change(square, "$Nodes\n5", "$Nodes\n-5"), "negative record count");
    reject(change(square, "$EndNodes", "$BrokenEnd"), "missing $EndNodes");
    reject(square.substr(0, square.find("$Elements")), "missing Nodes or Elements");
    reject(square.substr(0, square.find("$EndElements")), "unexpected end of file");
    reject(square + "$Nodes\n0\n$EndNodes\n", "duplicate section");
    reject(change(minimal, "1 2 0 30 10 20", "1 15 0 30"), "no linear triangles");
    reject("", "unexpected end of file");

    // Integration check against the team's committed Gmsh mesh.
    const auto repository_mesh = fem::read_gmsh_msh2(std::string(argv[1]));
    check(repository_mesh.points.size() == 142, "repository square point count");
    check(repository_mesh.triangles.size() == 242, "repository square triangle count");
    check(std::count(repository_mesh.on_boundary.begin(), repository_mesh.on_boundary.end(), true) == 40,
          "repository square boundary count");
    std::cout << checks << " reader checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
