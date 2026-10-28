#pragma once
#include <fstream>
#include <stdexcept>
#include <string>

#include "fem/mesh.hh"
#include "hdnum.hh"

namespace fem {

// Legacy ASCII VTK with the nodal solution as point data (open in ParaView).
inline void write_vtk(const std::string& filename, const Mesh& mesh,
                      const hdnum::Vector<double>& u, const std::string& name = "u") {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Cannot write: " + filename);
  out << "# vtk DataFile Version 3.0\nfem2d\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << mesh.points.size() << " double\n";
  for (auto& p : mesh.points) out << p[0] << " " << p[1] << " 0\n";
  out << "CELLS " << mesh.triangles.size() << " " << 4 * mesh.triangles.size() << "\n";
  for (auto& t : mesh.triangles) out << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
  out << "CELL_TYPES " << mesh.triangles.size() << "\n";
  for (std::size_t i = 0; i < mesh.triangles.size(); ++i) out << "5\n";
  out << "POINT_DATA " << mesh.points.size() << "\n";
  out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
  for (std::size_t i = 0; i < mesh.points.size(); ++i) out << u[i] << "\n";
}

}  // namespace fem
