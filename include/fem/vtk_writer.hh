#pragma once
#include <fstream>
#include <stdexcept>
#include <string>

#include "fem/mesh.hh"
#include "hdnum.hh"

// This header contains the code which writes the mesh and the solution from main.cc into a .vtk-file, a file which
// ParaView can read and use to generate the desired visualization

// The following code has been generated using Copilot and the following documentation has been generated with Claude.

namespace fem {

// Legacy ASCII VTK with the nodal solution as point data (open in ParaView).
// u contains one value per point (u[i] belongs to mesh.points[i]), name is the label ParaView shows for it.
inline void write_vtk(const std::string& filename, const Mesh& mesh,
                      const hdnum::Vector<double>& u, const std::string& name = "u") {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Cannot write: " + filename);

  // file header: VTK version, title, text format and type of the grid (unstructured = arbitrary triangles)
  out << "# vtk DataFile Version 3.0\nfem2d\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  
  // 1. all points, VTK always expects 3D coordinates, so z = 0
  out << "POINTS " << mesh.points.size() << " double\n";
  for (auto& p : mesh.points) out << p[0] << " " << p[1] << " 0\n";
  
  // 2. all triangles as cells, each line is "3 i j k" (number of corners, then the point indices),
  //    so the total number of entries is 4 per triangle
  out << "CELLS " << mesh.triangles.size() << " " << 4 * mesh.triangles.size() << "\n";
  for (auto& t : mesh.triangles) out << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
  
  // 3. type of each cell, 5 is the VTK code for a triangle
  out << "CELL_TYPES " << mesh.triangles.size() << "\n";
  for (std::size_t i = 0; i < mesh.triangles.size(); ++i) out << "5\n";
  
  // 4. the solution: one value per point, ParaView interpolates linearly in between
  out << "POINT_DATA " << mesh.points.size() << "\n";
  out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
  for (std::size_t i = 0; i < mesh.points.size(); ++i) out << u[i] << "\n";
}

} 
