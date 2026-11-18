add_executable(test_assembly test_assembly.cc)
target_link_libraries(test_assembly PRIVATE fem)
add_test(NAME assembly COMMAND test_assembly)

add_executable(test_mesh_reader test_mesh_reader.cc)
target_link_libraries(test_mesh_reader PRIVATE fem)
add_test(
  NAME mesh_reader
  COMMAND test_mesh_reader
    ${CMAKE_SOURCE_DIR}/meshes/square.msh
    ${CMAKE_CURRENT_SOURCE_DIR}/data/minimal_v4.msh
)
