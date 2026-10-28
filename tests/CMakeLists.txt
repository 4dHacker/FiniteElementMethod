add_executable(test_assembly test_assembly.cc)
target_link_libraries(test_assembly PRIVATE fem)
add_test(NAME assembly COMMAND test_assembly)
