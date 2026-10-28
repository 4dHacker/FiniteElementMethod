cmake_minimum_required(VERSION 3.16)
project(fem2d VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/external/hdnum/hdnum.hh)
  message(FATAL_ERROR "HDNUM missing. Run: git submodule update --init --recursive")
endif()

# HDNUM is header-only; its headers include "src/..." relative to its root.
add_library(fem INTERFACE)
target_include_directories(fem INTERFACE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(fem SYSTEM INTERFACE ${CMAKE_SOURCE_DIR}/external/hdnum)
target_compile_options(fem INTERFACE -Wall -Wextra)

add_executable(fem2d src/main.cc)
target_link_libraries(fem2d PRIVATE fem)

enable_testing()
add_subdirectory(tests)
