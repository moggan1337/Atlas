cmake_minimum_required(VERSION 3.15)
project(Atlas VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Warning flags
if(MSVC)
    add_compile_options(/W4 /WX-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter)
endif()

# Find OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Source files
set(ATLAS_SOURCES
    src/core/vec3.cpp
    src/core/quaternion.cpp
    src/core/matrix4x4.cpp
    src/core/transform.cpp
    src/physics/rigid_body/rigid_body.cpp
    src/physics/rigid_body/rigid_body_world.cpp
    src/physics/soft_body/soft_body.cpp
    src/physics/fluid/sph_solver.cpp
    src/physics/fluid/navier_stokes.cpp
    src/physics/cloth/cloth.cpp
    src/physics/hair/hair.cpp
    src/physics/terrain/terrain.cpp
    src/physics/raycasting/raycasting.cpp
    src/physics/constraints/constraints.cpp
)

set(ATLAS_HEADERS
    src/atlas.hpp
    src/core/vec3.hpp
    src/core/quaternion.hpp
    src/core/matrix4x4.hpp
    src/core/transform.hpp
    src/physics/rigid_body/rigid_body.hpp
    src/physics/rigid_body/rigid_body_world.hpp
    src/physics/soft_body/soft_body.hpp
    src/physics/fluid/sph_solver.hpp
    src/physics/fluid/navier_stokes.hpp
    src/physics/cloth/cloth.hpp
    src/physics/hair/hair.hpp
    src/physics/terrain/terrain.hpp
    src/physics/raycasting/raycasting.hpp
    src/physics/constraints/constraints.hpp
)

# Main library
add_library(Atlas STATIC ${ATLAS_SOURCES} ${ATLAS_HEADERS})
target_include_directories(Atlas PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Set optimization flags
target_compile_options(Atlas PRIVATE 
    $<$<CONFIG:Release>:-O3 -march=native -ffast-math>
)

# Unit tests executable
add_executable(atlas_tests 
    tests/test_core.cpp
    tests/test_rigid_body.cpp
    tests/test_soft_body.cpp
    tests/test_collision.cpp
)
target_link_libraries(atlas_tests PRIVATE Atlas)
target_include_directories(atlas_tests PRIVATE tests)

# Benchmark executable
add_executable(atlas_benchmarks 
    benchmarks/benchmark_main.cpp
)
target_link_libraries(atlas_benchmarks PRIVATE Atlas)
target_include_directories(atlas_benchmarks PRIVATE benchmarks)

# Example executable
add_executable(atlas_example 
    examples/example_main.cpp
)
target_link_libraries(atlas_example PRIVATE Atlas)
target_include_directories(atlas_example PRIVATE examples)

# Installation
install(TARGETS Atlas
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/atlas
)

install(DIRECTORY src/
    DESTINATION include/atlas
    FILES_MATCHING PATTERN "*.hpp"
)

# Testing
enable_testing()
add_test(NAME atlas_tests COMMAND atlas_tests)

# CTest options
set(CTEST_OUTPUT_ON_FAILURE ON)
