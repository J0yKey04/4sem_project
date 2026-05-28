cmake_minimum_required(VERSION 3.16)

project(rainbow_simulation LANGUAGES CXX)

# =========================
# C++ standard
# =========================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# =========================
# Main executable
# =========================

add_executable(rainbow

    # main
    src/main.cpp

    # math
    src/math/Vec2.cpp
    src/math/Ray.cpp
    src/math/Geometry.cpp

    # physics
    src/physics/Optics.cpp
    src/physics/Dispersion.cpp

    # scene
    src/scene/Droplet.cpp
    src/scene/LightSource.cpp
    src/scene/SimulationScene.cpp

    # simulation
    src/simulation/RayPath.cpp
    src/simulation/SimulationResult.cpp
    src/simulation/GeometricSimulation.cpp

    # render
    src/render/Renderer.cpp
)

# =========================
# Include directories
# =========================

target_include_directories(rainbow PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# =========================
# Compiler warnings
# =========================

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rainbow PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

# =========================
# Tests
# =========================

enable_testing()

add_executable(test_reflection
    tests/test_reflection.cpp

    src/math/Vec2.cpp
    src/math/Ray.cpp
    src/math/Geometry.cpp
    src/physics/Optics.cpp
)

target_include_directories(test_reflection PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_test(NAME ReflectionTest COMMAND test_reflection)

# -------------------------

add_executable(test_refraction
    tests/test_refraction.cpp

    src/math/Vec2.cpp
    src/math/Ray.cpp
    src/math/Geometry.cpp
    src/physics/Optics.cpp
)

target_include_directories(test_refraction PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_test(NAME RefractionTest COMMAND test_refraction)

# -------------------------

add_executable(test_intersection
    tests/test_intersection.cpp

    src/math/Vec2.cpp
    src/math/Ray.cpp
    src/math/Geometry.cpp
)

target_include_directories(test_intersection PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_test(NAME IntersectionTest COMMAND test_intersection)

