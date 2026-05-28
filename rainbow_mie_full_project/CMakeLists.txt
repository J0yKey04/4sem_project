cmake_minimum_required(VERSION 3.16)
project(rainbow_mie_full LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_executable(rainbow_mie
    src/main.cpp
    src/mie.cpp
    src/dispersion.cpp
    src/image.cpp
)
target_include_directories(rainbow_mie PRIVATE include)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rainbow_mie PRIVATE -O2 -Wall -Wextra -pedantic)
endif()
