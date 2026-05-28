cmake_minimum_required(VERSION 3.16)

project(rainbow_mie_single_drop LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(rainbow_mie
    src/main.cpp
    src/MieTable.cpp
)

target_include_directories(rainbow_mie PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(APPLE)
    target_link_libraries(rainbow_mie
        "-framework OpenGL"
        "-framework GLUT"
    )
else()
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)

    target_link_libraries(rainbow_mie
        OpenGL::GL
        GLUT::GLUT
    )
endif()


