cmake_minimum_required(VERSION 3.1)
project(mesh-point-triangle-mapper)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(LIBIGL REQUIRED QUIET)

include_directories(${LIBIGL_DIR}/include)

# Add your project files
<<<<<<< HEAD
add_executable(${PROJECT_NAME} mapper.cpp main.cpp)
target_link_libraries(${PROJECT_NAME} igl::opengl_glfw)
=======
add_executable(${PROJECT_NAME} main.cpp mapper.cpp)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw)
>>>>>>> 868e6ca65c79c604bec4076f4fc68d5c6743a22d
