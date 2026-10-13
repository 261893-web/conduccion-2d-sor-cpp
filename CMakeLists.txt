cmake_minimum_required(VERSION 3.16)
project(conduccion2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_library(laplace2d src/laplace2d.cpp)
target_include_directories(laplace2d PUBLIC include)

add_executable(conduccion2d src/main.cpp)
target_link_libraries(conduccion2d PRIVATE laplace2d)

enable_testing()
add_executable(pruebas tests/test_laplace2d.cpp)
target_link_libraries(pruebas PRIVATE laplace2d)
add_test(NAME pruebas COMMAND pruebas)
