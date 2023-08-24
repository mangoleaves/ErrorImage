cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

if (MSVC)
  add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
endif()

add_subdirectory(GCLF-ThirdPartyLibs)
add_subdirectory(GCLF-Utils)
add_subdirectory(GCLF-Geometry)
add_subdirectory(GCLF-ImageTriSimp)