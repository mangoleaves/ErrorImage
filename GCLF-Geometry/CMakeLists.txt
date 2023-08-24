cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 17)

# set the project name
project(GCLF-Geometry)

# search all source files in current directory
file(GLOB SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/Basic/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Exact/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Curved/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/AABB/*"
)

# OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP Not Found!")
    endif()
endif()

# add the library
add_library(${PROJECT_NAME} STATIC
  ${SRCS}
)

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${EIGEN3_INCLUDE_DIR}
)

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
  indirectPredicates
  OpenMP::OpenMP_CXX
  GCLF-Utils
)

# Compiler-specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # -frounding-math grant IEEE 754 compliance
  target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)

  # -g generate debug information
  target_compile_options(${PROJECT_NAME} PUBLIC -g)
# use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
# target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # grant IEEE 754 compliance
  target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")

  # generate debug information
  target_compile_options(${PROJECT_NAME} PUBLIC "/Zi")

  # use intrinsic functions
  target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")

  # turn off annoying warnings
  target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")

  # turn on multiprocessor compile.
  target_compile_options(${PROJECT_NAME} PUBLIC "/MP")
endif()
