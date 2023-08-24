cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

# set the project name
project(GCLF-Utils)

# search all source files in current directory
file(GLOB SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/Logger/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Utils/*"
)

# add the library
add_library(${PROJECT_NAME} STATIC
  ${SRCS}
)

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
  spdlog::spdlog
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
