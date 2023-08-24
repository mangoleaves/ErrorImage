cmake_minimum_required(VERSION 3.10)

# set the project name
project(indirectPredicates)

# search all source files in current directory
file(GLOB SRCS
${CMAKE_CURRENT_SOURCE_DIR}
*.h *.cpp
)

# add the executable
add_library(${PROJECT_NAME} STATIC
    ${SRCS}
)
# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Compiler-specific options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
# use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
#	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
endif()
