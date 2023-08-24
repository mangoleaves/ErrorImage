cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

project(GCLF-ImageTriSimp)

# OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP Not Found!")
    endif()
endif()

cmake_policy(SET CMP0074 NEW)
# Boost
# set(Boost_DEBUG ON)
set(Boost_NO_WARN_NEW_VERSIONS ON)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.80.0
    REQUIRED
    COMPONENTS filesystem json
)

# MKL (Intel oneAPI Math Kernel Library)
# find_package(MKL CONFIG)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Predicates2D SRC_Predicates2D)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/LinearSimplification SRC_LinearSimplification)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/LinearSimplification/LocalOperations
    SRC_LinearSimplification_LocalOperations)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/CurvedSimplification SRC_CurvedSimplification)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/CurvedSimplification/LocalOperations
    SRC_CurvedSimplification_LocalOperations)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Mesh SRC_Mesh)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Image SRC_Image)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/ColorAndError SRC_ColorAndError)

add_library(${PROJECT_NAME} STATIC
    ${SRC_Predicates2D}
    ${SRC_LinearSimplification}
    ${SRC_LinearSimplification_LocalOperations}
    ${SRC_CurvedSimplification}
    ${SRC_CurvedSimplification_LocalOperations}
    ${SRC_Mesh}
    ${SRC_Image}
    ${SRC_ColorAndError}
    ParamParser.cpp
)
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    # ${MKL_INCLUDE}
)
target_link_libraries(${PROJECT_NAME} PUBLIC
    GCLF-Geometry
    OpenMP::OpenMP_CXX
    HLBFGS
    ${Boost_LIBRARIES}
    # $<LINK_ONLY:MKL::MKL>
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

#QT
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui)
get_target_property(_qmake_executable Qt5::qmake IMPORTED_LOCATION)
get_filename_component(_qt_bin_dir "${_qmake_executable}" DIRECTORY)
find_program(DEPLOYQT_EXECUTABLE NAMES windeployqt macdeployqt HINTS "${_qt_bin_dir}")

add_executable(exe${PROJECT_NAME} main.cpp)
target_include_directories(exe${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(exe${PROJECT_NAME} PUBLIC
    ${PROJECT_NAME}
    Qt5::Core
    Qt5::Widgets
    Qt5::Gui
)
target_link_options(exe${PROJECT_NAME}  PUBLIC
    "/DEBUG"
)

#if(WIN32)
#    add_custom_command(TARGET exe${PROJECT_NAME}
#    POST_BUILD
#    COMMAND ${DEPLOYQT_EXECUTABLE} "$<TARGET_FILE:exe${PROJECT_NAME}>"
#        --$<LOWER_CASE:$<CONFIG>> --qmldir ${CMAKE_SOURCE_DIR}
#    )
#endif()