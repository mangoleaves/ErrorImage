# OpenMesh
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/OpenMesh)
link_libraries(${OPENMESH_LIBRARIES})
# Indirect Predicates
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IndirectPredicates)
# fmt
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fmt)
# spdlog
set(SPDLOG_FMT_EXTERNAL ON)
set(SPDLOG_FMT_EXTERNAL_HO OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/spdlog)
# HLBFGS
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/HLBFGS)
# Eigen
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.4.0 PARENT_SCOPE)