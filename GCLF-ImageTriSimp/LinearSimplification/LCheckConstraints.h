#pragma once
#include "Basic/SimpleTypes.h"
#include "Predicates2D/Predicates2D.h"
#include "Mesh/LinearMeshDefinition.h"
#include "LTopoOperations.h"

namespace GCLF
{
namespace ImageTriSimp
{

/// @brief check if there is any flipped triangles.
/// @return true if there are flipped triangles.
bool check_flip(
  const LMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const Vec3d& new_point
);

/// @brief check if there is any flipped triangles.
/// @return true if there are flipped triangles.
bool check_flip(const LMeshT* mesh);

/// @brief check if there is any small or large angle.
/// @param angle_threshold cos(small_angle) = cos(\pi - small_angle),
/// if angle < small_angle or angle > \pi - small_angle,
/// then cos(angle) > cos(small_angle) = cos(\pi - small_angle).
/// @return true if there are small or large angles.
bool check_small_large_angle(
  LMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const Vec3d& new_point,
  double angle_threshold = 0.866
);

double calc_max_cos(LMeshT* mesh, FaceHandle fh);

double calc_max_cos(LMeshT* mesh);

}// namespace ImageTriSimp
}// namespace GCLF