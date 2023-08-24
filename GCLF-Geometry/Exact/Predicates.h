#pragma once

#include "Basic/Types.h"
#include "Basic/Triangle.h"
#include "indirect_predicates.h"
#include "ip_filtered_ex.h"

namespace GCLF
{
namespace Geometry
{
/*********************************/
/********* Orientation ***********/
/*********************************/

int orient3d_wrapped(const Vec3d& p, const Vec3d& q, const Vec3d& r, const Vec3d& query);

int orient3d_wrapped(const Triangle& tri, const Vec3d& point);

bool are_points_colinear(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2);

int coplanar_orient2d(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2);

int collinear_are_ordered_along_line(const Vec3d& p, const Vec3d& q, const Vec3d& r);

/*********************************/
/************* Sign **************/
/*********************************/

inline int dot_sign(const Vec3d vec1, const Vec3d vec2)
{
  return dot_product_sign(
    vec1.x(), vec1.y(), vec1.z(),
    vec2.x(), vec2.y(), vec2.z());
}

/*********************************/
/********** Degenerate ***********/
/*********************************/

enum class TriDeType
{
  NOT_DE,
  SEG_DE,
  POINT_DE
};

TriDeType triangle_degenerate_type(const Triangle& triangle);

/*********************************/
/*********** Accurate ************/
/*********************************/

Vec3d accurate_triangle_normal(const Triangle& triangle);

Vec3d accurate_cross_normalized(const Vec3d firstBegin, const Vec3d firstEnd, const Vec3d secondBegin, const Vec3d secondEnd);
}// namespace Geometry
}// namespace GCLF