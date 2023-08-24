#pragma once
#include "Basic/SimpleTypes.h"
#include "Basic/Triangle.h"
#include "indirect_predicates.h"
#include "ip_filtered_ex.h"

namespace GCLF
{
using namespace Geometry;

namespace ImageTriSimp
{

// NOTE: Currently we use Vec3d to represent a point on 2D.
// TODO: Build a geometry library based on 2D space.

int orient2d_wrapped(const Vec3d& p0, const Vec3d& p1, const Vec3d& q);

bool is_point_in_triangle_2(const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c);
bool is_point_in_triangle_2(const Vec3d& p, const Triangle& tri);

bool are_points_colinear_2(const Vec3d& a, const Vec3d& b, const Vec3d& c);

}// namespace ImageTriSimp
}// namespace GCLF