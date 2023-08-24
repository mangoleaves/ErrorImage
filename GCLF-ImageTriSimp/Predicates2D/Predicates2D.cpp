#include "Predicates2D.h"

namespace GCLF
{
namespace ImageTriSimp
{

/// @brief The orientation of query point "q" with respect to segment "p0->p1"
/// @return 1: left, 0:on, -1:right.
int orient2d_wrapped(const Vec3d& p0, const Vec3d& p1, const Vec3d& q)
{
  return orient2d(p0.x(), p0.y(), p1.x(), p1.y(), q.x(), q.y());
}

/// @brief check if query point "p" is inside the triangle "abc".
/// The order of a,b and c should be couter-clockwise.
/// @return true if inside the triangle.
bool is_point_in_triangle_2(const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
  return orient2d_wrapped(a, b, p) >= 0 && orient2d_wrapped(b, c, p) >= 0 && orient2d_wrapped(c, a, p) >= 0;
}

bool is_point_in_triangle_2(const Vec3d& p, const Triangle& tri)
{
  return is_point_in_triangle_2(p, tri.ver0, tri.ver1, tri.ver2);
}

bool are_points_colinear_2(const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
  return orient2d(a.x(), a.y(), b.x(), b.y(), c.x(), c.y()) != ZERO;
}

}// namespace ImageTriSimp
}// namespace GCLF