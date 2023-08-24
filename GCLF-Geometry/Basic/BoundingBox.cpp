#include "BoundingBox.h"

namespace GCLF
{
namespace Geometry
{

BoundingBox BoundingBox::operator+(const BoundingBox& b) const
{
  BoundingBox result = *this;
  result.minBound.minimize(b.minBound);
  result.maxBound.maximize(b.maxBound);
  return result;
}

BoundingBox& BoundingBox::operator+=(const BoundingBox& b)
{
  minBound.minimize(b.minBound);
  maxBound.maximize(b.maxBound);
  return *this;
}

BoundingBox BoundingBox::operator+(const Vec3d& p)const
{
  BoundingBox result = *this;
  result.minBound.minimize(p);
  result.maxBound.maximize(p);
  return result;
}

BoundingBox& BoundingBox::operator+=(const Vec3d& p)
{
  minBound.minimize(p);
  maxBound.maximize(p);
  return *this;
}

void BoundingBox::enlarge(const Vec3d& offset)
{
  minBound -= offset;
  maxBound += offset;
}

size_t BoundingBox::longest_axis() const
{
  const double dx = maxBound.x() - minBound.x();
  const double dy = maxBound.y() - minBound.y();
  const double dz = maxBound.z() - minBound.z();

  if (dx >= dy)
  {
    if (dx >= dz)
      return 0;
    else // dz>dx and dx>=dy
      return 2;
  }
  else // dy>dx
  {
    if (dy >= dz)
      return 1;
    else  // dz>dy and dy>dx
      return 2;
  }
}

bool BoundingBox::do_intersect(const Vec3d& point) const
{
  return minBound.all_leq(point) && point.all_leq(maxBound);
}

bool BoundingBox::do_intersect(const BoundingBox& b) const
{
  return minBound.all_leq(b.maxBound) && b.minBound.all_leq(maxBound);
}

bool BoundingBox::do_intersect(const Segment& segment) const
{
  return do_intersect(BoundingBox(segment));
}

bool BoundingBox::do_intersect(const Sphere& sphere) const
{
  auto& center = sphere.center;
  double d = 0.0;
  double distance = 0.0;

  // for x
  if (center.x() < minBound.x())
  {
    d = minBound.x() - center.x();
    distance += d * d;
  }
  else if (center.x() > maxBound.x())
  {
    d = center.x() - maxBound.x();
    distance += d * d;
  }
  if (distance > sphere.squared_radius)
    return false;
  // for y
  if (center.y() < minBound.y())
  {
    d = minBound.y() - center.y();
    distance += d * d;
  }
  else if (center.y() > maxBound.y())
  {
    d = center.y() - maxBound.y();
    distance += d * d;
  }
  if (distance > sphere.squared_radius)
    return false;
  // for z
  if (center.z() < minBound.z())
  {
    d = minBound.z() - center.z();
    distance += d * d;
  }
  else if (center.z() > maxBound.z())
  {
    d = center.z() - maxBound.z();
    distance += d * d;
  }

  return distance <= sphere.squared_radius;
}

double BoundingBox::distance_to_point(const Vec3d& query, Vec3d& dists) const
{
  double distance = 0.0;
  for (size_t i = 0;i < 3;i++)
  {
    if (query[i] < minBound[i])
      dists[i] = minBound[i] - query[i];
    else if (query[i] > maxBound[i])
      dists[i] = query[i] - maxBound[i];
    else
      dists[i] = 0.0;
  }
  return dists.squaredNorm();
}

/// @details 
/// \par We mark eight points of bounding box with indices 1 to 8.
/// Index i means the i-th point lies in i-th octant.(box is moved to origin.)
/// \par Finding whether an implicit plane intersects a box
/// is same as finding whether two corner points of the box
/// are lies in different sides of the plane.
/// \par We chose proper corner points by normal of the plane.
/// For example, if normal lies in first octant or seventh octant,
/// we chose first point and seventh point as two corner points.
/// It is same with other octants.
/*
bool BoundingBox::do_intersect(const ImplicitPlane& plane)const
{
  static int swap_dimension[8] = { -1, 0, 2, 1, 2, 1, -1, 0 };

  int octant = find_octant(plane.normal);
  int swap_dim = swap_dimension[octant];

  return do_intersect(plane, swap_dim);
}

bool BoundingBox::do_intersect(const ImplicitPlane& plane, const int swap_dim)const
{
  if (swap_dim >= 0)
  {
    Vec3d corner0(min()), corner1(max());
    std::swap(corner0[swap_dim], corner1[swap_dim]);
    int ori0 = orient3d_wrapped(plane, corner0);
    int ori1 = orient3d_wrapped(plane, corner1);
    return ori0 * ori1 != 1;
  }
  else
  {
    int ori0 = orient3d_wrapped(plane, min());
    int ori1 = orient3d_wrapped(plane, max());
    return ori0 * ori1 != 1;
  }
}
*/

/// @brief After check box of triangle doesn't intersect with box,
/// call this function to check the rest of possible cases.
bool BoundingBox::do_intersect(const Triangle& tri)const
{
#define map_triangle_2D(dim0, dim1)\
  const double& x0=tri.ver0.dim0(); const double& x1=tri.ver1.dim0(); const double& x2=tri.ver2.dim0();\
  const double& y0=tri.ver0.dim1(); const double& y1=tri.ver1.dim1(); const double& y2=tri.ver2.dim1();

#define map_box_2D(dim0, dim1)\
  const double& xmin=min().dim0(); const double& xmax=max().dim0();\
  const double& ymin=min().dim1(); const double& ymax=max().dim1(); 

  // Consider a triangle on 2d subspace.
  // We use three segments from triangle as separating line.
  // For each separating line, if it separates the third point of triangle and box on 2d subspace,
  // triangle do not intersect with box, as well as in 3d space.

  // So, we test whether 2d triangle intersects with box on x-y, y-z, x-z subspace.
  // If they do not intersect in any subspace, they do not intersect in 3d space.
  // Otherwise, triangle may intersect with box.
  // We need to do triangle-triangle intersection test.

#define separate_2D(seg_x0, seg_y0, seg_x1, seg_y1, px, py, end)\
  int seg_x_diff_sign = seg_x1 > seg_x0 ? 1 : seg_x1 == seg_x0 ? 0 : -1;\
  int seg_y_diff_sign = seg_y1 > seg_y0 ? 1 : seg_y1 == seg_y0 ? 0 : -1;\
  if (seg_x_diff_sign * seg_y_diff_sign == 1)\
  {\
    int ori = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, px, py);\
    if (ori == 0)\
    {\
      int ori0 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymax);\
      int ori1 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymin);\
      if (ori0 * ori1 == 1) return false;\
    }\
    else\
    {\
      int ori0 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymax);\
      if (ori0 * ori >= 0) goto end;\
      int ori1 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymin);\
      if (ori1 * ori >= 0) goto end;\
      return false;\
    }\
  }\
  else if (seg_x_diff_sign * seg_y_diff_sign == -1)\
  {\
    int ori = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, px, py);\
    if (ori == 0)\
    {\
      int ori0 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymin);\
      int ori1 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymax);\
      if (ori0 * ori1 == 1) return false;\
    }\
    else\
    {\
      int ori0 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymin);\
      if (ori0 * ori >= 0) goto end;\
      int ori1 = orient2d(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymax);\
      if (ori1 * ori >= 0) goto end;\
      return false;\
    }\
  }\
  else if (seg_x_diff_sign == 0 && seg_y_diff_sign != 0)\
  {\
    int ori = px > seg_x0 ? 1 : px == seg_x0 ? 0 : -1;\
    if (ori == 0)\
    {\
      int ori0 = xmin > seg_x0 ? 1 : xmin == seg_x0 ? 0 : -1;\
      int ori1 = xmax > seg_x0 ? 1 : xmax == seg_x0 ? 0 : -1;\
      if (ori0 * ori1 == 1) return false;\
    }\
    else\
    {\
      int ori0 = xmin > seg_x0 ? 1 : xmin == seg_x0 ? 0 : -1;\
      if (ori0 * ori >= 0) goto end;\
      int ori1 = xmax > seg_x0 ? 1 : xmax == seg_x0 ? 0 : -1;\
      if (ori1 * ori >= 0) goto end;\
      return false;\
    }\
  }\
  else if (seg_y_diff_sign == 0 && seg_x_diff_sign != 0)\
  {\
    int ori = py > seg_y0 ? 1 : py == seg_y0 ? 0 : -1;\
    if (ori == 0)\
    {\
      int ori0 = ymin > seg_y0 ? 1 : ymin == seg_y0 ? 0 : -1;\
      int ori1 = ymax > seg_y0 ? 1 : ymax == seg_y0 ? 0 : -1;\
      if (ori0 * ori1 == 1) return false;\
    }\
    else\
    {\
      int ori0 = ymin > seg_y0 ? 1 : ymin == seg_y0 ? 0 : -1;\
      if (ori0 * ori >= 0) goto end;\
      int ori1 = ymax > seg_y0 ? 1 : ymax == seg_y0 ? 0 : -1;\
      if (ori1 * ori >= 0) goto end;\
      return false;\
    }\
  }

#define separate_subspace(dim0, dim1)\
{\
  map_triangle_2D(dim0, dim1);\
  map_box_2D(dim0, dim1);\
  {separate_2D(x0, y0, x1, y1, x2, y2, END_##dim0##_##dim1##_01);}\
END_##dim0##_##dim1##_01:;\
  {separate_2D(x1, y1, x2, y2, x0, y0, END_##dim0##_##dim1##_12);}\
END_##dim0##_##dim1##_12:;\
  {separate_2D(x2, y2, x0, y0, x1, y1, END_##dim0##_##dim1##_20);}\
END_##dim0##_##dim1##_20:;\
}
  // we find that about 97% cases report triangle intersects with box
  // (test CGAL::do_intersect(triangle, box) too).
  // skipping test triangle-box intersection 
  // and directly doing triangle-triangle intersection
  // saves more time.
#if 0
  separate_subspace(x, y);
  separate_subspace(y, z);
  separate_subspace(x, z);
#endif
  return true;
}
}// namespace Geometry
}// namespace GCLF