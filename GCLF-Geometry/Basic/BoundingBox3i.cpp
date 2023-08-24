#include "BoundingBox3i.h"

namespace GCLF
{
namespace Geometry
{
BoundingBox3i BoundingBox3i::operator+(const BoundingBox3i& b) const
{
  BoundingBox3i result = *this;
  result.minBound.minimize(b.minBound);
  result.maxBound.maximize(b.maxBound);
  return result;
}

BoundingBox3i& BoundingBox3i::operator+=(const BoundingBox3i& b)
{
  minBound.minimize(b.minBound);
  maxBound.maximize(b.maxBound);
  return *this;
}

BoundingBox3i BoundingBox3i::operator+(const Vec3i& p)const
{
  BoundingBox3i result = *this;
  result.minBound.minimize(p);
  result.maxBound.maximize(p);
  return result;
}

BoundingBox3i& BoundingBox3i::operator+=(const Vec3i& p)
{
  minBound.minimize(p);
  maxBound.maximize(p);
  return *this;
}

void BoundingBox3i::enlarge(const Vec3i& offset)
{
  minBound -= offset;
  maxBound += offset;
}

size_t BoundingBox3i::longest_axis() const
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

bool BoundingBox3i::do_intersect(const BoundingBox3i& box)const
{
  return minBound.all_leq(box.maxBound) && box.minBound.all_leq(maxBound);
}

void BoundingBox3i::intersect(const BoundingBox3i& box)
{
  minBound.maximize(box.minBound);
  maxBound.minimize(box.maxBound);
  // won't check validity.
}
}// namespace Geometry
}// namespace GCLF