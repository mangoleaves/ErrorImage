#pragma once

#include "Vector2d.h"
#include "Vector2i.h"
#include "Vector2u.h"
#include "Vector3d.h"
#include "Vector3i.h"
#include "Vector3u.h"

namespace GCLF
{
namespace Geometry
{

class Segment
{
public:
  Vec3d first;
  Vec3d second;

  Segment() {}
  Segment(const Vec3d& p0, const Vec3d& p1) :first(p0), second(p1) {}
};

class Sphere
{
public:
  Vec3d center;
  double squared_radius;

  Sphere() {}
  Sphere(const Vec3d& c, const double sr) :center(c), squared_radius(sr) {}
};

class ImplicitPlane
{
public:
  Vec3d point;
  Vec3d normal;

  ImplicitPlane() {}
  ImplicitPlane(const Vec3d& p, const Vec3d& n) :point(p), normal(n) {}
};
}// namespace Geometry
}// namespace GCLF