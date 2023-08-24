#pragma once

#include "SimpleTypes.h"

namespace GCLF
{
namespace Geometry
{
class Triangle
{
public:
  Vec3d ver0;
  Vec3d ver1;
  Vec3d ver2;
public:
  Triangle() {}
  Triangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3);

	Vec3d& operator[](size_t i);
	const Vec3d& operator[](size_t i) const;

  inline Vec3d centroid() { return (ver0 + ver1 + ver2) / 3.0; }
  inline Vec3d centroid() const { return (ver0 + ver1 + ver2) / 3.0; }
};

class IndexedTriangle : public Triangle
{
public:
  unsigned int index;

  IndexedTriangle() {}
  IndexedTriangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3)
  {
    ver0 = p1; ver1 = p2; ver2 = p3;
  }
  IndexedTriangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3, unsigned int fi)
    :IndexedTriangle(p1, p2, p3)
  {
    index = fi;
  }
};

class HeavyTriangle : public IndexedTriangle
{
public:
  // used for accelarating.
  Vec3d face_normal;
  Vec3d edge_vec01, edge_vec02, edge_vec12;
  Vec3d edge_nor01, edge_nor02, edge_nor12;
  double edge_sqrlen01, edge_sqrlen02, edge_sqrlen12;
  int8_t obtuse_angle;    // -1, 0, 1, 2
  uint8_t is_plane_degenerate;// true or false
public:
  HeavyTriangle() {}
  HeavyTriangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3);
  HeavyTriangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3, unsigned int fi)
    :HeavyTriangle(p1, p2, p3)
  {
    index = fi;
  }

  Vec3d& update_normal();
  inline Vec3d& normal() { return face_normal; }
  inline const Vec3d& normal() const { return face_normal; }

  double closest_point(const Vec3d& query, Vec3d& result) const;
  Vec3d closest_point(const Vec3d& query) const;
private:
  bool is_segment_degerate(std::pair<const Vec3d*, const Vec3d*> segment) const;
  std::pair<const Vec3d*, const Vec3d*> find_longest_segment()const;
  Vec3d project_to_segment(const std::pair<const Vec3d*, const Vec3d*> segment, const Vec3d& query)const;
};

class TriangleSplitPred
{
private:
  size_t split_dim;
public:
  TriangleSplitPred() :split_dim(0) {}
  TriangleSplitPred(size_t sd) :split_dim(sd) {}

  inline bool operator()(const Triangle& lhs, const Triangle& rhs)
  {
    return lhs.ver0[split_dim] < rhs.ver0[split_dim];
  }

  inline bool operator()(const IndexedTriangle& lhs, const IndexedTriangle& rhs)
  {
    return lhs.ver0[split_dim] < rhs.ver0[split_dim];
  }

  inline bool operator()(const HeavyTriangle& lhs, const HeavyTriangle& rhs)
  {
    return lhs.ver0[split_dim] < rhs.ver0[split_dim];
  }
};
}// namespace Geometry
}// namespace GCLF