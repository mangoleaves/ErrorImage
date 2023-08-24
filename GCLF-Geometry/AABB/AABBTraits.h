#pragma once

#include "Basic/Types.h"
#include "AABBTree.h"

namespace GCLF
{
namespace Geometry
{
class ProjectionTraits
{
private:
  typedef HeavyTriangle* HeavyTriPtr;
private:
  Vec3d m_query;
  Vec3d m_closest_point;
  double m_square_distance;
  HeavyTriPtr m_closest_triangle;
public:
  ProjectionTraits(const Vec3d& query, const Vec3d& hint, HeavyTriPtr hint_triangle)
    :m_query(query),
    m_closest_point(hint),
    m_closest_triangle(hint_triangle),
    m_square_distance((query - hint).squaredNorm())
  {}

  inline bool intersection(HeavyTriangle& tri)
  {
    Vec3d new_closest_point;
    double new_square_distance = tri.closest_point(m_query, new_closest_point);
    if (new_square_distance < m_square_distance)
    {
      m_closest_point = new_closest_point;
      m_square_distance = new_square_distance;
      m_closest_triangle = &tri;
    }
    return true;
  }

  inline bool do_inter(const BoundingBox& bbox) const
  {
    return bbox.do_intersect(Sphere(m_query, m_square_distance));
  }

  inline TraversalSequence which_traversal_first(const BoundingBox& left, const BoundingBox& right)
  {
    Vec3d dists;
    double left_dis = left.distance_to_point(m_query, dists);
    double right_dis = right.distance_to_point(m_query, dists);
    bool left_inter = left_dis < m_square_distance;
    bool right_inter = right_dis < m_square_distance;
    if (left_inter && right_inter)
      return left_dis < right_dis ? TraversalSequence::LEFT_THEN_RIGHT : TraversalSequence::RIGHT_THEN_LEFT;
    else if (left_inter)
      return TraversalSequence::ONLY_LEFT;
    else if (right_inter)
      return TraversalSequence::ONLY_RIGHT;
    else
      return TraversalSequence::NONE;
  }

  inline const Vec3d& closest_point() const { return m_closest_point; }
  inline HeavyTriPtr primitive()const { return m_closest_triangle; }
  inline double square_distance() const { return m_square_distance; }
};

template<typename Primitive>
class BoxInterTraits
{
private:
  typedef IndexedTriangle* TriPtr;
  typedef std::vector<int> Indices;
private:
  Primitive m_query;
  BoundingBox m_box_of_query;
  Indices m_indices;
public:
  BoxInterTraits(const Primitive& query)
    :m_query(query)
  {
    m_box_of_query = BoundingBox(m_query);
  }

  inline bool intersection(IndexedTriangle& tri) { m_indices.push_back(tri.index); return true; }

  inline bool do_inter(const BoundingBox& bbox) const
  {
    if (bbox.do_intersect(m_box_of_query))
      return bbox.do_intersect(m_query);
    else
      return false;
  }
  inline TraversalSequence which_traversal_first(const BoundingBox& left, const BoundingBox& right)
  {
    return TraversalSequence::I_DONT_KNOW;
  }

  inline const Indices& result() const { return m_indices; }
};

}// namespace Geometry
}// namespace GCLF