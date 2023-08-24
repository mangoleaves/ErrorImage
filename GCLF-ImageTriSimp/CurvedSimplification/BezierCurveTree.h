#pragma once

#include "AABB/AABBTraits.h"
#include "Curved/BezierCurve.h"
#include "Mesh/BezierMeshDefinition.h"

namespace GCLF
{
using namespace Geometry;
namespace ImageTriSimp
{

class BezierCurveTree;

class BezierCurveSplitPred
{
private:
  size_t split_dim;
public:
  BezierCurveSplitPred() :split_dim(0) {}
  BezierCurveSplitPred(size_t sd) :split_dim(sd) {}

  inline bool operator()(const BezierCurveImpl& lhs, const BezierCurveImpl& rhs)
  {
    return lhs.ctrl_points[0][split_dim] < rhs.ctrl_points[0][split_dim];
  }
};

class BezierCurveCalcBox
{
public:
  BoundingBox operator()(const BezierCurveImpl& curve)
  {
    BoundingBox box(curve.ctrl_points);
    return box;
  }
};

class BezierCurveTreeKernel
{
public:
  typedef IndexedBezierCurveImpl Primitive;
  typedef BezierCurveSplitPred SplitPred;
  typedef BezierCurveCalcBox CalcBox;
};

// TODO: complete the ray class?

class RayBezierCurveInterTraits
{
private:
  typedef IndexedBezierCurveImpl* CurvePtr;
private:
  // we only consider horizontal ray. ray_dir is (1, 0).
  Vec2d ray_start;
  std::vector<CurvePtr> intersected_curves;
  BezierCurveTree* tree;
public:
  RayBezierCurveInterTraits(BezierCurveTree* _tree, const Vec2d& _ray_start)
    :tree(_tree), ray_start(_ray_start)
  {}

  bool intersection(IndexedBezierCurveImpl& curve);

  bool do_inter(const BoundingBox& bbox) const
  {
    return bbox.max().x() >= ray_start.x() &&
      bbox.min().y() <= ray_start.y() &&
      ray_start.y() <= bbox.max().y();
  }

  inline TraversalSequence which_traversal_first(const BoundingBox& left, const BoundingBox& right)
  {
    return TraversalSequence::I_DONT_KNOW;
  }

  std::vector<CurvePtr> result()
  {
    return intersected_curves;
  }
};

/// Closest point search and intersection check.
class BezierCurveTree : public AABBTree<BezierCurveTreeKernel>
{
public:
  typedef BezierCurveTreeKernel K;
  typedef K::Primitive Primitive;
  typedef std::vector<Primitive> Curves;
  typedef Curves::iterator CurveIter;
  typedef Primitive* CurvePtr;
public:
  BezierCurveTree() = default;
  BezierCurveTree(BMeshT& mesh);

  std::vector<CurvePtr> intersect_curves(const Vec2d& ray_start);
public:
  // EdgeHandle.idx to Primitive.idx.
  std::vector<Index> global_to_local_idx;
  // We don't store tight box for primitives in tree nodes.
  // So construct tight boxes here.
  // found by Primitive.idx.
  std::vector<BoundingBox> tight_box;
};
}// namespace ImageTriSimp
}// namespace GCLF