#pragma once

#include "AABB/AABBTraits.h"
#include "Basic/SimpleTypes.h"
#include "Mesh/LinearMeshDefinition.h"

namespace GCLF
{
using namespace Geometry;

namespace Geometry
{

class IndexedSegment : public Segment
{
public:

  IndexedSegment() :Segment() {}
  IndexedSegment(const Segment& seg, int _idx) :
    Segment(seg), idx(_idx)
  {}
  IndexedSegment(Segment&& seg, int _idx) :
    Segment(std::move(seg)), idx(_idx)
  {}
  IndexedSegment(const IndexedSegment& rhs) noexcept
    :Segment(rhs), idx(rhs.idx)
  {}
  IndexedSegment(IndexedSegment&& rhs) noexcept
    :Segment(std::move(rhs)), idx(rhs.idx)
  {}
  IndexedSegment& operator=(const IndexedSegment& rhs) noexcept
  {
    first = rhs.first;
    second = rhs.second;
    idx = rhs.idx;
    return *this;
  }
  IndexedSegment& operator=(IndexedSegment&& rhs) noexcept
  {
    first = rhs.first;
    second = rhs.second;
    idx = rhs.idx;
    return *this;
  }

  /* Properties */
  int idx;
};

}// namespace Geometry

namespace ImageTriSimp
{

class SegmentTree;

class SegmentSplitPred
{
private:
  size_t split_dim;
public:
  SegmentSplitPred() :split_dim(0) {}
  SegmentSplitPred(size_t sd) :split_dim(sd) {}

  inline bool operator()(const Segment& lhs, const Segment& rhs)
  {
    return lhs.first[split_dim] < rhs.first[split_dim];
  }
};

class SegmentCalcBox
{
public:
  BoundingBox operator()(const Segment& seg)
  {
    return BoundingBox(seg);
  }
  BoundingBox operator()(const IndexedSegment& seg)
  {
    return BoundingBox(*static_cast<const Segment*>(&seg));
  }
};

class SegmentTreeKernel
{
public:
  typedef IndexedSegment Primitive;
  typedef SegmentSplitPred SplitPred;
  typedef SegmentCalcBox CalcBox;
};

// TODO: complete the ray class?

class RaySegmentInterTraits
{
private:
  typedef IndexedSegment* SegPtr;
private:
  // we only consider horizontal ray. ray_dir is (1, 0).
  Vec2d ray_start;
  std::vector<SegPtr> intersected_segs;
  SegmentTree* tree;
public:
  RaySegmentInterTraits(SegmentTree* _tree, const Vec2d& _ray_start)
    :tree(_tree), ray_start(_ray_start)
  {}

  bool intersection(IndexedSegment& seg);

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

  std::vector<SegPtr> result()
  {
    return intersected_segs;
  }
};

/// Closest point search and intersection check.
class SegmentTree : public AABBTree<SegmentTreeKernel>
{
public:
  typedef SegmentTreeKernel K;
  typedef K::Primitive Primitive;
  typedef std::vector<Primitive> Segs;
  typedef Segs::iterator SegIter;
  typedef Primitive* SegPtr;
public:
  SegmentTree() = default;
  SegmentTree(LMeshT& mesh);

  std::vector<SegPtr> intersect_segs(const Vec2d& ray_start);
public:
  // EdgeHandle.idx to Primitive.idx.
  std::vector<int> global_to_local_idx;
  // We don't store tight box for primitives in tree nodes.
  // So construct tight boxes here.
  // found by Primitive.idx.
  std::vector<BoundingBox> tight_box;
};
}// namespace ImageTriSimp
}// namespace GCLF