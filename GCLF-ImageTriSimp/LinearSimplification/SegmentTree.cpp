#include "SegmentTree.h"

namespace GCLF
{
namespace ImageTriSimp
{
bool RaySegmentInterTraits::intersection(IndexedSegment& seg)
{
  const BoundingBox& bbox = tree->tight_box[tree->global_to_local_idx[seg.idx]];
  if (bbox.max().x() >= ray_start.x() &&
    bbox.min().y() <= ray_start.y() &&
    ray_start.y() <= bbox.max().y())
    // the ray and seg's tight box are intersected,
    // put this seg into results.
    intersected_segs.push_back(&seg);
  return true; // continue to search more intersections.
}

SegmentTree::SegmentTree(LMeshT& mesh)
{
  Segs segs;
  segs.reserve(mesh.n_edges());
  global_to_local_idx.resize(mesh.n_edges(), -1);
  int prim_idx = 0;
  for (EdgeHandle eh : mesh.edges())
  {
    if (mesh.status(eh).deleted())
      continue;
    segs.emplace_back(mesh.seg_on_edge(eh), eh.idx());
    tight_box.emplace_back(segs.back());
    global_to_local_idx[eh.idx()] = prim_idx; prim_idx++;
  }
  insert(std::move(segs));
  build();
}

/// @brief find segs that intnersect with a ray.
/// @attention we only consider horizontal ray. ray_dir is (1, 0).
std::vector<SegmentTree::SegPtr> SegmentTree::intersect_segs(const Vec2d& ray_start)
{
  RaySegmentInterTraits trait(this, ray_start);
  traversal(trait);
  return trait.result();
}

}// namespace ImageTriSimp
}// namespace GCLF