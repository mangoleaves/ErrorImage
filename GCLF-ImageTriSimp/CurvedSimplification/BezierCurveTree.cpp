#include "BezierCurveTree.h"

namespace GCLF
{
namespace ImageTriSimp
{
bool RayBezierCurveInterTraits::intersection(IndexedBezierCurveImpl& curve)
{
  const BoundingBox& bbox = tree->tight_box[tree->global_to_local_idx[curve.idx]];
  if (bbox.max().x() >= ray_start.x() &&
    bbox.min().y() <= ray_start.y() &&
    ray_start.y() <= bbox.max().y())
    // the ray and curve's tight box are intersected,
    // put this curve into results.
    intersected_curves.push_back(&curve);
  return true; // continue to search more intersections.
}

BezierCurveTree::BezierCurveTree(BMeshT& mesh)
{
  Curves curves;
  curves.reserve(mesh.n_edges());
  global_to_local_idx.resize(mesh.n_edges(), -1);
  int prim_idx = 0;
  for (EdgeHandle eh : mesh.edges())
  {
    if (mesh.status(eh).deleted())
      continue;
    curves.emplace_back(mesh.curve_on_edge(eh), eh.idx());
    tight_box.emplace_back(curves.back().ctrl_points);
    global_to_local_idx[eh.idx()] = prim_idx; prim_idx++;
  }
  insert(std::move(curves));
  build();
}

/// @brief find curves that intnersect with a ray.
/// @attention we only consider horizontal ray. ray_dir is (1, 0).
std::vector<BezierCurveTree::CurvePtr> BezierCurveTree::intersect_curves(const Vec2d& ray_start)
{
  RayBezierCurveInterTraits trait(this, ray_start);
  traversal(trait);
  return trait.result();
}

}// namespace ImageTriSimp
}// namespace GCLF