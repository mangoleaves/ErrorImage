#include "LCheckConstraints.h"
#include "Mesh/OMUtils.h"

namespace GCLF
{
namespace ImageTriSimp
{

bool check_flip(
  const LMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const Vec3d& new_point)
{
  for (HalfedgeHandle h : halfedges)
  {
    // (1) find other two points to construct new triangle
    const Vec3d& from_point = mesh->point(mesh->from_vertex_handle(h));
    const Vec3d& to_point = mesh->point(mesh->to_vertex_handle(h));
    // (2) check
    if (orient2d_wrapped(from_point, to_point, new_point) <= 0)
    {
      // orientation == 0, degenerate.
      // orientation > 0, flip.
      return true;
    }
  }
  return false;
}

bool check_flip(const LMeshT* mesh)
{
  for (FaceHandle fh : mesh->faces())
  {
    if (mesh->status(fh).deleted())
      continue;
    // (1) find points
    auto [v0, v1, v2] = face_vertices(*mesh, fh);
    const Vec3d& p0 = mesh->point(v0);
    const Vec3d& p1 = mesh->point(v1);
    const Vec3d& p2 = mesh->point(v2);
    // (2) check
    if (orient2d_wrapped(p0, p1, p2) <= 0)
    {
      // orientation == 0, degenerate.
      // orientation > 0, flip.
      return true;
    }
  }
  return false;
}

bool check_small_large_angle(
  LMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const Vec3d& new_point,
  double angle_threshold
)
{
  for (HalfedgeHandle h : halfedges)
  {
    // (1) find other two points to construct new triangle
    Vec3d& from_point = mesh->point(mesh->from_vertex_handle(h));
    Vec3d& to_point = mesh->point(mesh->to_vertex_handle(h));
    // (2) calculate cos(sector angle) 
    double la2 = (from_point - to_point).sqrnorm();
    double lb2 = (to_point - new_point).sqrnorm();
    double lc2 = (from_point - new_point).sqrnorm();

    double la = sqrt(la2), lb = sqrt(lb2), lc = sqrt(lc2);

    double sa = (la2 + lb2 - lc2) / (2. * la * lb);
    double sb = (lb2 + lc2 - la2) / (2. * lb * lc);
    double sc = (la2 + lc2 - lb2) / (2. * la * lc);
    // (3) check, if cos(angle) > threshold, the angle is a small angle or large angle.
    if (sa > angle_threshold || sb > angle_threshold || sc > angle_threshold)
      return true;
  }
  return false;
}

double calc_max_cos(LMeshT* mesh, FaceHandle fh)
{
  // (1) find points
  auto [v0, v1, v2] = face_vertices(*mesh, fh);
  const Vec3d& p0 = mesh->point(v0);
  const Vec3d& p1 = mesh->point(v1);
  const Vec3d& p2 = mesh->point(v2);
  // (2) calculate cos(sector angle) 
  double la2 = (p1 - p0).sqrnorm();
  double lb2 = (p2 - p1).sqrnorm();
  double lc2 = (p0 - p2).sqrnorm();

  double la = sqrt(la2), lb = sqrt(lb2), lc = sqrt(lc2);

  double sa = (la2 + lb2 - lc2) / (2. * la * lb);
  double sb = (lb2 + lc2 - la2) / (2. * lb * lc);
  double sc = (la2 + lc2 - lb2) / (2. * la * lc);
  // (3) find max cos(angle)
  double face_max_cos = sa > sb ? (sa > sc ? sa : sc) : (sb > sc ? sb : sc);
  return face_max_cos;
}

double calc_max_cos(LMeshT* mesh)
{
  double max_cos = -1.0;
  for (FaceHandle fh : mesh->faces())
  {
    if (mesh->status(fh).deleted())
      continue;
    double face_max_cos = calc_max_cos(mesh, fh);
    max_cos = face_max_cos > max_cos ? face_max_cos : max_cos;
  }
  return max_cos;
}

}// namespace ImageTriSimp
}// namespace GCLF