#include "BezierFaceDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{

BezierFace::BezierFace(const BezierFace& bface)
{
  degree = bface.degree;
  vert_idx = bface.vert_idx;
  ctrl_pnt_idx = bface.ctrl_pnt_idx;
  ctrl_pnt = bface.ctrl_pnt;
}
BezierFace::BezierFace(BezierFace&& bface)noexcept
{
  degree = bface.degree;
  vert_idx = bface.vert_idx;
  ctrl_pnt_idx = std::move(bface.ctrl_pnt_idx);
  ctrl_pnt = std::move(bface.ctrl_pnt);
}
BezierFace& BezierFace::operator=(const BezierFace& bface)
{
  degree = bface.degree;
  vert_idx = bface.vert_idx;
  ctrl_pnt_idx = bface.ctrl_pnt_idx;
  ctrl_pnt = bface.ctrl_pnt;
  return *this;
}
BezierFace& BezierFace::operator=(BezierFace&& bface)noexcept
{
  degree = bface.degree;
  vert_idx = bface.vert_idx;
  ctrl_pnt_idx = std::move(bface.ctrl_pnt_idx);
  ctrl_pnt = std::move(bface.ctrl_pnt);
  return *this;
}

/// @param v VertexHandle.idx()
/// @return The v is 0th, 1st or 2nd conner vertex of this triangle.
Index BezierFace::find_local_vert_idx(Index v) const
{
  if (vert_idx[0] == v) return 0;
  if (vert_idx[1] == v) return 1;
  if (vert_idx[2] == v) return 2;

  ASSERT(false, "can't find the expected vertex.");
  return -1;
}

Index BezierFace::find_local_edge_idx(Index v_from, Index v_to) const
{
  if (vert_idx[0] == v_from && vert_idx[1] == v_to)
    return 0;
  if (vert_idx[1] == v_from && vert_idx[2] == v_to)
    return 1;
  if (vert_idx[2] == v_from && vert_idx[0] == v_to)
    return 2;

  ASSERT(false, "can't find the expected edge.");
  return -1;
}

}// namespace ImageTriSimp
}// namespace GCLF