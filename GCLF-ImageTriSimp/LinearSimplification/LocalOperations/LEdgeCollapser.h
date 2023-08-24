#pragma once
#include "LLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{
class LEdgeCollapser : public LLocalOperation
{
public:
  LEdgeCollapser(LinearSimplifier* _simplifier);

  bool init(EdgeHandle _e, size_t Np);
  void clear();

  uint32_t valence_after_collapse();
  VertexHandle collapse_vertex();
private:
  /***** data that describe topology around the edge *****/
  /* NOTE: these data are **SAFE** for parallely finding optimal target point for relocation. */

  /// collapse edge
  HalfedgeHandle collapse_he;
  HalfedgeHandle collapse_he_opp;
  VertexHandle from_vh, to_vh;
  /// halfedges(with respect to faces) opposite to vertex.
  std::vector<HalfedgeHandle> halfedges;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Vec3d& new_point);
};
}// namespace ImageTriSimp
}// namespace GCLF