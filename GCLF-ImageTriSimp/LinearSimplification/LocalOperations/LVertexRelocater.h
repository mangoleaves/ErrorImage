#pragma once
#include "LLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{
class LVertexRelocater : public LLocalOperation
{
public:
  LVertexRelocater(LinearSimplifier* _simplifier);

  bool init(VertexHandle _v, size_t Np);
  void clear();
private:
  /***** data that describe topology in the local region *****/
  /* NOTE: these data are **SAFE** for parallel. */

  /// relocate vertex
  VertexHandle relocate_vh;
  /// halfedges(with respect to faces) opposite to vertex.
  std::vector<HalfedgeHandle> halfedges;
  bool do_project_x, do_project_y;
  double x_project_to, y_project_to;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update local mesh *****/

  virtual void update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point);

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Vec3d& new_point);
  virtual void post_update(); // post update after updating global mesh
};
}// namespace ImageTriSimp
}// namespace GCLF