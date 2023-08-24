#pragma once
#include "LLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class LBoundaryEdgeCollapser : public LLocalOperation
{
public:
  LBoundaryEdgeCollapser(LinearSimplifier* _simplifier);

  bool init(EdgeHandle _e, size_t Np);
  void clear();

  uint32_t valence_after_collapse();
  bool is_collapsed_to_corner();
  Vec3d corner_position();
  VertexHandle collapse_vertex();

  bool check_constraints(const LMeshT& new_local_mesh);
  bool check_constraints();

  double error_after(size_t pidx);
  double error_after();
  double quality_after(size_t pidx);
  double quality_after();

  void just_do_it(const Vec3d& new_point);
private:
  /***** data that describe topology in the local region *****/
  /* NOTE: these data are **SAFE** for parallel. */

  /// edge to collapse
  HalfedgeHandle collapse_he;
  HalfedgeHandle collapse_he_opp;
  /// vertices on the edge to collapse
  VertexHandle from_vh, to_vh;
  size_t n_faces_after_collapsing;
  /// is collapse vertex a corner vertex?
  /// if it is, we should fix it.
  bool to_vh_is_corner;
  /// After collapsing, there are two boundary vertices
  /// adjacent to the collapse_vh.
  /// We name them forward_end and backward_end depending on
  /// the orientation of halfedges.
  VertexHandle forward_end_vh, backward_end_vh;
  Vec3d forward_end_p, backward_end_p;
  /// project
  bool do_project_x, do_project_y;
  double x_project_to, y_project_to;
  /// Local Mesh
  VertexHandle local_from_vh, local_to_vh;
  VertexHandle local_forward_end_vh;
  VertexHandle local_backward_end_vh;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update local mesh *****/

  virtual void update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point);

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Vec3d& new_point);
};

}// namespace ImageTriSimp
}// namespace GCLF