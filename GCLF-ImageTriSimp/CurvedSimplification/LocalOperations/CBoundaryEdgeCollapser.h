#pragma once
#include "CLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class CBoundaryEdgeCollapser : public CLocalOperation
{
public:
  CBoundaryEdgeCollapser() = delete;
  CBoundaryEdgeCollapser(CurvedSimplifier* _simplifier);

  bool init(EdgeHandle _e, size_t Np);
  void clear();

  uint32_t valence_after_collapse();
  VertexHandle collapse_vertex();

  virtual std::vector<Points> get_init_population(
    std::default_random_engine& init_generator,
    std::uniform_real_distribution<double>& init_dis,
    size_t Np)const;
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
  Index local_to_cpi, local_forward_cpi, local_backward_cpi;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update local mesh *****/

  void update_local_ctrlpnts(BMeshT& new_local_mesh, const Points& new_ctrlpnts);

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Points& new_ctrlpnts);
};

}// namespace ImageTriSimp
}// namespace GCLF