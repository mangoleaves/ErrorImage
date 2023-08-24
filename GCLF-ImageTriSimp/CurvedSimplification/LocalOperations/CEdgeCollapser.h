#pragma once
#include "CLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class CEdgeCollapser : public CLocalOperation
{
public:
  CEdgeCollapser() = delete;
  CEdgeCollapser(CurvedSimplifier* _simplifier);

  bool init(HalfedgeHandle _hh, size_t Np);
  void clear();

  uint32_t valence_after_collapse();
  VertexHandle collapse_vertex();

  Vec3d find_smooth_target(const Vec3d& original_point)const;

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

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Points& new_ctrlpnts);
};


}// namespace ImageTriSimp
}// namespace GCLF