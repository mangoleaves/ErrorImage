#pragma once
#include "CLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class CEdgeRelocater : public CLocalOperation
{
public:
  CEdgeRelocater() = delete;
  CEdgeRelocater(CurvedSimplifier* _simplifier);

  bool init(EdgeHandle _e, size_t Np);
  void clear();

  virtual std::vector<Points> get_init_population(
    std::default_random_engine& init_generator,
    std::uniform_real_distribution<double>& init_dis,
    size_t Np)const;
private:
  /***** data that describe topology in the local region *****/
  /* NOTE: these data are **SAFE** for parallel. */

  /// relocate control point on edge
  EdgeHandle relocate_eh;
  HalfedgeHandle relocate_hh, relocate_hh_opp;
  VertexHandle from_vh, to_vh;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Points& new_ctrlpnts);
};
}// namespace ImageTriSimp
}// namespace GCLF