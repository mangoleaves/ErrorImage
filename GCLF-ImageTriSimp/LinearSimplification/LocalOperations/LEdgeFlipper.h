#pragma once
#include "LLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{
class LEdgeFlipper : public LLocalOperation
{
public:
  LEdgeFlipper(LinearSimplifier* _simplifier);

  bool init(EdgeHandle _e);
  void clear();

  bool check_constraints()const;
  bool flip_would_decrease_valence();
  bool flip_would_exceed_max_valence(uint32_t max_valence);

  double error_after();
  double quality_after();

  void just_do_it();
private:
  /***** data that describe topology around the edge *****/
  /* NOTE: these data are **SAFE** for parallely finding optimal target point for relocation. */

  /// edge to flip
  EdgeHandle flip_eh;
  /// Useful handles for flip
  struct UsefulHandles
  {
    HalfedgeHandle a0, b0, a1, a2, b1, b2;
    VertexHandle   va0, va1, vb0, vb1;
    FaceHandle     fa, fb;
  }ghs, lhs;  // global handles, local handles

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update global mesh *****/

  void update_global_mesh();
};
}// namespace ImageTriSimp
}// namespace GCLF