#pragma once
#include "CLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class CEdgeFlipper : public CLocalOperation
{
public:
  CEdgeFlipper() = delete;
  CEdgeFlipper(CurvedSimplifier* _simplifier);

  bool init(EdgeHandle _e, size_t Np);
  void clear();

  bool flip_would_decrease_valence();
  bool flip_would_exceed_max_valence(uint32_t max_valence);

  virtual std::vector<Points> get_init_population(
    std::default_random_engine& init_generator,
    std::uniform_real_distribution<double>& init_dis,
    size_t Np)const;
private:
  /***** data that describe topology in the local region *****/
  /* NOTE: these data are **SAFE** for parallel. */

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

  virtual void update_global_mesh(const Points& new_ctrlpnts);
};

}// namespace ImageTriSimp
}// namespace GCLF