#pragma once
#include "LLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{
// forward declaration
class LinearSimplifier;

class LFaceSplitter : public LLocalOperation
{
public:
  LFaceSplitter(LinearSimplifier* _simplifier);

  bool init(FaceHandle _f, size_t Np);
  void clear();

  bool split_would_cause_over_valence(uint32_t max_valence);
  VertexHandle split_vertex() {return split_v;}
private:
  /***** data that describe topology around the face *****/
  /* NOTE: these data are **SAFE** for parallely finding optimal target point for relocation. */

  FaceHandle split_f, local_split_f;
  // ajdacent handles of split edge
  // notations come from source code of OpenMesh
  VertexHandle split_v, local_split_v;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Vec3d& new_point);
};
}// namespace ImageTriSimp
}// namespace GCLF