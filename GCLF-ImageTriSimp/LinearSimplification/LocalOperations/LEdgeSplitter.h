#pragma once
#include "LLocalOperation.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{
// forward declaration
class LinearSimplifier;

class LEdgeSplitter : public LLocalOperation
{
public:
  LEdgeSplitter(LinearSimplifier* _simplifier);

  bool init(EdgeHandle _e, size_t Np);
  void clear();

  virtual std::vector<Vec3d> get_init_population(
    std::default_random_engine& init_generator,
    std::uniform_real_distribution<double>& init_dis,
    size_t Np)const;

  bool split_would_cause_over_valence(uint32_t max_valence);
  VertexHandle split_vertex() {return split_v;}
private:
  /***** data that describe topology around the edge *****/
  /* NOTE: these data are **SAFE** for parallely finding optimal target point for relocation. */

  EdgeHandle split_e, local_split_e;
  // ajdacent handles of split edge
  // notations come from source code of OpenMesh
  VertexHandle split_v, local_split_v;
  HalfedgeHandle h0, o0;
  // HalfedgeHandle o1, o2, h1, h2;
  bool do_project_x, do_project_y;
  double x_project_to, y_project_to;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels();
  virtual void initialize_local_mesh();
  
  /***** functions that update local mesh *****/

  virtual void update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point);

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Vec3d& new_point);
  virtual void post_update();
};
}// namespace ImageTriSimp
}// namespace GCLF