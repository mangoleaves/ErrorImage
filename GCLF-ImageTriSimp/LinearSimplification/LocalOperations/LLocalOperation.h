#pragma once

#include <random>
#include "Image/ImageDefinition.h"
#include "Mesh/LinearMeshDefinition.h"
#include "LinearSimplification/LTopoOperations.h"
#include "LinearSimplification/LCheckConstraints.h"

namespace GCLF
{
namespace ImageTriSimp
{

// forward declaration
class LinearSimplifier;

class LLocalOperation
{
public:
  LLocalOperation(LinearSimplifier* _simplifier);

  void clear();

  virtual BoundingBox bounding_box()const;
  virtual Vec3d get_variable()const;
  virtual std::vector<Vec3d> get_init_population(
    std::default_random_engine& init_generator,
    std::uniform_real_distribution<double>& init_dis,
    size_t Np)const;

  virtual double error_before()const;
  virtual double quality_before()const;

  virtual bool update_variable(const Vec3d& new_point, size_t pidx);

  virtual double error_after(size_t pidx);
  virtual double quality_after(size_t pidx);

  virtual void just_do_it(const Vec3d&);
protected:
  /***** input data *****/
  LinearSimplifier* simplifier;
  ImageT* image;
  LMeshT* mesh;

  /***** data that describe topology in the local region *****/

  /* NOTE: below data are **SAFE** for parallel. */

  bool initialized;
  /// faces of the local region.
  std::set<FaceHandle> global_affected_faces;
  /// pixels on one_ring_faces
  std::vector<Vec2i> affected_pixel_coords;
  /// local mesh used for virtually updating.
  std::unique_ptr<LMeshT> local_mesh;
  /// map global vertices to local vertices.
  std::map<VertexHandle, VertexHandle> v2v;
  /// map global faces to local faces.
  std::map<FaceHandle, FaceHandle> f2f;
  /// map local vertices to global vertices (reversed v2v).
  std::vector<VertexHandle> rv2v;
  /// map local faces to global faces (reversed f2f).
  std::vector<FaceHandle> rf2f;
  /// control points to optimize.
  VertexHandle local_optimizing_vertex;
  /// pixel assigned (only pixels on one_ring_faces are unassigned.)
  Eigen::MatrixXi pixel_assigned;

  /* NOTE: below data are **NOT SAFE** for parallely finding optimal target point for relocation. */

  /// data used to virtually update
  struct virtual_update_struct
  {
    Eigen::MatrixXi vt_pixel_assigned;
    LMeshT vt_local_mesh;
    bool vt_updated;
  };
  std::vector<virtual_update_struct> vus;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels() = 0;
  virtual void initialize_local_mesh() = 0;

  /***** functions that update local mesh *****/

  virtual void update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point);
  virtual void update_pixel_color_error(LMeshT& local_mesh);

  /***** functions that update global mesh *****/

  virtual bool check_constraints(const LMeshT& new_local_mesh)const;
  virtual void update_global_mesh(const Vec3d& new_point) {}
  virtual void post_update() {}; // post update after updating global mesh

  /***** functions that check constraints *****/

  virtual bool operation_would_cause_flip(const LMeshT& new_local_mesh)const;

  /***** auxiliary functions *****/

  bool is_corner_vertex(VertexHandle vh);
};

}// namespace ImageTriSimp
}// namespace GCLF