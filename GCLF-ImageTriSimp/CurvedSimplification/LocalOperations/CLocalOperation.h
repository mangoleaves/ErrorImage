#pragma once

#include <random>
#include "Mesh/LinearMeshDefinition.h"
#include "Mesh/BezierMeshDefinition.h"
#include "CurvedSimplification/CTopoOperations.h"
#include "CurvedSimplification/CCheckConstraints.h"

namespace GCLF
{
namespace ImageTriSimp
{

#define APPROX_CURVATURE_MAX

// forward declaration
class CurvedSimplifier;

class CLocalOperation
{
public:
  CLocalOperation(CurvedSimplifier* _simplifier);

  void clear();

  virtual BoundingBox bounding_box()const;
  virtual Points get_variables()const;
  virtual std::vector<Points> get_init_population(
    std::default_random_engine& init_generator,
    std::uniform_real_distribution<double>& init_dis,
    size_t Np)const = 0;

  virtual double error_before()const;
  virtual double quality_before()const;

  virtual bool update_variables(const Points& new_ctrlpnts, size_t pidx);

  virtual double error_after(size_t pidx);
  virtual double quality_after(size_t pidx);

  virtual bool is_quality_bounded(size_t pidx, double curvature_bound, double min_angle_bound);

  virtual void just_do_it(const Points&);
protected:
  /***** input data *****/
  CurvedSimplifier* simplifier;
  ImageT* image;
  BMeshT* mesh;

  /***** data that describe topology in the local region *****/

  /* NOTE: below data are **SAFE** for parallel. */

  bool initialized;
  /// faces of the local region.
  std::set<FaceHandle> global_affected_faces;
  /// pixels on one_ring_faces
  std::vector<Vec2i> affected_pixel_coords;
  /// local mesh used for virtually updating.
  std::unique_ptr<BMeshT> local_mesh;
  /// map global vertices to local vertices.
  std::map<VertexHandle, VertexHandle> v2v;
  /// map global faces to local faces.
  std::map<FaceHandle, FaceHandle> f2f;
  /// map local vertices to global vertices (reversed v2v).
  std::vector<VertexHandle> rv2v;
  /// map local faces to global faces (reversed f2f).
  std::vector<FaceHandle> rf2f;
  /// control points to optimize.
  std::vector<Index> local_optimizing_cpi;
  /// affected local edges.
  std::vector<EdgeHandle> local_affected_edges;
  /// pixel assigned (only pixels on one_ring_faces are unassigned.)
  Eigen::MatrixXi pixel_assigned;

  /* NOTE: below data are **NOT SAFE** for parallely finding optimal target point for relocation. */

  /// data used to virtually update
  struct virtual_update_struct
  {
    Eigen::MatrixXi vt_pixel_assigned;
    BMeshT vt_local_mesh;
    bool vt_updated;
  };
  std::vector<virtual_update_struct> vus;

  /***** functions related to above data *****/

  virtual void find_affected_faces_pixels() = 0;
  virtual void initialize_local_mesh() = 0;

  /***** functions that update local mesh *****/

  virtual void update_local_ctrlpnts(BMeshT& new_local_mesh, const Points& new_ctrlpnts);
  virtual void update_pixel_color_error(BMeshT& local_mesh);

  /***** functions that update global mesh *****/

  virtual void update_global_mesh(const Points& new_ctrlpnts) = 0;

  /***** functions that check constraints *****/

  virtual bool check_constraints(const BMeshT& new_local_mesh)const;
  virtual bool operation_would_cause_flip(const BMeshT& new_local_mesh)const;

  /***** auxiliary functions *****/

  bool is_corner_vertex(VertexHandle vh);
};

}// namespace ImageTriSimp
}// namespace GCLF