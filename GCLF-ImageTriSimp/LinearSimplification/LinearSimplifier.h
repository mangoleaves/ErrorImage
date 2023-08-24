#pragma once
#include <random>
#include <queue>
#include "ParamParser.h"
#include "Predicates2D/Predicates2D.h"
#include "Image/ImageDefinition.h"
#include "Mesh/LinearMeshDefinition.h"
#include "Mesh/OMUtils.h"
#include "ColorAndError/CalcColorOnFace.h"
#include "ColorAndError/CalcErrorOnFace.h"
#include "ColorAndError/CalcErrorOnMesh.h"
#include "ColorAndError/MeshToImage.h"
#include "LAssignPixelsToFace.h"
#include "LocalOperations/LVertexRelocater.h"
#include "LocalOperations/LEdgeCollapser.h"
#include "LocalOperations/LBoundaryEdgeCollapser.h"
#include "LocalOperations/LBoundaryVertexCollapser.h"
#include "LocalOperations/LEdgeFlipper.h"
#include "LocalOperations/LEdgeSplitter.h"
#include "LocalOperations/LFaceSplitter.h"
#include "DE_LOptimizer.h"

namespace GCLF
{
namespace ImageTriSimp
{

#define ONLY_SPLIT_HIGH_ERROR_FACE
#define ONLY_SPLIT_HIGH_ERROR_EDGE
#define ONLY_RELOCATE_HIGH_ERROR_VERTEX

// forward declaration.
class LinearSimplifier;

class LinearSimplifierBehaviors
{
public:
  typedef CalcColorOnFace<LinearSimplifier> LCalcColorOnFace;
  typedef CalcErrorOnFace<LinearSimplifier> LCalcErrorOnFace;
  typedef CalcErrorOnMesh<LinearSimplifier> LCalcErrorOnMesh;
  typedef MeshToImage<LinearSimplifier> LMeshToImage;
public:
  LinearSimplifierBehaviors() = delete;
  LinearSimplifierBehaviors(LinearSimplifier* _simplifier) :
    calc_color_on_face(_simplifier),
    calc_error_on_face(_simplifier),
    calc_error_on_mesh(_simplifier),
    mesh_to_image(_simplifier)
  {
    // default behaviors.
    calc_color_on_face.type = LCalcColorOnFace::Type::Constant;
    calc_error_on_face.type = LCalcErrorOnFace::Type::L2Norm;
    calc_error_on_mesh.p = 2.0;
  }

  LAssignPixelsToFace assign_pixels_to_face;
  LCalcColorOnFace calc_color_on_face;
  LCalcErrorOnFace calc_error_on_face;
  LCalcErrorOnMesh calc_error_on_mesh;
  LMeshToImage mesh_to_image;
};

class LinearSimplifier
{
public:
  struct Config
  {
    //  color type
    //  "constant", "linear", "quadratic"
    std::string color_type;
    // error type
    // "l1", "l2", "l4"
    std::string error_type;
    // error bound
    double error_bound;
    // quality bound
    double quality_bound;
    double split_quality_bound;
    // simplification
    size_t max_simplify_iter;
    size_t convergence_collapse_number;
    // collapse
    uint32_t max_valence;
  }config;

  enum class OptStrategy
  {
    EBQO,   // error bounded, quality optimized
    EOQB    // error optimized, quality bounded
  };

  LinearSimplifierBehaviors behavior;

  std::unique_ptr<ImageT> image;
  std::unique_ptr<LMeshT> mesh;

  std::unique_ptr<ImageT> out_image;
public:
  LinearSimplifier() :behavior(this) {}

  void initialize(ImageT& input_image, LMeshT& input_mesh, ParamTriangulator& param);

  void set_color_type(const std::string& color_type);
  void set_error_type(const std::string& error_type);
  void set_error_bound(double rmse);
  void set_max_valence(uint32_t mv);
  void set_quality_bound(double min_angle, double split_min_angle);

  void split_to_error_bound();

  void simplify();

  size_t relocate(OptStrategy opt);
  size_t collapse(OptStrategy opt);
  size_t flip(OptStrategy opt);
  size_t split_edges(OptStrategy opt);
  size_t split_faces(OptStrategy opt);

  size_t relocate_with_priority();
  size_t collapse_with_priority();
  size_t split_edges_with_priority();
  size_t split_faces_with_priority();

  double lp_error() { return behavior.calc_error_on_mesh.lp_error(); }
  double mean_lp_error() { return behavior.calc_error_on_mesh.mean_lp_error(); }
  double root_mean_lp_error() { return behavior.calc_error_on_mesh.root_mean_lp_error(); }

  /* Test Functions */
  std::vector<Vec3d> intersect_points(const Vec2d& ray_start);

private:
  /****** Optimization ********/
  double current_error;

  DE_LOptimizer optimizer;

  /****** Color & Error *******/

  void update_color_error();
  void initialize_error_on_verts();
  void initialize_error_on_edges();

  /******* Priority Queue for Collapse *******/
  template<typename HandleT>
  class Priority
  {
  public:
    HandleT handle;       // the handle to operate.
    double error_change;  // the change in error (positive -> increase, negtive -> decrease)
    Vec3d pos;            // the target position.
    uint32_t timestamp;   // if timestamp is older than current timestamp, this item is invalid.

    Priority(HandleT _handle, double _error_change, const Vec3d& _pos, uint32_t _timestamp) :
      handle(_handle), error_change(_error_change), pos(_pos), timestamp(_timestamp)
    {}
    bool operator<(const Priority& rhs) const
    {
      return error_change > rhs.error_change;  // minimal first.
    }
  };
  typedef std::priority_queue<Priority<VertexHandle>> VertexQueue;
  typedef std::priority_queue<Priority<EdgeHandle>> EdgeQueue;
  typedef std::priority_queue<Priority<FaceHandle>> FaceQueue;
  std::vector<uint32_t> verts_timestamp;
  std::vector<uint32_t> edges_timestamp;
  std::vector<uint32_t> faces_timestamp;
  VertexQueue verts_queue;
  EdgeQueue edges_queue;
  FaceQueue faces_queue;
  double verts_error_threshold;
  double edges_error_threshold;
  double faces_error_threshold;

  void init_verts_queue_to_relocate();
  // void find_verts_affected_by_relocate(VertexHandle relocate_vh, std::vector<VertexHandle>& verts);
  void update_verts_to_relocate(VertexHandle vh);

  void init_edges_queue_to_collapse();
  void find_edges_affected_by_collapse(VertexHandle collapse_vh, std::set<EdgeHandle>& edges);
  void update_edge_to_collapse(EdgeHandle eh);

  void init_edges_queue_to_split();
  void find_edges_affected_by_split(VertexHandle split_vh, std::vector<EdgeHandle>& edges);
  void update_edge_to_split(EdgeHandle eh);

  void init_faces_queue_to_split();
  // void find_faces_affected_by_split(VertexHandle split_vh, std::vector<FaceHandle>& faces);
  void update_face_to_split(FaceHandle fh);

  /****** Local Operators *******/

  bool relocate_vertex_ebqo(VertexHandle vh);
  bool collapse_edge_ebqo(EdgeHandle eh);
  bool collapse_boundary_edge_ebqo(EdgeHandle eh);
  bool collapse_boundary_vertex_ebqo(EdgeHandle eh);
  bool flip_edge_ebqo(EdgeHandle eh);

  bool relocate_vertex_eoqb(VertexHandle vh);
  bool collapse_edge_eoqb(EdgeHandle eh);
  bool collapse_boundary_edge_eoqb(EdgeHandle eh);
  bool collapse_boundary_vertex_eoqb(EdgeHandle eh);
  bool flip_edge_eoqb(EdgeHandle eh);
  bool split_edge_eoqb(EdgeHandle eh);
  bool split_face_eoqb(FaceHandle fh);

  LVertexRelocater new_relocater() { return LVertexRelocater(this); }
  LEdgeCollapser new_collapser() { return LEdgeCollapser(this); }
  LBoundaryVertexCollapser new_bv_collapser() { return LBoundaryVertexCollapser(this); }
  LBoundaryEdgeCollapser new_be_collapser() { return LBoundaryEdgeCollapser(this); }
  LEdgeFlipper new_flipper() { return LEdgeFlipper(this); }
  LEdgeSplitter new_edge_splitter() { return LEdgeSplitter(this); }
  LFaceSplitter new_face_splitter() { return LFaceSplitter(this); }

  // TEST Output files
  bool error_to_file;
  size_t operation_cnt;
  FILE* error_fout;
};

}// namespace ImageTriSimp
}// namespace GCLF