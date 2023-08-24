#pragma once

#include <random>
#include "ParamParser.h"
#include "Predicates2D/Predicates2D.h"
#include "Image/ImageDefinition.h"
#include "Mesh/BezierMeshDefinition.h"
#include "Mesh/OMUtils.h"
#include "ColorAndError/CalcColorOnFace.h"
#include "ColorAndError/CalcErrorOnFace.h"
#include "ColorAndError/CalcErrorOnMesh.h"
#include "ColorAndError/MeshToImage.h"
#include "ColorAndError/ColorMap.h"
#include "CAssignPixelsToFace.h"
#include "LocalOperations/CVertexRelocater.h"
#include "LocalOperations/CEdgeRelocater.h"
#include "LocalOperations/CEdgeFlipper.h"
#include "LocalOperations/CEdgeCollapser.h"
#include "LocalOperations/CBoundaryVertexCollapser.h"
#include "LocalOperations/CBoundaryEdgeCollapser.h"
#include "DE_COptimizer.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

#define ERROR_BOUNDED_QUALITY_OPTIMIZED

class CurvedSimplifier;

class CurvedSimplifierBehaviors
{
public:
  typedef CalcColorOnFace<CurvedSimplifier> CCalcColorOnFace;
  typedef CalcErrorOnFace<CurvedSimplifier> CCalcErrorOnFace;
  typedef CalcErrorOnMesh<CurvedSimplifier> CCalcErrorOnMesh;
  typedef MeshToImage<CurvedSimplifier> CMeshToImage;
public:
  CurvedSimplifierBehaviors() = delete;
  CurvedSimplifierBehaviors(CurvedSimplifier* _simplifier) :
    calc_color_on_face(_simplifier),
    calc_error_on_face(_simplifier),
    calc_error_on_mesh(_simplifier),
    mesh_to_image(_simplifier)
  {
    // default behaviors.
    calc_color_on_face.type = CCalcColorOnFace::Type::Constant;
    calc_error_on_face.type = CCalcErrorOnFace::Type::L2Norm;
    calc_error_on_mesh.p = 2.0;
  }

  CAssignPixelsToFace assign_pixels_to_face;
  CCalcColorOnFace calc_color_on_face;
  CCalcErrorOnFace calc_error_on_face;
  CCalcErrorOnMesh calc_error_on_mesh;
  CMeshToImage mesh_to_image;
};

class CurvedSimplifier
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
    double min_angle;
    double max_curvature;
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

  CurvedSimplifierBehaviors behavior;

  std::unique_ptr<ImageT> image;
  std::unique_ptr<BMeshT> mesh;

  std::unique_ptr<ImageT> out_image;
public:
  CurvedSimplifier() : behavior(this) {}

  void initialize(ImageT& input_image, LMeshT& input_mesh, ParamTriangulator& param);
  void initialize(ImageT& input_image, BMeshT& input_mesh, ParamTriangulator& param);

  void set_color_type(const std::string& color_type);
  void set_error_type(const std::string& error_type);
  void set_error_bound(double rmse);
  void set_max_valence(uint32_t mv);
  void set_min_angle(double min_angle);
  void set_max_curvature(double min_angle);

  void simplify();

  size_t relocate_vertices(OptStrategy opt);
  size_t relocate_edges(OptStrategy opt);
  size_t collapse(OptStrategy opt);
  size_t flip(OptStrategy opt);

  size_t collapse_with_priority();

  double lp_error() { return behavior.calc_error_on_mesh.lp_error(); }
  double mean_lp_error() { return behavior.calc_error_on_mesh.mean_lp_error(); }
  double root_mean_lp_error() { return behavior.calc_error_on_mesh.root_mean_lp_error(); }

  /* TEST: Functions */
  std::vector<Vec3d> intersect_points(const Vec2d& ray_start);
  void visualize_error();
  void visualize_image();
private:
  /****** Optimization *******/
  double current_error;

  DE_COptimizer optimizer;

  /****** Color & Error *******/

  void update_color_error();

  /******* Priority Queue for Collapse *******/
  class EdgePriority
  {
  public:
    EdgeHandle eh;        // the edge to collapse.
    HalfedgeHandle hh;    // the halfedge to collapse.
    double error_change;  // the change in error (positive -> increase, negtive -> decrease)
    Points pos;           // the edge is collapsed to this position.
    uint32_t timestamp;   // if timestamp is older than current timestamp, this item is invalid.

    EdgePriority(EdgeHandle _eh, HalfedgeHandle _hh, double _error_change, Points&& _pos, uint32_t _timestamp) noexcept
      :eh(_eh), hh(_hh), error_change(_error_change), pos(std::move(_pos)), timestamp(_timestamp)
    {}
    EdgePriority(const EdgePriority& rhs) noexcept
      :eh(rhs.eh), hh(rhs.hh), error_change(rhs.error_change), pos(rhs.pos), timestamp(rhs.timestamp)
    {}
    EdgePriority(EdgePriority&& rhs) noexcept
      :eh(rhs.eh), hh(rhs.hh), error_change(rhs.error_change), pos(std::move(rhs.pos)), timestamp(rhs.timestamp)
    {}
    EdgePriority& operator=(const EdgePriority& rhs) noexcept
    {
      eh = rhs.eh; hh = rhs.hh; error_change = rhs.error_change;
      pos = rhs.pos; timestamp = rhs.timestamp;
      return *this;
    }
    EdgePriority& operator=(EdgePriority&& rhs) noexcept
    {
      eh = rhs.eh; hh = rhs.hh; error_change = rhs.error_change;
      pos = std::move(rhs.pos); timestamp = rhs.timestamp;
      return *this;
    }
    bool operator<(const EdgePriority& rhs) const
    {
      return error_change > rhs.error_change;  // minimal first.
    }
  };
  std::vector<uint32_t> edges_timestamp;
  std::priority_queue<EdgePriority> edges_queue;

  void init_edges_queue_to_collapse();
  void find_edges_affected_by_collapse(VertexHandle collapse_vh, std::set<EdgeHandle>& edges);
  void update_edge_in_queue(EdgeHandle eh);

  /****** Local Operators *******/

  CVertexRelocater new_v_relocater() { return CVertexRelocater(this); }
  CEdgeRelocater new_e_relocater() { return CEdgeRelocater(this); }
  CEdgeFlipper new_flipper() { return CEdgeFlipper(this); }
  CEdgeCollapser new_collapser() { return CEdgeCollapser(this); }
  CBoundaryVertexCollapser new_bv_collapser() { return CBoundaryVertexCollapser(this); }
  CBoundaryEdgeCollapser new_be_collapser() { return CBoundaryEdgeCollapser(this); }
};
}// namespace ImageTriSimp
}// namespace GCLF