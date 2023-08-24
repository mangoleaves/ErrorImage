#pragma once
#include "Basic/BoundingBox3i.h"
#include "Image/ImageDefinition.h"
#include "Mesh/BezierMeshDefinition.h"
#include "Predicates2D/Predicates2D.h"
#include "BezierCurveTree.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class CurvedSimplifier;

class CAssignPixelsToFace
{
  friend class CurvedSimplifier;
public:
  CAssignPixelsToFace() = default;
  CAssignPixelsToFace(const CAssignPixelsToFace& rhs) :
    mesh(rhs.mesh), image(rhs.image), tree(nullptr)
  {}
  CAssignPixelsToFace& operator=(const CAssignPixelsToFace& rhs)
  {
    mesh = rhs.mesh; image = rhs.image; tree = nullptr;
    return *this;
  }

  bool operator()(
    BMeshT* _mesh, ImageT* _image,
    Eigen::MatrixXi& assigned_mask,
    const std::vector<Vec2i>& affected_pixels);
private:
  void initialize(
    BMeshT* _mesh, ImageT* _image,
    Eigen::MatrixXi& assigned_mask,
    const std::vector<Vec2i>& affected_pixels);

  void assign_pixels(
    Eigen::MatrixXi& assigned_mask,
    const Vec2d& ray_start);

  BoundingBox3i box_for_mesh();

  enum class RayBezierCurveInterType
  {
    NoIntersect,
    NearlyParallel,
    Tangent,
    One,
    Two
  };

  void ray_Bezier_curve_intersect(
    const Vec2d& ray_start, const BezierCurveImpl& bcurve,
    RayBezierCurveInterType& inter_type, std::vector<double>& inter_params
  );

  class InterPoint
  {
  public:
    enum class Type
    {
      ThroughEdge,
      ThroughVertex
    }type;
    Vec3d point;
    EdgeHandle eh;
    VertexHandle vh;
    /*
      If type is ThroughEdge:
        eh is the handle of the edge, vh is the closest vertex handle on the edge.
      If type is ThroughVertex:
        eh is invalidate, vh is the handle of the vertex.
    */
    InterPoint() = default;
    InterPoint(Type t, const Vec3d& p, EdgeHandle e, VertexHandle v) : type(t), point(p), eh(e), vh(v) {}
    bool operator==(const InterPoint& rhs) { return point.x() == rhs.point.x(); }
    bool operator<(const InterPoint& rhs) { return point.x() < rhs.point.x(); }
    bool operator<=(const InterPoint& rhs) { return point.x() <= rhs.point.x(); }
  };
  std::vector<InterPoint> find_intersection_points(
    const Vec2d& ray_start,
    const std::vector<IndexedBezierCurveImpl*>& possible_curves);

  void cluster_and_collapse_inter_points(std::vector<InterPoint>& inter_points);

  void enter_faces_along_the_ray(
    const int y, const std::vector<InterPoint>& inter_points, Eigen::MatrixXi& assigned_mask);

  bool enter_by_edge(
    const InterPoint& inter_point,
    FaceHandle current_fh, FaceHandle& next_fh);
  bool enter_by_vertex(
    const InterPoint& inter_point, FaceHandle& next_fh);
  void restart_enter_at_edge(
    const InterPoint& inter_point, FaceHandle& next_fh);

  ImageT::Color random_color();
private:
  BMeshT* mesh;
  ImageT* image;
  std::unique_ptr<BezierCurveTree> tree;
};

}// namespace ImageTriSimp
}// namespace GCLF