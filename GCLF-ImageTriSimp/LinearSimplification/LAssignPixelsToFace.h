#pragma once
#include "Basic/BoundingBox3i.h"
#include "Image/ImageDefinition.h"
#include "Mesh/LinearMeshDefinition.h"
#include "Predicates2D/Predicates2D.h"
#include "SegmentTree.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

class LinearSimplifier;

class LAssignPixelsToFace
{
  friend class LinearSimplifier;
public:
  LAssignPixelsToFace() = default;
  LAssignPixelsToFace(const LAssignPixelsToFace& rhs) :
    mesh(rhs.mesh), image(rhs.image), tree(nullptr)
  {}
  LAssignPixelsToFace& operator=(const LAssignPixelsToFace& rhs)
  {
    mesh = rhs.mesh; image = rhs.image; tree = nullptr;
    return *this;
  }

  bool operator()(
    LMeshT* _mesh, ImageT* _image,
    Eigen::MatrixXi& assigned_mask,
    const std::vector<Vec2i>& affected_pixels);
private:
  void initialize(
    LMeshT* _mesh, ImageT* _image,
    Eigen::MatrixXi& assigned_mask,
    const std::vector<Vec2i>& affected_pixels);

  void assign_pixels(
    Eigen::MatrixXi& assigned_mask,
    const Vec2d& ray_start);

  BoundingBox3i box_for_mesh();

  enum class RaySegmentInterType
  {
    NoIntersect,
    Parallel,
    Intersect,
  };

  void ray_segment_intersect(
    const Vec2d& ray_start, const Segment& seg,
    RaySegmentInterType& inter_type, double& inter_params
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
    const std::vector<IndexedSegment*>& possible_segs);

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
  LMeshT* mesh;
  ImageT* image;
  std::unique_ptr<SegmentTree> tree;
};

}// namespace ImageTriSimp
}// namespace GCLF