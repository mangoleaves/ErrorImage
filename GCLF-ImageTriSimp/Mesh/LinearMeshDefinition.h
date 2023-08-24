#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "Image/ImageDefinition.h"
#include "Basic/Types.h"
#include "Logger/Logger.h"

namespace GCLF
{
namespace ImageTriSimp
{
namespace OM = OpenMesh;
using namespace Utils;
using namespace Geometry;

struct MeshTraits : public OpenMesh::DefaultTraits
{
  typedef Geometry::Vec3d Point;
  typedef Geometry::Vec3d Normal;

  VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  EdgeAttributes(OpenMesh::Attributes::Status);
  HalfedgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::PrevHalfedge);

  FaceTraits
  {
  public:
    std::vector<Vec2i> pixels;
    FaceColor color;
    double error;
  };

  EdgeTraits {
  public:
    double error;
  };

  HalfedgeTraits {};

  VertexTraits {
  public:
    double error;
  };
};

// bezier mesh.
class BMeshT;

class LMeshT: public OpenMesh::TriMesh_ArrayKernelT<MeshTraits>
{
public:
  Segment seg_on_edge(EdgeHandle eh);

  Vec3d tangent_vec(HalfedgeHandle hh);

  void local_mesh(
    const std::set<FaceHandle>& local_faces, LMeshT& local,
    std::map<VertexHandle, VertexHandle>& v2v, std::map<FaceHandle, FaceHandle>& f2f,
    std::vector<VertexHandle>& rv2v, std::vector<FaceHandle>& rf2f)const;

  void scale_mesh(const ImageT& image);
};

using OpenMesh::VertexHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::FaceHandle;

using OpenMesh::SmartVertexHandle;
using OpenMesh::SmartEdgeHandle;
using OpenMesh::SmartHalfedgeHandle;
using OpenMesh::SmartFaceHandle;
}// namespace ImageTriSimp
}// namespace GCLF