#include "LinearMeshDefinition.h"
#include "OMUtils.h"
#include "LinearSimplification/LTopoOperations.h"

namespace GCLF
{
namespace ImageTriSimp
{

Segment LMeshT::seg_on_edge(EdgeHandle eh)
{
  HalfedgeHandle hh = halfedge_handle(eh, 0);
  VertexHandle from_vh = from_vertex_handle(hh);
  VertexHandle to_vh = to_vertex_handle(hh);
  if (from_vh.idx() > to_vh.idx())
    std::swap(from_vh, to_vh);
  return Segment(point(from_vh), point(to_vh));
}

Vec3d LMeshT::tangent_vec(HalfedgeHandle hh)
{
  VertexHandle from_vh = from_vertex_handle(hh);
  VertexHandle to_vh = to_vertex_handle(hh);
  return point(to_vh) - point(from_vh);
}

/// @brief Extract a local region to form a mesh.
/// @param [in] local_faces Faces contained by the local region.
/// @param [out] local The local mesh.
/// @param [out] v2v map global vertices to local vertices. 
/// @param [out] f2f map global faces to local faces.
/// @param [out] rv2v map local vertices to global vertices (reversed v2v).
/// @param [out] rf2f map local faces to global faces (reversed f2f).
void LMeshT::local_mesh(
  const std::set<FaceHandle>& local_faces, LMeshT& local,
  std::map<VertexHandle, VertexHandle>& v2v, std::map<FaceHandle, FaceHandle>& f2f,
  std::vector<VertexHandle>& rv2v, std::vector<FaceHandle>& rf2f)const
{
  // add vertices to local mesh.
  std::set<VertexHandle> local_vertices;
  vertices_around_faces(this, local_faces, local_vertices);
  rv2v.resize(local_vertices.size());
  for (VertexHandle vh : local_vertices)
  {
    VertexHandle local_vh = local.add_vertex(point(vh));
    v2v[vh] = local_vh;
    rv2v[local_vh.idx()] = vh;
  }

  // add faces to local mesh.
  rf2f.resize(local_faces.size());
  for (FaceHandle fh : local_faces)
  {
    auto [v0, v1, v2] = face_vertices(*this, fh);
    FaceHandle local_fh = local.add_face(v2v.at(v0), v2v.at(v1), v2v.at(v2));
    f2f[fh] = local_fh;
    rf2f[local_fh.idx()] = fh;
  }
}

void LMeshT::scale_mesh(const ImageT& image)
{
  // calculate the bounding box for mesh
  BoundingBox box;
  for (VertexHandle vh : vertices())
  {
    box += point(vh);
  }

  ASSERT(box.min().x() == 0. && box.min().y() == 0., "mesh does not start at (0, 0).");

  int width = image.width;
  int height = image.height;

  double mesh_width = box.max().x();
  double mesh_height = box.max().y();

  double scale_width = (double)width/ mesh_width ;
  double scale_height = (double)height / mesh_height;

  for (VertexHandle vh : vertices())
  {
    Vec3d& p = point(vh);
    p.x() *= scale_width;
    p.y() *= scale_height;
  }

  // project boundary vertices to integer.
  for (VertexHandle vh : vertices())
  {
    if (!is_boundary(vh))
      continue;
    HalfedgeHandle hh = halfedge_handle(vh);
    Vec3d& cur_p = point(vh);
    const Vec3d& next_p = point(to_vertex_handle(hh));
    const Vec3d& prev_p = point(from_vertex_handle(prev_halfedge_handle(hh)));
    bool fix_x = fabs(cur_p.x() - next_p.x()) < 1e-6 && fabs(cur_p.x() - prev_p.x()) < 1e-6;
    bool fix_y = fabs(cur_p.y() - next_p.y()) < 1e-6 && fabs(cur_p.y() - prev_p.y()) < 1e-6;
    if (!fix_x && !fix_y)
    {
      // fix it to corner
      if (fabs(cur_p.x() - 0.0) < 1e-6)
        cur_p.x() = 0.0;
      else if (fabs(cur_p.x() - width) < 1e-6)
        cur_p.x() = width;

      if (fabs(cur_p.y() - 0.0) < 1e-6)
        cur_p.y() = 0.0;
      else if (fabs(cur_p.y() - height) < 1e-6)
        cur_p.y() = height;
    }
    else if (fix_x)
    {
      if (fabs(cur_p.x() - 0.0) < 1e-6)
        cur_p.x() = 0.0;
      else if (fabs(cur_p.x() - width) < 1e-6)
        cur_p.x() = width;
    }
    else if (fix_y)
    {
      if (fabs(cur_p.y() - 0.0) < 1e-6)
        cur_p.y() = 0.0;
      else if (fabs(cur_p.y() - height) < 1e-6)
        cur_p.y() = height;
    }
    else
    {
      ASSERT(false, "fail to find boundary.");
    }
  }
}

}// namespace ImageTriSimp
}// namespace GCLF