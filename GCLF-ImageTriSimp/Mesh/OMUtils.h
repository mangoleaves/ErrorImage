#pragma once

#include <tuple>
#include <OpenMesh/Core/Mesh/TriConnectivity.hh>
#include <OpenMesh/Core/Utils/Property.hh>

namespace GCLF
{
namespace ImageTriSimp
{

namespace OM = OpenMesh;
using namespace Utils;
using namespace Geometry;


using OpenMesh::VertexHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::FaceHandle;

using OpenMesh::SmartVertexHandle;
using OpenMesh::SmartEdgeHandle;
using OpenMesh::SmartHalfedgeHandle;
using OpenMesh::SmartFaceHandle;

template<typename MeshT>
inline std::tuple<VertexHandle, VertexHandle, VertexHandle> face_vertices(MeshT& mesh, FaceHandle fh)
{
  auto fv_iter = mesh.cfv_begin(fh);
  VertexHandle v0 = *fv_iter; fv_iter++;
  VertexHandle v1 = *fv_iter; fv_iter++;
  VertexHandle v2 = *fv_iter;
  return std::make_tuple(v0, v1, v2);
}

template<typename MeshT>
inline auto face_points(MeshT& mesh, FaceHandle fh) ->
std::tuple<
  decltype(mesh.point(VertexHandle())),
  decltype(mesh.point(VertexHandle())),
  decltype(mesh.point(VertexHandle()))
>
{
  typedef decltype(mesh.point(VertexHandle())) PT;
  auto fv_iter = mesh.cfv_begin(fh);
  VertexHandle v0 = *fv_iter; fv_iter++;
  VertexHandle v1 = *fv_iter; fv_iter++;
  VertexHandle v2 = *fv_iter;
  return std::tuple<PT, PT, PT>(mesh.point(v0), mesh.point(v1), mesh.point(v2));
}

template<typename MeshT>
inline auto face_points(MeshT& mesh, FaceHandle fh, VertexHandle vh) ->
std::tuple<
  decltype(mesh.point(VertexHandle())),
  decltype(mesh.point(VertexHandle())),
  decltype(mesh.point(VertexHandle()))
>
{
  typedef decltype(mesh.point(VertexHandle())) PT;
  // returned points start with vh.
  HalfedgeHandle hh = mesh.halfedge_handle(fh);
  if (mesh.from_vertex_handle(hh) == vh)
    hh = mesh.prev_halfedge_handle(hh);
  else if (mesh.to_vertex_handle(hh) != vh)
    hh = mesh.next_halfedge_handle(hh);
  VertexHandle v0 = mesh.to_vertex_handle(hh); hh = mesh.next_halfedge_handle(hh);
  VertexHandle v1 = mesh.to_vertex_handle(hh); hh = mesh.next_halfedge_handle(hh);
  VertexHandle v2 = mesh.to_vertex_handle(hh);
  return std::tuple<PT, PT, PT>(mesh.point(v0), mesh.point(v1), mesh.point(v2));
}

template<typename MeshT>
inline std::tuple<VertexHandle, VertexHandle> from_to_vh(MeshT& mesh, HalfedgeHandle hh)
{
  return std::make_tuple(mesh.from_vertex_handle(hh), mesh.to_vertex_handle(hh));
}

template<typename MeshT>
inline auto from_to_p(MeshT& mesh, HalfedgeHandle hh)
-> std::tuple<decltype(mesh.point(VertexHandle())), decltype(mesh.point(VertexHandle()))>
{
  return std::make_tuple(
    mesh.point(mesh.from_vertex_handle(hh)),
    mesh.point(mesh.to_vertex_handle(hh)));
}

template<typename MeshT>
inline HalfedgeHandle opposite_halfedge_handle(MeshT& mesh, FaceHandle fh, VertexHandle vh)
{
  HalfedgeHandle hh = mesh.halfedge_handle(fh);
  if (mesh.to_vertex_handle(hh) == vh)
    return mesh.prev_halfedge_handle(hh);
  else if (mesh.from_vertex_handle(hh) == vh)
    return mesh.next_halfedge_handle(hh);
  else
    return hh;
}

}// namespace ImageTriSimp
}// namespace GCLF