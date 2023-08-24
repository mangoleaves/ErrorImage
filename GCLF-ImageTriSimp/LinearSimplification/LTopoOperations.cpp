#pragma once
#include "LTopoOperations.h"

namespace GCLF
{
namespace ImageTriSimp
{

std::vector<FaceHandle> faces_adjacent_to_edge(const LMeshT* mesh, EdgeHandle e)
{
  std::vector<FaceHandle> one_ring_faces;
  HalfedgeHandle he = mesh->halfedge_handle(e, 0);
  if (mesh->is_boundary(e))
  {
    if (mesh->is_boundary(he))
      he = mesh->opposite_halfedge_handle(he);
    one_ring_faces = { mesh->face_handle(he) };
  }
  else
    one_ring_faces = { mesh->face_handle(he), mesh->opposite_face_handle(he) };
  return one_ring_faces;
}

void one_ring_faces_around_vertex(
  const LMeshT* mesh,
  VertexHandle vh,
  std::set<FaceHandle>& faces)
{
  for (FaceHandle fh : mesh->vf_range(vh))
  {
    faces.insert(fh);
  }
}

void one_ring_faces_around_edge(
  const LMeshT* mesh,
  HalfedgeHandle hh,
  std::set<FaceHandle>& faces)
{
  faces.clear();
  VertexHandle v_from = mesh->from_vertex_handle(hh);
  VertexHandle v_to = mesh->to_vertex_handle(hh);
  one_ring_faces_around_vertex(mesh, v_from, faces);
  one_ring_faces_around_vertex(mesh, v_to, faces);
}

void vertices_around_faces(
  const LMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<VertexHandle>& vertices)
{
  for (FaceHandle fh : faces)
  {
    HalfedgeHandle hh = mesh->halfedge_handle(fh);
    vertices.insert(mesh->from_vertex_handle(hh));
    vertices.insert(mesh->to_vertex_handle(hh));
    vertices.insert(mesh->opposite_vh(hh));
  }
}

void edges_around_faces(
  const LMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<EdgeHandle>& edges)
{
  for (FaceHandle fh : faces)
  {
    HalfedgeHandle hh = mesh->halfedge_handle(fh);
    edges.insert(mesh->edge_handle(hh));
    edges.insert(mesh->edge_handle(mesh->next_halfedge_handle(hh)));
    edges.insert(mesh->edge_handle(mesh->prev_halfedge_handle(hh)));
  }
}

void extend_faces_by_one_ring(
  const LMeshT* mesh,
  std::set<FaceHandle>& faces)
{
  std::set<VertexHandle> vertices;
  vertices_around_faces(mesh, faces, vertices);
  for (VertexHandle vh : vertices)
  {
    one_ring_faces_around_vertex(mesh, vh, faces);
  }
}

/// @brief given faces, collect n ring extended-faces.
/// @param [in] mesh
/// @param [in] faces input faces
/// @param [in] stencil_ring_size n ring, can be 1,2,3...
/// @param [out] extended_faces n ring faces, including input faces.
void extend_faces(
  const LMeshT* mesh,
  const std::set<FaceHandle>& faces,
  int stencil_ring_size,
  std::set<FaceHandle>& extended_faces)
{
  extended_faces.clear();
  extended_faces = faces;
  for (int i = 0; i < stencil_ring_size; ++i)
  {
    extend_faces_by_one_ring(mesh, extended_faces);
  }
}

/// @brief find two ring vertices of given vertex handle.
/// then find one ring edges of two ring vertices.
std::vector<EdgeHandle> find_1rv_1re(const LMeshT* mesh, VertexHandle v)
{
  std::set<EdgeHandle> surrounding_edges;
  for (VertexHandle vv : mesh->vv_range(v))
  {
    for (EdgeHandle ve : mesh->ve_range(vv))
      surrounding_edges.insert(ve);
  }
  return std::vector<EdgeHandle>(surrounding_edges.begin(), surrounding_edges.end());
}

/*
void init_one_ring_faces(LMeshT* mesh, FaceHandle fh)
{
  std::set<FaceHandle>& one_ring_faces = mesh->data(fh).one_ring_faces;
  one_ring_faces.clear();
  for (VertexHandle fvh : mesh->fv_range(fh))
    for (FaceHandle vfh : mesh->vf_range(fvh))
      one_ring_faces.insert(vfh);
}

void init_one_ring_faces(LMeshT* mesh)
{
  for (FaceHandle fh : mesh->faces())
    init_one_ring_faces(mesh, fh);
}
*/
}// namespace ImageTriSimp
}// namespace GCLF