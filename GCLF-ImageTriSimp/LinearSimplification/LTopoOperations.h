#pragma once
#include "Mesh/LinearMeshDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{
/********************************/
/*      neighbor search         */
/********************************/

std::vector<FaceHandle> faces_adjacent_to_edge(const LMeshT* mesh, EdgeHandle e);

void one_ring_faces_around_vertex(
  const LMeshT* mesh,
  VertexHandle vh,
  std::set<FaceHandle>& faces);

void one_ring_faces_around_edge(
  const LMeshT* mesh,
  HalfedgeHandle hh,
  std::set<FaceHandle>& faces);

void vertices_around_faces(
  const LMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<VertexHandle>& vertices);

void edges_around_faces(
  const LMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<EdgeHandle>& edges);

void extend_faces_by_one_ring(
  const LMeshT* mesh,
  std::set<FaceHandle>& faces);

void extend_faces(
  const LMeshT* mesh,
  const std::set<FaceHandle>& faces,
  int stencil_ring_size,
  std::set<FaceHandle>& extended_faces);

std::vector<EdgeHandle> find_1rv_1re(const LMeshT* mesh, VertexHandle v);

// void init_one_ring_faces(LMeshT* mesh, FaceHandle fh);
// void init_one_ring_faces(LMeshT* mesh);

}// namespace ImageTriSimp
}// namespace GCLF