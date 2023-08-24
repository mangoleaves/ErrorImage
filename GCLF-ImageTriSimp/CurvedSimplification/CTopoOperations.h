#pragma once
#include "Mesh/BezierMeshDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{
using OpenMesh::TriConnectivity;
/********************************/
/*      neighbor search         */
/********************************/

std::vector<FaceHandle> faces_adjacent_to_edge(const BMeshT* mesh, EdgeHandle e);

void one_ring_faces_around_vertex(
  const BMeshT* mesh,
  VertexHandle vh,
  std::set<FaceHandle>& faces);

void one_ring_faces_around_edge(
  const BMeshT* mesh,
  HalfedgeHandle hh,
  std::set<FaceHandle>& faces);

void vertices_around_faces(
  const BMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<VertexHandle>& vertices);

void edges_around_faces(
  const BMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<EdgeHandle>& edges);

void extend_faces_by_one_ring(
  const BMeshT* mesh,
  std::set<FaceHandle>& faces);

void extend_faces(
  const BMeshT* mesh,
  const std::set<FaceHandle>& faces,
  int stencil_ring_size,
  std::set<FaceHandle>& extended_faces);

std::vector<EdgeHandle> find_1rv_1re(const BMeshT* mesh, VertexHandle v);
}// namespace ImageTriSimp
}// namespace GCLF