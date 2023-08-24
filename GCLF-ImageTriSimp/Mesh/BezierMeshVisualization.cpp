#include "BezierMeshVisualization.h"

namespace GCLF
{
namespace ImageTriSimp
{

void linear_approximation(const BMeshT& bmesh, uint32_t density, LMeshT& lmesh)
{
  if (bmesh.n_faces() == 0)
    return;

  lmesh.clear();

  // If a vertex of bmesh is added, we store the index in lmesh.
  std::vector<VertexHandle> added_vertices(bmesh.n_vertices(), lmesh.InvalidVertexHandle);
  // If an edge of bmesh is added, we store the indices in lmesh.
  // (An edge of bmesh contains several points.)
  std::vector<std::vector<VertexHandle>> added_edges(bmesh.n_halfedges());

  // The topology of linear approximation of each Bezier triangle is same.
  // We first calculate the triangles.
  std::vector<Index3> triangles;
  Points points;
  {
    const BezierFace& first_bface = bmesh.data(FaceHandle(0)).bface;
    BezierTriangle::linear_approximation(
      first_bface.ctrl_pnt, first_bface.degree,
      density, true, points, triangles);
  }

  // We construct linear approximation triangle by triangle.
  // But we need to connect these triangles to get a connected surface.

  std::vector<VertexHandle> map_to_global_indices; // per face, re-use it.
  for (FaceHandle b_fh : bmesh.faces())
  {
    // if (!bmesh.is_valid_handle(b_fh))
    //  continue;
    
    const BezierFace& bface = bmesh.data(b_fh).bface;
    BezierTriangle::linear_approximation(
      bface.ctrl_pnt, bface.degree,
      density, false, points, triangles);
    map_to_global_indices.clear();
    map_to_global_indices.resize(points.size(), lmesh.InvalidVertexHandle);

    // Add all local points to lmesh and store corresponding indices into map_to_global_indices.
    // Note that we should avoid duplicate vertex, thus we check all vertices and halfedges.
    for(VertexHandle vh : bmesh.fv_range(b_fh))
    {
      Index vert_idx = bface.find_local_vert_idx(vh.idx());
      Index local_idx = BezierTriangle::index_of_vert(density, vert_idx);
      if (!added_vertices[vh.idx()].is_valid())
      {
        VertexHandle added_vh = lmesh.add_vertex(points[local_idx]);
        added_vertices[vh.idx()] = added_vh;
        map_to_global_indices[local_idx] = added_vh;
      }
      else
      {
        map_to_global_indices[local_idx] = added_vertices[vh.idx()];
      }
    }

    for (HalfedgeHandle heh : bmesh.fh_range(b_fh))
    {
      VertexHandle v_from = bmesh.from_vertex_handle(heh);
      VertexHandle v_to = bmesh.to_vertex_handle(heh);
      // We check if the opposite halfedge is added,
      // assuming that the surface is manifold.
      HalfedgeHandle opp_heh = bmesh.opposite_halfedge_handle(heh);
      // We find the local edge index of the halfedge in Bezier triangle.
      Index bedge_idx = bface.find_local_edge_idx(v_from.idx(), v_to.idx());
      // We find the local vertex indices of the local edge in Bezier triangle.
      Indices pnts_idx_on_bedge = BezierTriangle::indices_on_edge(density, bedge_idx);
      if (added_edges[opp_heh.idx()].empty())
      {
        // Add the points in the edge and store the global indices.
        auto& added_edge = added_edges[heh.idx()];    added_edge.reserve(pnts_idx_on_bedge.size() - 2);
        for (size_t i = pnts_idx_on_bedge.size() - 2;i >= 1;i--)  // reverse order!
        {
          Index local_idx = pnts_idx_on_bedge[i];
          VertexHandle added_vh = lmesh.add_vertex(points[local_idx]);
          added_edge.push_back(added_vh);
          map_to_global_indices[local_idx] = added_vh;
        }
      }
      else
      {
        // find the indices of the points in the edge.
        auto& opp_added_edge = added_edges[opp_heh.idx()];
        for (size_t i = 1;i < pnts_idx_on_bedge.size() - 1;i++)
        {
          map_to_global_indices[pnts_idx_on_bedge[i]] = opp_added_edge[i - 1];
        }
      }
    }
    // Add all interior points of linear approximation of Bezier triangle to lmesh.
    Indices indices_interior = BezierTriangle::indices_interior(density);
    for(Index local_idx : indices_interior)
    {
      map_to_global_indices[local_idx] = lmesh.add_vertex(points[local_idx]);
    }

    // Finally, we add triangles into lmesh.
    for(const Index3& tri : triangles)
    {
      VertexHandle v0 = map_to_global_indices[tri[0]];
      VertexHandle v1 = map_to_global_indices[tri[1]];
      VertexHandle v2 = map_to_global_indices[tri[2]];
      lmesh.add_face(v0, v1, v2);
    }
  }
}

/// @brief Get the control mesh of a Bezier mesh.
/// @param [in] bmesh Bezier mesh.
/// @param [out] cmesh Control mesh, containing triangles.
/// @param [out] jmesh Junction mesh, containing edges that are junction of two triangles.
void get_ctrl_mesh(const BMeshT& bmesh, LMeshT& cmesh, LMeshT& jmesh)
{
  cmesh.clear();
  jmesh.clear();
  // Add all control points to cmesh.
  for (const Vec3d& cp : bmesh.ctrl_pnt)
  {
    cmesh.add_vertex(cp);
    jmesh.add_vertex(cp);
  }

  std::vector<Index3> triangles = BezierTriangle::get_ctrl_triangles(bmesh.bezier_degree);
  std::vector<Index2> junctions = BezierTriangle::get_junctions_of_ctrl_triangles(bmesh.bezier_degree);

  // Add all control trianlges to cmesh.
  for (FaceHandle fh : bmesh.faces())
  {
    const BezierFace& bface = bmesh.data(fh).bface;

    for(const Index3& tri : triangles)
    {
      VertexHandle vh0(bface.ctrl_pnt_idx[tri[0]]);
      VertexHandle vh1(bface.ctrl_pnt_idx[tri[1]]);
      VertexHandle vh2(bface.ctrl_pnt_idx[tri[2]]);

      cmesh.add_face(vh0, vh1, vh2);
    }

    for(const Index2& junc : junctions)
    {
      VertexHandle vh0(bface.ctrl_pnt_idx[junc[0]]);
      VertexHandle vh1(bface.ctrl_pnt_idx[junc[1]]);
      VertexHandle mid_vh = jmesh.add_vertex((bface.ctrl_pnt[junc[0]] + bface.ctrl_pnt[junc[1]]) * 0.5);

      jmesh.add_face(vh0, vh1, mid_vh);
    }
  }
}

}// namespace ImageTriSimp
}// namespace GCLF