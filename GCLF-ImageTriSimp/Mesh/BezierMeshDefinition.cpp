#include "BezierMeshDefinition.h"
#include "LinearMeshDefinition.h"
#include "CurvedSimplification/CTopoOperations.h"
#include <fstream>

namespace GCLF
{
namespace ImageTriSimp
{

BMeshT::FaceData::FaceData(const FaceData& rhs)
{
  bface = rhs.bface;
  color = rhs.color;
  pixels = rhs.pixels;
  error = rhs.error;
}

BMeshT::FaceData::FaceData(FaceData&& rhs) noexcept
{
  bface = std::move(rhs.bface);
  color = std::move(rhs.color);
  pixels = std::move(rhs.pixels);
  error = rhs.error;
}

BMeshT::FaceData& BMeshT::FaceData::operator=(const FaceData& rhs)
{
  bface = rhs.bface;
  color = rhs.color;
  pixels = rhs.pixels;
  error = rhs.error;
  return *this;
}

BMeshT::FaceData& BMeshT::FaceData::operator=(FaceData&& rhs) noexcept
{
  bface = std::move(rhs.bface);
  color = std::move(rhs.color);
  pixels = std::move(rhs.pixels);
  error = rhs.error;
  return *this;
}

BMeshT::BMeshT()
{
  request_vertex_status();
  request_edge_status();
  request_face_status();
  request_halfedge_status();
  this->add_property(data_vpph_);
  this->add_property(data_hpph_);
  this->add_property(data_epph_);
  this->add_property(data_fpph_);
}

BMeshT::BMeshT(const BMeshT& rhs)
  :TriConnectivity(*static_cast<const OM::TriConnectivity*>(&rhs))
{
  copy_from(rhs);
}

BMeshT::BMeshT(BMeshT&& rhs)noexcept
  :TriConnectivity(std::move(*static_cast<const OM::TriConnectivity*>(&rhs)))
{
  move_from(std::move(rhs));
}

BMeshT& BMeshT::operator=(const BMeshT& rhs)
{
  *static_cast<OM::TriConnectivity*>(this) = *static_cast<const OM::TriConnectivity*>(&rhs);
  copy_from(rhs);
  return *this;
}

BMeshT& BMeshT::operator=(BMeshT&& rhs)noexcept
{
  *static_cast<OM::TriConnectivity*>(this) = std::move(*static_cast<OM::TriConnectivity*>(&rhs));
  move_from(std::move(rhs));
  return *this;
}

Index BMeshT::add_ctrl_pnt(const Vec3d& p)
{
  ctrl_pnt.push_back(p);
  ctrl_pnt_deleted.push_back(false);
  return (Index)(ctrl_pnt.size() - 1);
}

/// @brief Mark a control point as deleted.
/// @param ctrl_pnt_idx We assume that this control point is isolated.
void BMeshT::delete_ctrl_pnt(Index ctrl_pnt_idx)
{
  ctrl_pnt_deleted[ctrl_pnt_idx] = true;
}

/// @brief Initialize the second order Bezier mesh from a piecewise linear mesh.
/// @param lmesh A piecewise linear mesh.
void BMeshT::initialize(LMeshT* lmesh)
{
  // make sure that two meshes are clean.
  clear();
  lmesh->garbage_collection();

  // 1. We initialize a bezier mesh whose degree is 1.
  bezier_degree = 1;
  // (1.1) Add vertices of linear mesh into bezier mesh.
  //       Now indices of vertices of linear mesh, contrl points and
  //       vertices of bezier mesh are same.
  for (VertexHandle vh : lmesh->vertices())
  {
    Index ctrl_idx = add_ctrl_pnt(lmesh->point(vh));
    VertexHandle new_vh = TriConnectivity::add_vertex();
  }
  // (1.2) Add faces of linear mesh into bezier mesh.
  for (FaceHandle fh : lmesh->faces())
  {
    auto [v0, v1, v2] = face_vertices(*lmesh, fh);
    FaceHandle fh = add_face(v0, v1, v2);
    BezierFace& bface = data(fh).bface;
    bface.degree = 1;
    bface.vert_idx = Index3(v0.idx(), v1.idx(), v2.idx());
    bface.ctrl_pnt_idx = { v0.idx(), v1.idx(), v2.idx() };
    bface.ctrl_pnt = { ctrl_pnt[v0.idx()], ctrl_pnt[v1.idx()], ctrl_pnt[v2.idx()] };
  }

  // 2. We elevate the degree of the Bezier mesh to the expected degree.
  degree_elevation();

  store_extra_cpi_to_ve();
}

/// @brief Elevate the degree of mesh by one.
void BMeshT::degree_elevation()
{
  clear_ctrl_pnt();

  Indices added_vertices(n_vertices(), -1);
  std::vector<Indices> added_edges(n_halfedges());

  // We elevate degree of Bezier triangles respectively.
  for (FaceHandle fh : faces())
  {
    BezierFace& bface = data(fh).bface;
    // Elevate the degree and calculate new control points.
    bface.ctrl_pnt = BezierTriangle::degree_elevation(bface.ctrl_pnt, bface.degree);
    bface.degree++;
    // Find new indices for control points.
    bface.ctrl_pnt_idx.clear();
    bface.ctrl_pnt_idx.resize(bface.ctrl_pnt.size());

    // (1) Find new indices for corner control points.
    for (VertexHandle vh : fv_range(fh))
    {
      if (added_vertices[vh.idx()] == -1)
        added_vertices[vh.idx()] = add_ctrl_pnt(bface.global_corner_ctrlpnt(vh.idx()));
      bface.global_corner_ctrlpnt_idx(vh.idx()) = added_vertices[vh.idx()];
    }
    // (2) Find new indices for control points on edges.
    for (HalfedgeHandle heh : fh_range(fh))
    {
      VertexHandle v_from = from_vertex_handle(heh);
      VertexHandle v_to = to_vertex_handle(heh);
      // We check if the opposite halfedge is added,
      // assuming that the surface is manifold.
      HalfedgeHandle opp_heh = opposite_halfedge_handle(heh);

      if (added_edges[opp_heh.idx()].empty())
      {
        // Add the points in the edge and store the global indices.
        Indices& added_edge = added_edges[heh.idx()];    added_edge.reserve(bface.degree - 1);
        for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_mr(v_from.idx(), v_to.idx()))
        {
          cp_idx = add_ctrl_pnt(cp);
          added_edge.push_back(cp_idx);
        }
      }
      else
      {
        // find the indices of the points in the edge.
        Indices& opp_added_edge = added_edges[opp_heh.idx()];
        size_t i = 0;
        for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()))
        {
          cp_idx = opp_added_edge[i++];
          cp = ctrl_pnt[cp_idx];  // avoid numeric error.
        }
      }
    }

    // (3) Add all interior control points and store the indices.
    for (auto [cp_idx, cp] : bface.global_interior_ctrlpnt_ms())
    {
      cp_idx = add_ctrl_pnt(cp);
    }
  }
  bezier_degree++;
}

/// @brief When we need to build (or rebuild) relation between control points of each faces and
/// the global data ctrl_pnt, we call this function.
/// @note we also make the common control points between faces be same to avoid numeric errors.
void BMeshT::build_ctrlpnt_idx()
{
  // We assume that control points have been added into all faces
  // and the mesh is clean (garbage_collection operated).
  clear_ctrl_pnt();

  Indices added_vertices(n_vertices(), -1);
  std::vector<Indices> added_edges(n_halfedges());

  for (FaceHandle fh : faces())
  {
    BezierFace& bface = data(fh).bface;
    // Find new indices for control points.
    bface.ctrl_pnt_idx.clear();
    bface.ctrl_pnt_idx.resize(bface.ctrl_pnt.size());

    // (1) Find new indices for corner control points.
    for (VertexHandle vh : fv_range(fh))
    {
      if (added_vertices[vh.idx()] == -1)
        added_vertices[vh.idx()] = add_ctrl_pnt(bface.global_corner_ctrlpnt(vh.idx()));
      bface.global_corner_ctrlpnt_idx(vh.idx()) = added_vertices[vh.idx()];
    }
    // (2) Find new indices for control points on edges.
    for (HalfedgeHandle heh : fh_range(fh))
    {
      VertexHandle v_from = from_vertex_handle(heh);
      VertexHandle v_to = to_vertex_handle(heh);
      // We check if the opposite halfedge is added,
      // assuming that the surface is manifold.
      HalfedgeHandle opp_heh = opposite_halfedge_handle(heh);

      if (added_edges[opp_heh.idx()].empty())
      {
        // Add the points in the edge and store the global indices.
        Indices& added_edge = added_edges[heh.idx()];    added_edge.reserve(bface.degree - 1);
        for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_mr(v_from.idx(), v_to.idx()))
        {
          cp_idx = add_ctrl_pnt(cp);
          added_edge.push_back(cp_idx);
        }
      }
      else
      {
        // find the indices of the points in the edge.
        Indices& opp_added_edge = added_edges[opp_heh.idx()];
        size_t i = 0;
        for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()))
        {
          cp_idx = opp_added_edge[i++];
          cp = ctrl_pnt[cp_idx];  // avoid numeric error.
        }
      }
    }

    // (3) Add all interior control points and store the indices.
    for (auto [cp_idx, cp] : bface.global_interior_ctrlpnt_ms())
    {
      cp_idx = add_ctrl_pnt(cp);
    }
  }

  store_extra_cpi_to_ve();
}

/// @brief After updating global control points array,
/// we need to update local control points in each Bezier face.
void BMeshT::sync_local_ctrlpnt()
{
  for (FaceHandle fh : faces())
  {
    BezierFace& bface = data(fh).bface;
    for (auto [cp_idx, cp] : bface.global_ctrlpnt_ms())
      cp = ctrl_pnt[cp_idx];
  }
}

void BMeshT::store_extra_cpi_to_ve()
{
  ASSERT(bezier_degree == 2, "wrong bezier degree.");

  for (VertexHandle vh : vertices())
  {
    data(vh).cpi = data(*vf_begin(vh)).bface.global_corner_ctrlpnt_idx(vh.idx());
  }
  for (EdgeHandle eh : edges())
  {
    HalfedgeHandle hh = halfedge_handle(eh, 0);
    if (is_boundary(hh))
      hh = opposite_halfedge_handle(hh);
    FaceHandle fh = face_handle(hh);
    VertexHandle v_from = from_vertex_handle(hh);
    VertexHandle v_to = to_vertex_handle(hh);
    for (auto [cpi, cp] : data(fh).bface.global_edge_interior_ctrlpnt_cs(v_from.idx(), v_to.idx()))
      data(eh).cpi = cpi;
  }
}

/// @brief Split a face of the Bezier mesh.
/// @param uvw Barycentric coordinates.
/// @return The split vertex.
/// @note Ask Guo Jiapeng for a figure that illustrates the handles and their relations.
/// @bug Update cpi in VertexData and EdgeData.
VertexHandle BMeshT::split_face(FaceHandle fh_to_split, Vec3d uvw)
{
  ASSERT(TriConnectivity::is_valid_handle(fh_to_split), "Invalid face.");

  // Split the triangle to three sub-triangles geometrically.
  BezierFace old_bface = data(fh_to_split).bface;
  std::array<Points, 3> sub_ctrl_pnts;
  Vec3d sub_center_ctrl_pnt;
  BezierTriangle::subdivision_1_3(old_bface.ctrl_pnt, bezier_degree, uvw, sub_ctrl_pnts, sub_center_ctrl_pnt);

  // Find the 0-th halfedge.
  HalfedgeHandle hh_0th = halfedge_handle(fh_to_split);
  {
    Index local_edge_idx = old_bface.find_local_edge_idx(
      from_vertex_handle(hh_0th).idx(),
      to_vertex_handle(hh_0th).idx());

    // Rotate to the 0-th edge.
    if (local_edge_idx == 1)
      hh_0th = next_halfedge_handle(next_halfedge_handle(hh_0th));
    if (local_edge_idx == 2)
      hh_0th = next_halfedge_handle(hh_0th);
  }

  // Split the face topologically.
  VertexHandle center_vh = new_vertex();
  split(fh_to_split, center_vh);

  // Put the new control points to three sub-triangles 
  HalfedgeHandle hh = hh_0th;
  for (Index local_edge_idx = 0;local_edge_idx < 3;local_edge_idx++)
  {
    FaceHandle fh_of_hh = face_handle(hh);
    BezierFace& bface = data(fh_of_hh).bface;
    bface.vert_idx = old_bface.vert_idx;
    bface.ctrl_pnt = std::move(sub_ctrl_pnts[local_edge_idx]);
    bface.ctrl_pnt_idx = old_bface.ctrl_pnt_idx;
    // Note: only indices on old edges are right, others are wrong.
    // We will fix them later.
    hh = next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(hh)));
  }

  // Delete old interior control points.
  for (auto [cp_idx, cp] : old_bface.global_interior_ctrlpnt_ms())
    delete_ctrl_pnt(cp_idx);

  // Add the center ctrl pnt into global array.
  Index global_center_ctrl_pnt_idx = add_ctrl_pnt(sub_center_ctrl_pnt);

  // Fix the indices and add new control points to global array.
  for (HalfedgeHandle hh : vih_range(center_vh))
  {
    VertexHandle v_from = from_vertex_handle(hh);
    BezierFace& bface = data(face_handle(hh)).bface;
    bface.degree = bezier_degree;
    // Fix the index on the subdivide vh.
    {
      Index local_subdivide_vidx = (bface.find_local_vert_idx(v_from.idx()) + 1) % 3;
      bface.vert_idx[local_subdivide_vidx] = center_vh.idx();
      Index local_ctrl_pnt_idx = BezierTriangle::index_of_vert(bezier_degree, local_subdivide_vidx);
      bface.ctrl_pnt_idx[local_ctrl_pnt_idx] = global_center_ctrl_pnt_idx;
    }
    // Fix the indices in the triangle.
    {
      for (auto [cp_idx, cp] : bface.global_interior_ctrlpnt_ms())
        cp_idx = add_ctrl_pnt(cp);
    }
  }

  for (HalfedgeHandle hh : vih_range(center_vh))
  {
    HalfedgeHandle opp_hh = opposite_halfedge_handle(hh);
    VertexHandle v_from = from_vertex_handle(hh);
    VertexHandle v_to = to_vertex_handle(hh);   // v_to is also center_vh
    BezierFace& bface = data(face_handle(hh)).bface;
    BezierFace& opp_bface = data(face_handle(opp_hh)).bface;

    // Fix the indices on the common edge between two triangles.
    for (auto& [cpi, opp_cpi] : BF_Zip(bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()),
      opp_bface.global_edge_interior_ctrlpnt_mr(v_to.idx(), v_from.idx())))
    {
      auto [idx, cp] = cpi;
      auto [opp_idx, opp_cp] = opp_cpi;
      opp_idx = idx = add_ctrl_pnt(cp);
    }
  }
  return center_vh;
}

/// @brief Split a Bezier mesh at an edge.
/// @param hh_to_split A halfedge. from_vh -> u=1; to_vh -> v=1.
/// @param uv Barycentric coordinates.
/// @return The split vertex.
/// @note Ask Guo Jiapeng for a figure that illustrates the handles and their relations.
/// @bug Update cpi in VertexData and EdgeData.
VertexHandle BMeshT::split_edge(HalfedgeHandle hh_to_split, Vec2d uv)
{
  ASSERT(TriConnectivity::is_valid_handle(hh_to_split), "Invalid edge.");

  // Adjust the handles.
  EdgeHandle eh = edge_handle(hh_to_split);
  if (halfedge_handle(eh, 1) == hh_to_split)
  {
    hh_to_split = opposite_halfedge_handle(hh_to_split);
    uv = { uv[1], uv[0] };
  }

  // Pre-declarations.
  HalfedgeHandle hh = hh_to_split;
  HalfedgeHandle opp_hh = opposite_halfedge_handle(hh_to_split);
  FaceHandle fh = face_handle(hh);
  FaceHandle opp_fh = face_handle(opp_hh);
  VertexHandle v_from = from_vertex_handle(hh);
  VertexHandle v_to = to_vertex_handle(hh);

  BezierFace old_bface, old_opp_bface;
  if (!is_boundary(hh)) old_bface = data(fh).bface;
  if (!is_boundary(opp_hh)) old_opp_bface = data(opp_fh).bface;

  Vec2d vu = { uv[1], uv[0] };
  Vec3d sub_center_ctrl_pnt;

  // Split the edge topologically.
  VertexHandle center_vh = new_vertex();
  split(eh, center_vh);

  // Split the edge geometrically.
  if (!is_boundary(hh))
  {
    // The variable names come from OpenMesh.
    // Two sub-edges of the original edge: h0 and t1.
    HalfedgeHandle h0 = hh;
    HalfedgeHandle t1 = prev_halfedge_handle(opposite_halfedge_handle(prev_halfedge_handle(h0)));
    FaceHandle f0 = face_handle(h0);  // the face contains v_to.
    FaceHandle f1 = face_handle(t1);  // the face contains v_from.

    // Subdivide the face at a point (uv) on the edge (h0).
    Index local_edge_idx = old_bface.find_local_edge_idx(v_from.idx(), v_to.idx());
    std::array<Points, 2> sub_ctrl_pnts;
    BezierTriangle::subdivision_1_2(
      old_bface.ctrl_pnt, bezier_degree, local_edge_idx, uv,
      sub_ctrl_pnts, sub_center_ctrl_pnt);

    // Put the new control points into two sub-triangles.
    BezierFace& bface0 = data(f0).bface;
    BezierFace& bface1 = data(f1).bface;
    bface0.vert_idx = old_bface.vert_idx;
    bface0.ctrl_pnt = std::move(sub_ctrl_pnts[1]);
    bface0.ctrl_pnt_idx = old_bface.ctrl_pnt_idx;
    bface1.vert_idx = old_bface.vert_idx;
    bface1.ctrl_pnt = std::move(sub_ctrl_pnts[0]);
    bface1.ctrl_pnt_idx = old_bface.ctrl_pnt_idx;
    // Note: only indices on old edges are right, others are wrong.
    // We will fix them later.

    // Delete old control points.
    for (auto [cp_idx, cp] : old_bface.global_interior_ctrlpnt_cs())
      delete_ctrl_pnt(cp_idx);
    for (auto [cp_idx, cp] : old_bface.global_edge_interior_ctrlpnt_cs(local_edge_idx))
      delete_ctrl_pnt(cp_idx);
  }
  if (!is_boundary(opp_hh))
  {
    // The variable names come from OpenMesh.
    // Two sub-edges of the original edge.
    HalfedgeHandle o0 = opp_hh;
    HalfedgeHandle e1 = next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(o0)));
    FaceHandle f3 = face_handle(o0);  // the face contains v_to.
    FaceHandle f2 = face_handle(e1);  // the face contains v_from.

    // Subdivide the face at a point (vu) on the edge (o0).
    Index local_edge_idx = old_opp_bface.find_local_edge_idx(v_to.idx(), v_from.idx());
    std::array<Points, 2> sub_ctrl_pnts;
    BezierTriangle::subdivision_1_2(
      old_opp_bface.ctrl_pnt, bezier_degree, local_edge_idx, vu,
      sub_ctrl_pnts, sub_center_ctrl_pnt);

    // Put the new control points into two sub-triangles.
    BezierFace& bface3 = data(f3).bface;
    BezierFace& bface2 = data(f2).bface;
    bface3.vert_idx = old_opp_bface.vert_idx;
    bface3.ctrl_pnt = std::move(sub_ctrl_pnts[0]);
    bface3.ctrl_pnt_idx = old_opp_bface.ctrl_pnt_idx;
    bface2.vert_idx = old_opp_bface.vert_idx;
    bface2.ctrl_pnt = std::move(sub_ctrl_pnts[1]);
    bface2.ctrl_pnt_idx = old_opp_bface.ctrl_pnt_idx;
    // Note: only indices on old edges are right, others are wrong.
    // We will fix them later.

    // Delete old control points.
    for (auto [cp_idx, cp] : old_opp_bface.global_interior_ctrlpnt_cs())
      delete_ctrl_pnt(cp_idx);
    for (auto [cp_idx, cp] : old_opp_bface.global_edge_interior_ctrlpnt_cs(local_edge_idx))
      delete_ctrl_pnt(cp_idx);
  }

  // Add the center ctrl pnt into global array.
  Index global_center_ctrl_pnt_idx = add_ctrl_pnt(sub_center_ctrl_pnt);

  // Fix the indices and add new control points to global array.
  for (HalfedgeHandle hh : vih_range(center_vh))
  {
    if (is_boundary(hh))
      continue;
    VertexHandle v_from = from_vertex_handle(hh);
    BezierFace& bface = data(face_handle(hh)).bface;
    bface.degree = bezier_degree;
    // Fix the index on the subdivide center (center_vh).
    {
      Index local_subdivide_vidx = (bface.find_local_vert_idx(v_from.idx()) + 1) % 3;
      bface.vert_idx[local_subdivide_vidx] = center_vh.idx();
      Index local_ctrl_pnt_idx = BezierTriangle::index_of_vert(bezier_degree, local_subdivide_vidx);
      bface.ctrl_pnt_idx[local_ctrl_pnt_idx] = global_center_ctrl_pnt_idx;
    }
    // Fix the indices inside the triangle.
    {
      for (auto [cp_idx, cp] : bface.global_interior_ctrlpnt_ms())
        cp_idx = add_ctrl_pnt(cp);
    }
  }
  for (HalfedgeHandle hh : vih_range(center_vh))
  {
    if (!is_boundary(edge_handle(hh)))
    {
      // We simultaneously process two halfedges.
      HalfedgeHandle opp_hh = opposite_halfedge_handle(hh);
      VertexHandle v_from = from_vertex_handle(hh);
      VertexHandle v_to = to_vertex_handle(hh);   // v_to is also center_vh
      BezierFace& bface = data(face_handle(hh)).bface;
      BezierFace& opp_bface = data(face_handle(opp_hh)).bface;
      // Fix the indices on the common edge between two triangles.
      for (auto& [cpi, opp_cpi] : BF_Zip(bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()),
        opp_bface.global_edge_interior_ctrlpnt_mr(v_to.idx(), v_from.idx())))
      {
        auto [cp_idx, cp] = cpi;
        auto [opp_cp_idx, opp_cp] = opp_cpi;
        opp_cp_idx = cp_idx = add_ctrl_pnt(cp);
      }
    }
    else
    {
      // We only process one halfedge that is not boundary.
      HalfedgeHandle hh_not_boundary = is_boundary(hh) ? opposite_halfedge_handle(hh) : hh;
      VertexHandle v_from = from_vertex_handle(hh_not_boundary);
      VertexHandle v_to = to_vertex_handle(hh_not_boundary);   // v_to is also center_vh
      BezierFace& bface = data(face_handle(hh_not_boundary)).bface;
      // Fix the indices on the edge.
      for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()))
        cp_idx = add_ctrl_pnt(cp);
    }
  }
  return center_vh;
}

/// @brief Collapse an edge. The topology variation is same with the OpenMesh.
/// @param hh The edge to collapse.
/// @return Succeed?
/// @note Ask Guo Jiapeng for a figure that illustrates the handles and their relations.
bool BMeshT::collapse(HalfedgeHandle hh)
{
  if (!is_collapse_ok(hh))
    return false;

  // Pre-declarations.
  HalfedgeHandle phh = prev_halfedge_handle(hh);
  HalfedgeHandle nhh = next_halfedge_handle(hh);
  HalfedgeHandle ohh = opposite_halfedge_handle(hh);
  HalfedgeHandle pohh = prev_halfedge_handle(ohh);
  HalfedgeHandle nohh = next_halfedge_handle(ohh);

  VertexHandle v_from = from_vertex_handle(hh);
  VertexHandle v_to = to_vertex_handle(hh);
  Vec3d ctrl_pnt_of_v_to;
  Index ctrl_pnt_idx_of_v_from, ctrl_pnt_idx_of_v_to;
  // Collapse edge geometrically.

  // Delete control points on face_handle(hh)
  if (!is_boundary(hh))
  {
    FaceHandle fh = face_handle(hh);
    BezierFace& bface = data(fh).bface;
    Index local_edge_idx = bface.find_local_edge_idx(v_from.idx(), v_to.idx());

    // Delete control points on edge_handle(hh).
    for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_cs(local_edge_idx))
      delete_ctrl_pnt(cp_idx);

    ctrl_pnt_idx_of_v_from = bface.global_corner_ctrlpnt_idx(v_from.idx());
    ctrl_pnt_idx_of_v_to = bface.global_corner_ctrlpnt_idx(v_to.idx());
    ctrl_pnt_of_v_to = bface.global_corner_ctrlpnt(v_to.idx());

    // Delete control points on edge_handle(phh).
    for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_cs((local_edge_idx + 2) % 3))
      delete_ctrl_pnt(cp_idx);

    // Delete control points inside the face_handle(hh).
    for (auto [cp_idx, cp] : bface.global_interior_ctrlpnt_cs())
      delete_ctrl_pnt(cp_idx);
  }
  // Delete control points on face_handle(ohh)
  if (!is_boundary(ohh))
  {
    FaceHandle ofh = face_handle(ohh);
    BezierFace& bface = data(ofh).bface;
    Index local_edge_idx = bface.find_local_edge_idx(v_to.idx(), v_from.idx());

    if (is_boundary(hh))
    {
      for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_cs(local_edge_idx))
        delete_ctrl_pnt(cp_idx);

      ctrl_pnt_idx_of_v_from = bface.global_corner_ctrlpnt_idx(v_from.idx());
      ctrl_pnt_idx_of_v_to = bface.global_corner_ctrlpnt_idx(v_to.idx());
      ctrl_pnt_of_v_to = bface.global_corner_ctrlpnt(v_to.idx());
    }
    // else control points on edge_handle(ohh) have been deleted.

    // Delete control points on edge_handle(nohh).
    for (auto [cp_idx, cp] : bface.global_edge_interior_ctrlpnt_cs((local_edge_idx + 1) % 3))
      delete_ctrl_pnt(cp_idx);

    // Delete control points inside the face_handle(ohh).
    for (auto [cp_idx, cp] : bface.global_interior_ctrlpnt_cs())
      delete_ctrl_pnt(cp_idx);
  }
  // Update control points on opposite(phh).
  // The control points on opposite(phh) are changed to be same with control points on nhh.
  HalfedgeHandle ophh = opposite_halfedge_handle(phh);
  if (!is_boundary(hh) && !is_boundary(ophh))
  {
    BezierFace& bface_update = data(face_handle(ophh)).bface;
    BezierFace& bface_ref = data(face_handle(hh)).bface;

    for (auto& [update_cpi, ref_cpi] : BF_Zip(
      bface_update.global_edge_interior_ctrlpnt_ms(from_vertex_handle(ophh).idx(), to_vertex_handle(ophh).idx()),
      bface_ref.global_edge_interior_ctrlpnt_cs(from_vertex_handle(nhh).idx(), to_vertex_handle(nhh).idx())))
    {
      auto [update_cp_idx, update_cp] = update_cpi;
      auto [ref_cp_idx, ref_cp] = ref_cpi;
      update_cp_idx = ref_cp_idx;
      update_cp = ref_cp;
    }
  }
  // Update control points on opposite(nohh).
  // The control points on opposite(nohh) are changed to be same with control points on pohh.
  HalfedgeHandle onohh = opposite_halfedge_handle(nohh);
  if (!is_boundary(ohh) && !is_boundary(onohh))
  {
    BezierFace& bface_update = data(face_handle(onohh)).bface;
    BezierFace& bface_ref = data(face_handle(ohh)).bface;

    for (auto& [update_cpi, ref_cpi] : BF_Zip(
      bface_update.global_edge_interior_ctrlpnt_ms(from_vertex_handle(onohh).idx(), to_vertex_handle(onohh).idx()),
      bface_ref.global_edge_interior_ctrlpnt_cs(from_vertex_handle(pohh).idx(), to_vertex_handle(pohh).idx())))
    {
      auto [update_cp_idx, update_cp] = update_cpi;
      auto [ref_cp_idx, ref_cp] = ref_cpi;
      update_cp_idx = ref_cp_idx;
      update_cp = ref_cp;
    }
  }
  // Delete v_from and change v_from to v_to in one ring faces of v_from.
  delete_ctrl_pnt(ctrl_pnt_idx_of_v_from);
  for (HalfedgeHandle vih : vih_range(v_from))
  {
    if (is_boundary(vih))
      continue;

    BezierFace& bface = data(face_handle(vih)).bface;

    Index local_vert_idx = bface.find_local_vert_idx(v_from.idx());
    bface.vert_idx[local_vert_idx] = v_to.idx();

    Index cp_idx = BezierTriangle::index_of_vert(bezier_degree, local_vert_idx);
    bface.ctrl_pnt_idx[cp_idx] = ctrl_pnt_idx_of_v_to;
    bface.ctrl_pnt[cp_idx] = ctrl_pnt_of_v_to;
  }
  // Collapse the edge topologically.
  TriConnectivity::collapse(hh);

  return true;
}

/// @brief  Flip only for the second order Bezier mesh.
/// The second order Bezier triangle has no interior control point,
/// so it is easy to flip an edge.
bool BMeshT::flip_2(EdgeHandle eh)
{
  ASSERT(bezier_degree == 2, "This is not a second order Bezier mesh.");

  // Check if flip ok.
  if (!is_flip_ok(eh))
    return false;

  // Initialize relative handles.
  HalfedgeHandle a0, b0, a1, a2, b1, b2;
  VertexHandle   va0, va1, vb0, vb1;
  FaceHandle     fa, fb;

  a0 = halfedge_handle(eh, 0);
  b0 = halfedge_handle(eh, 1);

  a1 = next_halfedge_handle(a0);
  a2 = next_halfedge_handle(a1);

  b1 = next_halfedge_handle(b0);
  b2 = next_halfedge_handle(b1);

  va0 = to_vertex_handle(a0);
  va1 = to_vertex_handle(a1);

  vb0 = to_vertex_handle(b0);
  vb1 = to_vertex_handle(b1);

  fa = face_handle(a0);
  fb = face_handle(b0);

  // Do real flip
  TriConnectivity::flip(eh);
  // Update control points on two faces.
  // fa: va1, a2, a0, vb0, b1, vb1
  data(fa).bface.ctrl_pnt = {
    point(va1),
    point(edge_handle(a2)), point(edge_handle(a0)),
    point(vb0), point(edge_handle(b1)), point(vb1) };
  data(fa).bface.ctrl_pnt_idx = {
    data(va1).cpi,
    data(edge_handle(a2)).cpi, data(edge_handle(a0)).cpi,
    data(vb0).cpi, data(edge_handle(b1)).cpi, data(vb1).cpi };
  data(fa).bface.vert_idx = { va1.idx(), vb0.idx(), vb1.idx() };
  // fb: va1, b0, a1, vb1, b2, va0
  data(fb).bface.ctrl_pnt = {
    point(va1),
    point(edge_handle(b0)), point(edge_handle(a1)),
    point(vb1), point(edge_handle(b2)), point(va0) };
  data(fb).bface.ctrl_pnt_idx = {
    data(va1).cpi,
    data(edge_handle(b0)).cpi, data(edge_handle(a1)).cpi,
    data(vb1).cpi, data(edge_handle(b2)).cpi, data(va0).cpi };
  data(fb).bface.vert_idx = { va1.idx(), vb1.idx(), va0.idx() };
  return true;
}

void BMeshT::set_point(VertexHandle vh, const Vec3d& p)
{
  ctrl_pnt[data(vh).cpi] = p;
  for (FaceHandle fh : vf_range(vh))
  {
    BezierFace& bface = data(fh).bface;
    bface.global_corner_ctrlpnt(vh.idx()) = p;
  }
}

void BMeshT::set_point(EdgeHandle eh, const Vec3d& p)
{
  ctrl_pnt[data(eh).cpi] = p;
  HalfedgeHandle hh = halfedge_handle(eh, 0);
  if (!is_boundary(hh))
  {
    auto [v_from, v_to] = from_to_vh(*this, hh);
    BezierFace& bface = data(face_handle(hh)).bface;
    for (auto [cpi, cp] : bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()))
      cp = p;
    return;
  }
  hh = opposite_halfedge_handle(hh);
  if (!is_boundary(hh))
  {
    auto [v_from, v_to] = from_to_vh(*this, hh);
    BezierFace& bface = data(face_handle(hh)).bface;
    for (auto [cpi, cp] : bface.global_edge_interior_ctrlpnt_ms(v_from.idx(), v_to.idx()))
      cp = p;
    return;
  }
  ASSERT(false, "isolated edge.");
}


/// @brief Calculate the tangent vector along the halfedge hh,
/// starting from from_vertex_handle(hh).
Vec3d BMeshT::tangent_vec(HalfedgeHandle hh)const
{
  auto [v_from, v_to] = from_to_vh(*this, hh);
  if (!is_boundary(hh))
  {
    const BezierFace& bface = data(face_handle(hh)).bface;
    auto edge_ctrlpnts = bface.global_edge_ctrlpnt_cs(v_from.idx(), v_to.idx());
    auto iter = edge_ctrlpnts.begin();
    auto [first_i, first] = *iter;
    iter++;
    auto [second_i, second] = *iter;
    return second - first;
  }
  else
  {
    HalfedgeHandle opp_hh = opposite_halfedge_handle(hh);
    const BezierFace& bface = data(face_handle(opp_hh)).bface;
    auto edge_ctrlpnts = bface.global_edge_ctrlpnt_cr(v_to.idx(), v_from.idx());
    auto iter = edge_ctrlpnts.begin();
    auto [first_i, first] = *iter;
    iter++;
    auto [second_i, second] = *iter;
    return second - first;
  }
}

/// @brief Get the Bezier curve on edge eh.
/// @details The control points are in the order from v_from to v_to.
/// v_from = from_vertex_handle(halfedge_handle(eh, 0));
/// v_to = to_vertex_handle(halfedge_handle(eh, 0));
/// @param eh 
BezierCurveImpl BMeshT::curve_on_edge(EdgeHandle eh)const
{
  BezierCurveImpl c;

  HalfedgeHandle heh = halfedge_handle(eh, 0);
  VertexHandle v_from = from_vertex_handle(heh);
  VertexHandle v_to = to_vertex_handle(heh);

  if (v_from.idx() > v_to.idx())
  {
    heh = opposite_halfedge_handle(heh);
    std::swap(v_from, v_to);
  }
  // The control points are sorted from small index to large index
  // despite of the adjacent faces.
  if (is_boundary(heh))
  {
    heh = opposite_halfedge_handle(heh);
    ASSERT(!is_boundary(heh), "an isolated edge in Bezier mesh.");
    FaceHandle fh = face_handle(heh);
    const BezierFace& bface = data(fh).bface;
    c.degree = bface.degree;
    c.ctrl_points.reserve(c.degree + 1);
    for (auto [cpi, cp] : bface.global_edge_ctrlpnt_cr(v_to.idx(), v_from.idx()))
      c.ctrl_points.push_back(cp);
  }
  else
  {
    FaceHandle fh = face_handle(heh);
    const BezierFace& bface = data(fh).bface;
    c.degree = bface.degree;
    c.ctrl_points.reserve(c.degree + 1);
    for (auto [cpi, cp] : bface.global_edge_ctrlpnt_cs(v_from.idx(), v_to.idx()))
      c.ctrl_points.push_back(cp);
  }
  return c;
}

/// @brief Extract a local region to form a mesh.
/// @param [in] local_faces Faces contained by the local region.
/// @param [out] local The local mesh.
/// @param [out] v2v map global vertices to local vertices. 
/// @param [out] f2f map global faces to local faces.
/// @param [out] rv2v map local vertices to global vertices (reversed v2v).
/// @param [out] rf2f map local faces to global faces (reversed f2f).
void BMeshT::local_mesh(
  const std::set<FaceHandle>& local_faces, BMeshT& local,
  std::map<VertexHandle, VertexHandle>& v2v, std::map<FaceHandle, FaceHandle>& f2f,
  std::vector<VertexHandle>& rv2v, std::vector<FaceHandle>& rf2f)const
{
  // local.bezier_degree = bezier_degree;
  // add vertices to local mesh.
  std::set<VertexHandle> local_vertices;
  vertices_around_faces(this, local_faces, local_vertices);
  rv2v.resize(local_vertices.size());
  for (VertexHandle vh : local_vertices)
  {
    VertexHandle local_vh = local.add_vertex();
    v2v[vh] = local_vh;
    rv2v[local_vh.idx()] = vh;
  }

  // add faces and copy geometry to local mesh.
  std::vector<Index> ctrl_pnt_added(ctrl_pnt.size(), -1);
  rf2f.resize(local_faces.size());
  for (FaceHandle fh : local_faces)
  {
    auto [v0, v1, v2] = face_vertices(*this, fh);
    FaceHandle local_fh = local.add_face(v2v.at(v0), v2v.at(v1), v2v.at(v2));
    f2f[fh] = local_fh;
    rf2f[local_fh.idx()] = fh;

    const BezierFace& bface = data(fh).bface;
    BezierFace& local_bface = local.data(local_fh).bface;
    local_bface.degree = bface.degree;
    local_bface.ctrl_pnt = bface.ctrl_pnt;
    local_bface.vert_idx = {
      v2v.at(VertexHandle(bface.vert_idx[0])).idx(),
      v2v.at(VertexHandle(bface.vert_idx[1])).idx(),
      v2v.at(VertexHandle(bface.vert_idx[2])).idx()
    };
  }
  local.bezier_degree = bezier_degree;
  local.build_ctrlpnt_idx();
}

/// @brief linear approximation of all Bezier curves on edges.
/// @param density approximation density.
/// @return a set of polylines.
std::vector<std::vector<Vec3d>> BMeshT::linear_approx_curves(uint32_t density)const
{
  std::vector<Points> curves;

  for (EdgeHandle eh : edges())
  {
    BezierCurveImpl curve = curve_on_edge(eh);
    curves.push_back(BezierCurve::linear_approx(curve.ctrl_points, curve.degree, density));
  }

  return curves;
}

bool BMeshT::read_ascii(std::string file_name)
{
  std::fstream fin;
  fin.open(file_name, std::fstream::in);
  if (!fin.is_open())
    return false;

  clean();

  std::string line;
  size_t n_v, n_f, n_cp;

  // get vertices
  if (std::getline(fin, line))
  {
    line = line.substr(line.find(" "));
    n_v = std::stol(line);
  }
  else return false;
  // get faces
  if (std::getline(fin, line))
  {
    line = line.substr(line.find(" "));
    n_f = std::stol(line);
  }
  else return false;
  // get control points
  if (std::getline(fin, line))
  {
    line = line.substr(line.find(" "));
    n_cp = std::stol(line);
  }
  else return false;
  // get degree
  if (std::getline(fin, line))
  {
    line = line.substr(line.find(" "));
    bezier_degree = std::stol(line);
  }
  else return false;

  // Add vertices
  for (size_t i = 0;i < n_v;i++)
    add_vertex();

  auto split_string = [&](const std::string& line, char delim) {
    std::size_t previous = 0;
    std::size_t current = line.find(delim);
    std::vector<std::string> elems;
    while (current != std::string::npos)
    {
      if (current > previous)
      {
        elems.push_back(line.substr(previous, current - previous));
      }
      previous = current + 1;
      current = line.find(delim, previous);
    }
    if (previous != line.size())
    {
      elems.push_back(line.substr(previous));
    }
    return elems;
  };
  // get geometry of control points
  // get topology of mesh
  ctrl_pnt.reserve(n_cp);
  ctrl_pnt_deleted.resize(n_cp, false);
  while (std::getline(fin, line))
  {
    switch (line[0])
    {
    case 'c':
    {
      auto sub_strs = split_string(line, ' ');
      double x = std::stod(sub_strs[1]);
      double y = std::stod(sub_strs[2]);
      double z = std::stod(sub_strs[3]);
      ctrl_pnt.emplace_back(x, y, z);
    }
    break;
    case 'f':
    {
      // get face vertices
      auto sub_strs = split_string(line, ' ');
      Index v0 = std::stoi(sub_strs[1]);
      Index v1 = std::stoi(sub_strs[2]);
      Index v2 = std::stoi(sub_strs[3]);
      FaceHandle fh = add_face(VertexHandle(v0), VertexHandle(v1), VertexHandle(v2));
      // get face control points
      std::getline(fin, line);  // control point indices.
      sub_strs = split_string(line, ' ');
      BezierFace& bface = data(fh).bface;
      bface.degree = bezier_degree;
      bface.ctrl_pnt_idx.reserve(sub_strs.size());
      for (const auto& idx : sub_strs)
        bface.ctrl_pnt_idx.push_back(std::stoi(idx));
      bface.ctrl_pnt.resize(bface.ctrl_pnt_idx.size());
      bface.vert_idx = { v0,v1,v2 };
    }
    break;
    }
  }

  if (n_faces() != n_f)
    return false;
  if (ctrl_pnt.size() != n_cp)
    return false;

  sync_local_ctrlpnt();
  store_extra_cpi_to_ve();

  fin.close();
  return true;
}

/// @brief Write mesh into file in ascii mode.
/// Make sure that garbage is collected.
bool BMeshT::write_ascii(std::string file_name)
{
  std::fstream fout;
  fout.open(file_name, std::fstream::out);
  if (!fout.is_open())
    return false;

  // write size information.
  fout << "vertices " << n_vertices() << std::endl;
  fout << "faces " << n_faces() << std::endl;
  fout << "ctrlpnts " << ctrl_pnt.size() << std::endl;
  fout << "degree " << bezier_degree << std::endl;

  fout.precision(15);
  // write geometry information
  for (const Vec3d& cp : ctrl_pnt)
    fout << "c " << cp.x() << " " << cp.y() << " " << cp.z() << std::endl;

  // write topology information
  for (FaceHandle fh : faces())
  {
    BezierFace& bface = data(fh).bface;
    fout << "f " << bface.vert_idx[0] << " " << bface.vert_idx[1] << " " << bface.vert_idx[2] << std::endl;
    for (Index cpi : bface.ctrl_pnt_idx)
      fout << cpi << " ";
    fout << std::endl;
  }

  fout.close();
  return true;
}

void BMeshT::scale_to(ImageT& image)
{
  // calculate the bounding box for mesh
  BoundingBox box;
  for (VertexHandle vh : vertices())
  {
    box += point(vh);
  }
  box.min().z() = 0;
  box.max().z() = 0;

  // ASSERT(box.min().x() == 0. && box.min().y() == 0., "mesh does not start at (0, 0).");
  // For Xiao Yanyang's meshes, need to translation
  Vec3d translation = -box.min();

  int width = image.width;
  int height = image.height;

  double mesh_width = box.max().x();
  double mesh_height = box.max().y();

  double scale_width = (double)width/ mesh_width ;
  double scale_height = (double)height / mesh_height;

  for (VertexHandle vh : vertices())
  {
    Vec3d p = point(vh);
    p += translation;
    p.x() *= scale_width;
    p.y() *= scale_height;
    set_point(vh, p);
  }

  for (EdgeHandle eh : edges())
  {
    Vec3d p = point(eh);
    p += translation;
    p.x() *= scale_width;
    p.y() *= scale_height;
    set_point(eh, p);
  }

  // project boundary vertices to integer.
  for (VertexHandle vh : vertices())
  {
    if (!is_boundary(vh))
      continue;
    HalfedgeHandle hh = halfedge_handle(vh);
    Vec3d cur_p = point(vh);
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
    set_point(vh, cur_p);
  }
}

void BMeshT::delete_face(FaceHandle fh)
{
  for (VertexHandle fvh : fv_range(fh))
  {
    if (valence(fvh) != 2)
      continue;
    // fvh will be an isolated vertex.
    // remove its control point.
    delete_ctrl_pnt(data(fvh).cpi);
  }

  for (HalfedgeHandle hh : fh_range(fh))
  {
    if (!is_boundary(opposite_halfedge_handle(hh)))
      continue;
    // hh will be an isolated edge.
    // remove its interior control point.
  #if 0
    VertexHandle from_v = from_vertex_handle(hh), to_v = to_vertex_handle(hh);
    for (auto [cpi, cp] : data(fh).bface.global_edge_interior_ctrlpnt_cs(from_v.idx(), to_v.idx()))
      delete_ctrl_pnt(cpi);
  #endif
    // when degree is two, we have a more efficient implementation.
    delete_ctrl_pnt(data(edge_handle(hh)).cpi);
  }

  // when degree is two, there it no interior control point.
#if 0
  for (auto [cpi, cp] : data(fh).bface.global_interior_ctrlpnt_cs())
    delete_ctrl_pnt(cpi);
#endif

  TriConnectivity::delete_face(fh, true);
}

void BMeshT::garbage_collection()
{
  // Collect garbage of topology.
  // Trace the updated vertex handles.
  std::vector<VertexHandle> vhs(n_vertices());
  for (Index i = 0;i < n_vertices();i++)
    vhs[i] = VertexHandle(i);
  std::vector<VertexHandle*> vertices_to_update(n_vertices());
  for (Index i = 0;i < n_vertices();i++)
    vertices_to_update[i] = &vhs[i];

  std::vector<HalfedgeHandle*> empty_hh;
  std::vector<FaceHandle*> empty_fh;
  TriConnectivity::garbage_collection(vertices_to_update, empty_hh, empty_fh);

  // Collect garbage of geometry.
  // Trace the updated indices of control point.
  std::vector<Index> cp_map(ctrl_pnt.size());
  std::iota(cp_map.begin(), cp_map.end(), 0);

  if (!ctrl_pnt.empty())
  {
    Index i0 = 0, i1 = (Index)ctrl_pnt.size() - 1;
    // really delete control points.
    while (true)
    {
      // find 1st deleted and last un-deleted
      while (!ctrl_pnt_deleted[i0] && i0 < i1)  ++i0;
      while (ctrl_pnt_deleted[i1] && i0 < i1)  --i1;
      if (i0 >= i1) break;

      // We keep track of the control points for updates,
      // we need to have the opposite direction
      cp_map[i1] = i0;
      cp_map[i0] = -1;

      // swap all properties about i0 and i1
      std::swap(ctrl_pnt[i0], ctrl_pnt[i1]);
      std::swap(ctrl_pnt_deleted[i0], ctrl_pnt_deleted[i1]);
    }

    // resize all properties
    ctrl_pnt.resize(ctrl_pnt_deleted[i0] ? i0 : i0 + 1);
    ctrl_pnt_deleted.resize(ctrl_pnt.size());
  }

  // Update the indices in BezierFace::vert_idx and BezierFace::ctrl_pnt_idx.
  // Update the control points in BezierFace::ctrl_pnt.
  for (FaceHandle fh : faces())
  {
    BezierFace& bface = data(fh).bface;
    for (Index i = 0;i < 3;i++)
      bface.vert_idx[i] = vhs[bface.vert_idx[i]].idx();

    for (Index i = 0;i < bface.ctrl_pnt_idx.size();i++)
    {
      bface.ctrl_pnt_idx[i] = cp_map[bface.ctrl_pnt_idx[i]];
      bface.ctrl_pnt[i] = ctrl_pnt[bface.ctrl_pnt_idx[i]];
    }
  }

  store_extra_cpi_to_ve();
}

void BMeshT::clean()
{
  TriConnectivity::resize(0, 0, 0);
  ctrl_pnt.clear();
  ctrl_pnt_deleted.clear();
}

void BMeshT::copy_from(const BMeshT& rhs)
{
  bezier_degree = rhs.bezier_degree;
  ctrl_pnt = rhs.ctrl_pnt;
  ctrl_pnt_deleted = rhs.ctrl_pnt_deleted;
  data_vpph_ = rhs.data_vpph_;
  data_hpph_ = rhs.data_hpph_;
  data_epph_ = rhs.data_epph_;
  data_fpph_ = rhs.data_fpph_;
}

void BMeshT::move_from(BMeshT&& rhs)
{
  bezier_degree = rhs.bezier_degree;
  ctrl_pnt = std::move(rhs.ctrl_pnt);
  ctrl_pnt_deleted = std::move(rhs.ctrl_pnt_deleted);
  data_vpph_ = rhs.data_vpph_;
  data_hpph_ = rhs.data_hpph_;
  data_epph_ = rhs.data_epph_;
  data_fpph_ = rhs.data_fpph_;
}

}// namespace ImageTriSimp
}// namespace GCLF