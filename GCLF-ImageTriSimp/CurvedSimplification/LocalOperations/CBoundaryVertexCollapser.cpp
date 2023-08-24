#include "CBoundaryVertexCollapser.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

CBoundaryVertexCollapser::CBoundaryVertexCollapser(CurvedSimplifier* _simplifier)
  :CLocalOperation(_simplifier)
{}

bool CBoundaryVertexCollapser::init(EdgeHandle _e, size_t Np)
{
  clear();
  // Initialize handles
  collapse_he = mesh->halfedge_handle(_e, 0);
  VertexHandle from_vh = mesh->from_vertex_handle(collapse_he);
  VertexHandle to_vh = mesh->to_vertex_handle(collapse_he);
  // check if its the boundary case 
  if (!mesh->is_boundary(_e) && (mesh->is_boundary(from_vh) || mesh->is_boundary(to_vh)))
  {
    if (mesh->is_boundary(from_vh))
    {
      collapse_he_opp = collapse_he;
      collapse_he = mesh->opposite_halfedge_handle(collapse_he);
      boundary_vh = from_vh;
      not_boundary_vh = to_vh;
    }
    else
    {
      collapse_he_opp = mesh->opposite_halfedge_handle(collapse_he);
      boundary_vh = to_vh;
      not_boundary_vh = from_vh;
    }
  }
  else return false;
  // can't collapse edge
  if (!mesh->is_collapse_ok(collapse_he))
    return false;

  // Initialize neccessary topology data.
  find_affected_faces_pixels();

  // Initialize neccessary geometry data.
  initialize_local_mesh();

  // Initialize virtual local mesh used to virtually update.
  virtual_update_struct tmp_vus;
  tmp_vus.vt_local_mesh = *local_mesh;
  tmp_vus.vt_pixel_assigned = pixel_assigned;
  tmp_vus.vt_updated = false;
  vus.resize(Np, tmp_vus);

  initialized = true;
  return true;
}

void CBoundaryVertexCollapser::clear()
{
  CLocalOperation::clear();
  collapse_he.invalidate();
  collapse_he_opp.invalidate();
  boundary_vh.invalidate();
  not_boundary_vh.invalidate();
  faces_to_delete.clear();
  forward_end_vh.invalidate();
  backward_end_vh.invalidate();
}

uint32_t CBoundaryVertexCollapser::valence_after_collapse()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return local_mesh->valence(local_boundary_vh);
}

VertexHandle CBoundaryVertexCollapser::collapse_vertex()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return boundary_vh;
}

std::vector<Points> CBoundaryVertexCollapser::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  BoundingBox bbox = bounding_box();
  double radius = (bbox.max() - bbox.min()).length() * 0.1;
  std::vector<Points> population; population.resize(Np, get_variables());
  for (Points& ps : population)
  {
    for (Vec3d& p : ps)
    {
      // TODO: consider nearby feature points
      p.x() += (init_dis(init_generator) - 0.5) * radius * 2.0;
      p.y() += (init_dis(init_generator) - 0.5) * radius * 2.0;
    }
  }
  return population;
}

/// @brief Find affected faces on mesh and affected pixels on image.
void CBoundaryVertexCollapser::find_affected_faces_pixels()
{
  // pixels on one ring faces of boundary_vh and not_boundary_vh
  // will be re-assigned.
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  for (HalfedgeHandle voh : mesh->voh_range(boundary_vh))
  {
    if (!mesh->is_boundary(voh))
    {
      FaceHandle fh = mesh->face_handle(voh);
      global_affected_faces.insert(fh);
    }
  }
  for (HalfedgeHandle voh : mesh->voh_range(not_boundary_vh))
  {
    if (!mesh->is_boundary(voh))
    {
      FaceHandle fh = mesh->face_handle(voh);
      global_affected_faces.insert(fh);
    }
  }
  for (FaceHandle fh : global_affected_faces)
  {
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;

  boundary_vh_is_corner = is_corner_vertex(boundary_vh);

  n_faces_after_collapsing = global_affected_faces.size();
  need_collapse = 2;

  /// we gonna find faces_to_delete.
  /// see figure notes made by Guo Jia-Peng.
  HalfedgeHandle moving_hh;
  moving_hh = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(collapse_he));
  if (mesh->is_boundary(moving_hh))
  {
    need_collapse -= 1;
    while (mesh->opposite_he_opposite_vh(moving_hh) == not_boundary_vh)
    {
      faces_to_delete.push_back(mesh->opposite_face_handle(moving_hh));
      backward_end_vh = mesh->from_vertex_handle(moving_hh);
      if (is_corner_vertex(backward_end_vh))
        break;
      moving_hh = mesh->prev_halfedge_handle(moving_hh);
    }
  }
  else
  {
    moving_hh = mesh->prev_halfedge_handle(mesh->halfedge_handle(boundary_vh));
    ASSERT(mesh->to_vertex_handle(moving_hh) == boundary_vh, "boundary halfedge error.");
    backward_end_vh = mesh->from_vertex_handle(moving_hh);
    n_faces_after_collapsing -= 1;
  }
  moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(collapse_he_opp));
  if (mesh->is_boundary(moving_hh))
  {
    need_collapse -= 1;
    while (mesh->opposite_he_opposite_vh(moving_hh) == not_boundary_vh)
    {
      faces_to_delete.push_back(mesh->opposite_face_handle(moving_hh));
      forward_end_vh = mesh->to_vertex_handle(moving_hh);
      if (is_corner_vertex(forward_end_vh))
        break;
      moving_hh = mesh->next_halfedge_handle(moving_hh);
    }
  }
  else
  {
    moving_hh = mesh->halfedge_handle(boundary_vh);
    ASSERT(mesh->from_vertex_handle(moving_hh) == boundary_vh, "boundary halfedge error.");
    forward_end_vh = mesh->to_vertex_handle(moving_hh);
    n_faces_after_collapsing -= 1;
  }
  forward_end_p = mesh->point(forward_end_vh);
  backward_end_p = mesh->point(backward_end_vh);
  n_faces_after_collapsing -= faces_to_delete.size();
}

void CBoundaryVertexCollapser::initialize_local_mesh()
{
  local_mesh = std::make_unique<BMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  local_boundary_vh = v2v.at(boundary_vh);
  local_not_boundary_vh = v2v.at(not_boundary_vh);
  Vec3d boundary_point = local_mesh->point(local_boundary_vh);

  // Collapse the local mesh
  for (FaceHandle fh : faces_to_delete)
    local_mesh->delete_face(f2f.at(fh));

  // No matter we do collapse or not, we need to logically
  // delete not_boundary_vh, collapse the edge to boundary_vh.
  // That is, not_boundary_vh is deleted, boundary_vh is keeped
  // as the final boundary vertex.
  if (need_collapse)
  {
    HalfedgeHandle local_collapse_he =
      local_mesh->find_halfedge(local_not_boundary_vh, local_boundary_vh);
    local_mesh->collapse(local_collapse_he);
    // after collapsing, not_boundary_vh is deleted, boundary_vh is keeped.
  }
  else
  {
    // if we only delete faces, boundary_vh is deleted, not_boundary_vh is keeped,
    // thus we need to swap them and update geometry.
    std::swap(local_boundary_vh, local_not_boundary_vh);
    local_mesh->set_point(local_boundary_vh, boundary_point);
  }

  local_boundary_cpi = local_mesh->data(local_boundary_vh).cpi;
  local_forward_end_vh = v2v.at(forward_end_vh);
  local_backward_end_vh = v2v.at(backward_end_vh);
  {
    HalfedgeHandle hh = local_mesh->halfedge_handle(local_boundary_vh);
    local_forward_cpi = local_mesh->data(local_mesh->edge_handle(hh)).cpi;
    ASSERT(local_mesh->to_vertex_handle(hh) == local_forward_end_vh, "local mesh topology error.");
    hh = local_mesh->prev_halfedge_handle(hh);
    ASSERT(local_mesh->from_vertex_handle(hh) == local_backward_end_vh, "local mesh topology error.");
    local_backward_cpi = local_mesh->data(local_mesh->edge_handle(hh)).cpi;
  }

  // Find control points to optimize.
  if (boundary_vh_is_corner)
  {
    do_project_x = true;
    x_project_to = boundary_point.x();
    do_project_y = true;
    y_project_to = boundary_point.y();
  }
  else
  {
    const Vec3d& boundary_point = local_mesh->point(local_boundary_vh);
    if (fabs(forward_end_p.x() - backward_end_p.x()) < 1e-6)
    {
      do_project_x = true;
      x_project_to = boundary_point.x();
      do_project_y = false;
    }
    else if (fabs(forward_end_p.y() - backward_end_p.y()) < 1e-6)
    {
      do_project_x = false;
      do_project_y = true;
      y_project_to = boundary_point.y();
    }
    else
    {
      ASSERT(false, "project error.");
    }
  }

  local_optimizing_cpi.push_back(local_mesh->data(local_boundary_vh).cpi);
  for (HalfedgeHandle hh : local_mesh->voh_range(local_boundary_vh))
  {
    VertexHandle to_vh = local_mesh->to_vertex_handle(hh);
    if (to_vh != local_forward_end_vh && to_vh != local_backward_end_vh)
    {
      local_optimizing_cpi.push_back(local_mesh->data(local_mesh->edge_handle(hh)).cpi);
      local_affected_edges.push_back(local_mesh->edge_handle(hh));
    }
  }
}

void CBoundaryVertexCollapser::update_local_ctrlpnts(
  BMeshT& mesh_to_update,
  const Points& new_ctrlpnts)
{
  for (size_t i = 0;i < local_optimizing_cpi.size();i++)
    mesh_to_update.ctrl_pnt[local_optimizing_cpi[i]] = new_ctrlpnts[i];
  // project boundary vertex
  if (do_project_x)
    mesh_to_update.ctrl_pnt[local_boundary_cpi].x() = x_project_to;
  if (do_project_y)
    mesh_to_update.ctrl_pnt[local_boundary_cpi].y() = y_project_to;
  // fix two boundary control points
  mesh_to_update.ctrl_pnt[local_forward_cpi] =
    0.5 * (mesh_to_update.ctrl_pnt[local_boundary_cpi] + forward_end_p);
  mesh_to_update.ctrl_pnt[local_backward_cpi] =
    0.5 * (mesh_to_update.ctrl_pnt[local_boundary_cpi] + backward_end_p);

  mesh_to_update.sync_local_ctrlpnt();
  // extra cpi is not stored in VertexData and EdgeData.
}

void CBoundaryVertexCollapser::update_global_mesh(const Points& new_ctrlpnts)
{
  // update global mesh
  // do real collapse operation.
  for (FaceHandle fh : faces_to_delete)
    mesh->delete_face(fh);
  if (need_collapse)
    mesh->collapse(collapse_he);
  // else std::swap(boundary_vh, not_boundary_vh);

  // update control points in global mesh.
  for (FaceHandle local_fh : local_mesh->vf_range(local_boundary_vh))
  {
    FaceHandle global_fh = rf2f[local_fh.idx()];
    // update control points
    mesh->data(global_fh).bface.ctrl_pnt = local_mesh->data(local_fh).bface.ctrl_pnt;
    for (auto [cp_idx, cp] : mesh->data(global_fh).bface.global_ctrlpnt_cs())
      mesh->ctrl_pnt[cp_idx] = cp;
    // update pixels, color and error
    mesh->data(global_fh).pixels = local_mesh->data(local_fh).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fh).color;
    mesh->data(global_fh).error = local_mesh->data(local_fh).error;
  }
}

}// namespace ImageTriSimp
}// namespace GCLF