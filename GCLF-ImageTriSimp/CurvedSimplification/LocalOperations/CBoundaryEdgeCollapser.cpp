#include "CBoundaryEdgeCollapser.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

CBoundaryEdgeCollapser::CBoundaryEdgeCollapser(CurvedSimplifier* _simplifier)
  :CLocalOperation(_simplifier)
{
}

bool CBoundaryEdgeCollapser::init(EdgeHandle _e, size_t Np)
{
  clear();
  // Initialize handles
  if (!mesh->is_boundary(_e))
    return false;
  collapse_he = mesh->halfedge_handle(_e, 0);
  collapse_he_opp = mesh->opposite_halfedge_handle(collapse_he);
  from_vh = mesh->from_vertex_handle(collapse_he);
  to_vh = mesh->to_vertex_handle(collapse_he);
  // can't collapse edge
  if (!mesh->is_collapse_ok(collapse_he))
    return false;

  to_vh_is_corner = is_corner_vertex(to_vh);
  bool from_vh_is_corner = is_corner_vertex(from_vh);
  if (from_vh_is_corner && to_vh_is_corner)
  {
    return false;
    if (is_corner_vertex(from_vh))
      return false;
  }
  else if(from_vh_is_corner)
  {
    std::swap(collapse_he, collapse_he_opp);
    std::swap(from_vh, to_vh);
    std::swap(from_vh_is_corner, to_vh_is_corner);
  }

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

void CBoundaryEdgeCollapser::clear()
{
  CLocalOperation::clear();
  initialized = false;
  collapse_he.invalidate();
  collapse_he_opp.invalidate();
  from_vh.invalidate();
  to_vh.invalidate();
  forward_end_vh.invalidate();
  backward_end_vh.invalidate();
}

uint32_t CBoundaryEdgeCollapser::valence_after_collapse()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return local_mesh->valence(local_to_vh);
}

VertexHandle CBoundaryEdgeCollapser::collapse_vertex()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return to_vh;
}

std::vector<Points> CBoundaryEdgeCollapser::get_init_population(
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
void CBoundaryEdgeCollapser::find_affected_faces_pixels()
{
  // pixels on one ring faces of boundary_vh and not_boundary_vh
  // will be re-assigned.
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  for (HalfedgeHandle voh : mesh->voh_range(from_vh))
  {
    if (!mesh->is_boundary(voh))
    {
      FaceHandle fh = mesh->face_handle(voh);
      global_affected_faces.insert(fh);
    }
  }
  for (HalfedgeHandle voh : mesh->voh_range(to_vh))
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

  n_faces_after_collapsing = global_affected_faces.size() - 1;

  HalfedgeHandle forward_hh = mesh->halfedge_handle(to_vh);
  ASSERT(mesh->is_boundary(forward_hh) && mesh->from_vertex_handle(forward_hh) == to_vh,
    "forward halfedge error.");
  if (mesh->to_vertex_handle(forward_hh) != from_vh)
    forward_end_vh = mesh->to_vertex_handle(forward_hh);
  else
    forward_end_vh = mesh->from_vertex_handle(mesh->prev_halfedge_handle(forward_hh));
  forward_end_p = mesh->point(forward_end_vh);

  HalfedgeHandle backward_hh = mesh->prev_halfedge_handle(mesh->halfedge_handle(from_vh));
  ASSERT(mesh->is_boundary(backward_hh) && mesh->to_vertex_handle(backward_hh) == from_vh,
    "backward halfedge error.");
  if (mesh->from_vertex_handle(backward_hh) != to_vh)
    backward_end_vh = mesh->from_vertex_handle(backward_hh);
  else
    backward_end_vh = mesh->to_vertex_handle(mesh->next_halfedge_handle(backward_hh));
  backward_end_p = mesh->point(backward_end_vh);
}

void CBoundaryEdgeCollapser::initialize_local_mesh()
{
  local_mesh = std::make_unique<BMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  local_from_vh = v2v.at(from_vh);
  local_to_vh = v2v.at(to_vh);
  Vec3d to_point = local_mesh->point(local_to_vh);

  HalfedgeHandle hh = local_mesh->find_halfedge(local_from_vh, local_to_vh);
  local_mesh->collapse(hh);

  local_to_cpi = local_mesh->data(local_to_vh).cpi;
  local_forward_end_vh = v2v.at(forward_end_vh);
  local_backward_end_vh = v2v.at(backward_end_vh);
  {
    hh = local_mesh->find_halfedge(local_to_vh, local_forward_end_vh);
    local_forward_cpi = local_mesh->data(local_mesh->edge_handle(hh)).cpi;
    hh = local_mesh->find_halfedge(local_to_vh, local_backward_end_vh);
    local_backward_cpi = local_mesh->data(local_mesh->edge_handle(hh)).cpi;
  }

  // Find control points to optimize.
  if (to_vh_is_corner)
  {
    do_project_x = true;
    x_project_to = to_point.x();
    do_project_y = true;
    y_project_to = to_point.y();
  }
  else
  {
    const Vec3d& boundary_point = local_mesh->point(local_to_vh);
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

  local_optimizing_cpi.push_back(local_mesh->data(local_to_vh).cpi);
  for (HalfedgeHandle hh : local_mesh->voh_range(local_to_vh))
  {
    VertexHandle to_vh = local_mesh->to_vertex_handle(hh);
    if (to_vh != local_forward_end_vh && to_vh != local_backward_end_vh)
    {
      local_optimizing_cpi.push_back(local_mesh->data(local_mesh->edge_handle(hh)).cpi);
      local_affected_edges.push_back(local_mesh->edge_handle(hh));
    }
  }
}

void CBoundaryEdgeCollapser::update_local_ctrlpnts(
  BMeshT& mesh_to_update,
  const Points& new_ctrlpnts)
{
  for (size_t i = 0;i < local_optimizing_cpi.size();i++)
    mesh_to_update.ctrl_pnt[local_optimizing_cpi[i]] = new_ctrlpnts[i];
  // project boundary vertex
  if (do_project_x)
    mesh_to_update.ctrl_pnt[local_to_cpi].x() = x_project_to;
  if (do_project_y)
    mesh_to_update.ctrl_pnt[local_to_cpi].y() = y_project_to;
  // fix two boundary control points
  mesh_to_update.ctrl_pnt[local_forward_cpi] =
    0.5 * (mesh_to_update.ctrl_pnt[local_to_cpi] + forward_end_p);
  mesh_to_update.ctrl_pnt[local_backward_cpi] =
    0.5 * (mesh_to_update.ctrl_pnt[local_to_cpi] + backward_end_p);

  mesh_to_update.sync_local_ctrlpnt();
}

void CBoundaryEdgeCollapser::update_global_mesh(const Points& new_ctrlpnts)
{
  // update global mesh
  // do real collapse operation.
  mesh->collapse(collapse_he);

  // update control points in global mesh.
  for (FaceHandle local_fh : local_mesh->vf_range(local_to_vh))
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