#include "LBoundaryEdgeCollapser.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

LBoundaryEdgeCollapser::LBoundaryEdgeCollapser(LinearSimplifier* _simplifier)
  :LLocalOperation(_simplifier)
{}

bool LBoundaryEdgeCollapser::init(EdgeHandle _e, size_t Np)
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

void LBoundaryEdgeCollapser::clear()
{
  LLocalOperation::clear();
  initialized = false;
  collapse_he.invalidate();
  collapse_he_opp.invalidate();
  from_vh.invalidate();
  to_vh.invalidate();
  forward_end_vh.invalidate();
  backward_end_vh.invalidate();
  to_vh_is_corner = false;
}

uint32_t LBoundaryEdgeCollapser::valence_after_collapse()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return local_mesh->valence(local_to_vh);
}

bool LBoundaryEdgeCollapser::is_collapsed_to_corner()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return to_vh_is_corner;
}

Vec3d LBoundaryEdgeCollapser::corner_position()
{
  ASSERT(to_vh_is_corner, "Not a corner.");
  return mesh->point(to_vh);
}

VertexHandle LBoundaryEdgeCollapser::collapse_vertex()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return to_vh;
}

bool LBoundaryEdgeCollapser::check_constraints(const LMeshT& new_local_mesh)
{
  return LLocalOperation::check_constraints(new_local_mesh);
}

bool LBoundaryEdgeCollapser::check_constraints()
{
  ASSERT(initialized, "Collapser is not initialized.");
  ASSERT(to_vh_is_corner, "Do not collapse to corner.");
  if (operation_would_cause_flip(*local_mesh))
    return false;
  return true;
}

double LBoundaryEdgeCollapser::error_after(size_t pidx)
{
  return LLocalOperation::error_after(pidx);
}

double LBoundaryEdgeCollapser::error_after()
{
  ASSERT(initialized, "Collapser is not initialized.");
  ASSERT(to_vh_is_corner, "Do not collapse to corner.");

  update_pixel_color_error(*local_mesh);

  double color_error = std::accumulate(
    local_mesh->faces_begin(), local_mesh->faces_end(), 0.0,
    [&](double v, FaceHandle fh) {return v + local_mesh->data(fh).error;});
  return color_error;
}

double LBoundaryEdgeCollapser::quality_after(size_t pidx)
{
  return LLocalOperation::quality_after(pidx);
}

double LBoundaryEdgeCollapser::quality_after()
{
  ASSERT(initialized, "Collapser is not initialized.");
  ASSERT(to_vh_is_corner, "Do not collapse to corner.");

  return calc_max_cos(local_mesh.get());
}

void LBoundaryEdgeCollapser::just_do_it(const Vec3d& new_point)
{
  if (to_vh_is_corner)
  {
    update_pixel_color_error(*local_mesh);
    update_global_mesh(Vec3d());
  }
  else
  {
    LLocalOperation::just_do_it(new_point);
  }
}

/// @brief Find affected faces on mesh and affected pixels on image.
void LBoundaryEdgeCollapser::find_affected_faces_pixels()
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

void LBoundaryEdgeCollapser::initialize_local_mesh()
{
  local_mesh = std::make_unique<LMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  local_from_vh = v2v.at(from_vh);
  local_to_vh = v2v.at(to_vh);
  Vec3d to_point = local_mesh->point(local_to_vh);

  HalfedgeHandle hh = local_mesh->find_halfedge(local_from_vh, local_to_vh);
  local_mesh->collapse(hh);

  local_forward_end_vh = v2v.at(forward_end_vh);
  local_backward_end_vh = v2v.at(backward_end_vh);

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

  local_optimizing_vertex = v2v.at(to_vh);
}

void LBoundaryEdgeCollapser::update_local_mesh(
  LMeshT& mesh_to_update,
  const Vec3d& new_point)
{
  ASSERT(!to_vh_is_corner, "relocate a corner vertex.");
  mesh_to_update.set_point(local_optimizing_vertex, new_point);
  // project boundary vertex
  if (do_project_x)
    mesh_to_update.point(local_optimizing_vertex).x() = x_project_to;
  if (do_project_y)
    mesh_to_update.point(local_optimizing_vertex).y() = y_project_to;
}

void LBoundaryEdgeCollapser::update_global_mesh(const Vec3d& new_point)
{
  // update global mesh
  // do real collapse operation.
  mesh->collapse(collapse_he);

  if (!to_vh_is_corner)
  {
    mesh->set_point(to_vh, new_point);
    if (do_project_x)
      mesh->point(to_vh).x() = x_project_to;
    if (do_project_y)
      mesh->point(to_vh).y() = y_project_to;
  }
  // update control points in global mesh.
  for (FaceHandle local_fh : local_mesh->vf_range(local_to_vh))
  {
    FaceHandle global_fh = rf2f[local_fh.idx()];
    // update pixels, color and error
    mesh->data(global_fh).pixels = local_mesh->data(local_fh).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fh).color;
    mesh->data(global_fh).error = local_mesh->data(local_fh).error;
  }
}

}// namespace ImageTriSimp
}// namespace GCLF