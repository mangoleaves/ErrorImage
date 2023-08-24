#include "LEdgeCollapser.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

LEdgeCollapser::LEdgeCollapser(LinearSimplifier* _simplifier)
  :LLocalOperation(_simplifier)
{}

/// @brief Initialize neccessary data before collapsing.
/// @return true of collapse ok.
bool LEdgeCollapser::init(EdgeHandle _e, size_t Np)
{
  clear();
  collapse_he = mesh->halfedge_handle(_e, 0);
  collapse_he_opp = mesh->halfedge_handle(_e, 1);
  from_vh = mesh->from_vertex_handle(collapse_he);
  to_vh = mesh->to_vertex_handle(collapse_he);
  // can't collapse
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

void LEdgeCollapser::clear()
{
  LLocalOperation::clear();
  collapse_he.invalidate();
  collapse_he_opp.invalidate();
  halfedges.clear();
}

uint32_t LEdgeCollapser::valence_after_collapse()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return local_mesh->valence(local_optimizing_vertex);
}

VertexHandle LEdgeCollapser::collapse_vertex()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return to_vh;
}

/// @brief Find affected faces on mesh and affected pixels on image.
void LEdgeCollapser::find_affected_faces_pixels()
{
  // Initialize halfedges.
  int n_valance = mesh->valence(mesh->from_vertex_handle(collapse_he)) +
    mesh->valence(mesh->to_vertex_handle(collapse_he));
  halfedges.clear();
  halfedges.reserve(n_valance);

  HalfedgeHandle hh_next = mesh->next_halfedge_handle(collapse_he);
  HalfedgeHandle hh_oppo_next = mesh->next_halfedge_handle(collapse_he_opp);

  HalfedgeHandle moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(collapse_he));
  if (mesh->is_boundary(moving_hh) && moving_hh != hh_oppo_next)
    moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(moving_hh));
  while (moving_hh != hh_oppo_next)
  {
    moving_hh = mesh->next_halfedge_handle(moving_hh);
    if (!mesh->is_boundary(moving_hh))
      halfedges.push_back(moving_hh);
    moving_hh = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(moving_hh));
    if (mesh->is_boundary(moving_hh) && moving_hh != hh_oppo_next)
      moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(moving_hh));
  }

  moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(collapse_he_opp));
  if (mesh->is_boundary(moving_hh) && moving_hh != hh_next)
    moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(moving_hh));
  while (moving_hh != hh_next)
  {
    moving_hh = mesh->next_halfedge_handle(moving_hh);
    if (!mesh->is_boundary(moving_hh))
      halfedges.push_back(moving_hh);
    moving_hh = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(moving_hh));
    if (mesh->is_boundary(moving_hh) && moving_hh != hh_next)
      moving_hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(moving_hh));
  }
  // Initialize one ring faces
  for (HalfedgeHandle hh : halfedges)
    global_affected_faces.insert(mesh->face_handle(hh));
  if (!mesh->is_boundary(collapse_he))
    global_affected_faces.insert(mesh->face_handle(collapse_he));
  if (!mesh->is_boundary(collapse_he_opp))
    global_affected_faces.insert(mesh->face_handle(collapse_he_opp));
  // Initialize pixels
  for (FaceHandle fh : global_affected_faces)
  {
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  if (!mesh->is_boundary(collapse_he))
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(mesh->face_handle(collapse_he)).pixels.begin(),
      mesh->data(mesh->face_handle(collapse_he)).pixels.end());
  if (!mesh->is_boundary(collapse_he_opp))
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(mesh->face_handle(collapse_he_opp)).pixels.begin(),
      mesh->data(mesh->face_handle(collapse_he_opp)).pixels.end());

  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void LEdgeCollapser::initialize_local_mesh()
{
  local_mesh = std::make_unique<LMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Collapse the local mesh
  HalfedgeHandle local_collapse_he = local_mesh->find_halfedge(v2v.at(from_vh), v2v.at(to_vh));
  local_optimizing_vertex = local_mesh->to_vertex_handle(local_collapse_he);
  local_mesh->collapse(local_collapse_he);  // the edge is collapsed to v2v.at(to_vh).
}

void LEdgeCollapser::update_global_mesh(const Vec3d& new_point)
{
  // update global mesh
  mesh->collapse(collapse_he);
  mesh->set_point(to_vh, new_point);
  for (FaceHandle local_fh : local_mesh->vf_range(local_optimizing_vertex))
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