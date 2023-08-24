#include "LEdgeFlipper.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

LEdgeFlipper::LEdgeFlipper(LinearSimplifier* _simplifier)
  :LLocalOperation(_simplifier)
{
}

/// @brief Initialize neccessary data before flipping.
/// @return true if flip_ok.
bool LEdgeFlipper::init(EdgeHandle _e)
{
  clear();
  // can't flip edge
  if (mesh->is_boundary(_e))
    return false;
  if (!mesh->is_flip_ok(_e))
    return false;

  flip_eh = _e;

  // Initialize neccessary topology data.
  find_affected_faces_pixels();

  // Initialize neccessary geometry data.
  initialize_local_mesh();

  initialized = true;
  return true;
}

void LEdgeFlipper::clear()
{
  LLocalOperation::clear();
  initialized = false;
}

bool LEdgeFlipper::check_constraints()const
{
  return LLocalOperation::check_constraints(*local_mesh);
}

bool LEdgeFlipper::flip_would_decrease_valence()
{
  int target_val_a0 = mesh->is_boundary(ghs.va0) ? 4 : 6;
  int target_val_a1 = mesh->is_boundary(ghs.va1) ? 4 : 6;
  int target_val_b0 = mesh->is_boundary(ghs.vb0) ? 4 : 6;
  int target_val_b1 = mesh->is_boundary(ghs.vb1) ? 4 : 6;

  int deviation_pre =
    std::abs((int)(mesh->valence(ghs.va0)) - target_val_a0) +
    std::abs((int)(mesh->valence(ghs.va1)) - target_val_a1) +
    std::abs((int)(mesh->valence(ghs.vb0)) - target_val_b0) +
    std::abs((int)(mesh->valence(ghs.vb1)) - target_val_b1);
  int deviation_post =
    std::abs((int)(mesh->valence(ghs.va0)) - 1 - target_val_a0) +
    std::abs((int)(mesh->valence(ghs.va1)) + 1 - target_val_a1) +
    std::abs((int)(mesh->valence(ghs.vb0)) - 1 - target_val_b0) +
    std::abs((int)(mesh->valence(ghs.vb1)) + 1 - target_val_b1);

  return deviation_post < deviation_pre;
}

bool LEdgeFlipper::flip_would_exceed_max_valence(uint32_t max_valence)
{
  return mesh->valence(ghs.va1) + 1 > max_valence ||
    mesh->valence(ghs.vb1) + 1 > max_valence;
}

double LEdgeFlipper::error_after()
{
  auto& behavior = simplifier->behavior;
  // assign pixels to new local mesh.
  try
  {
    behavior.assign_pixels_to_face(local_mesh.get(), image, pixel_assigned, affected_pixel_coords);
  }
  catch (...)
  {
    return DBL_MAX;
  }

  // calculate color and error on new local mesh.
  for (FaceHandle fh : local_mesh->faces())
  {
    local_mesh->data(fh).color = behavior.calc_color_on_face(
      local_mesh->data(fh).pixels);
    local_mesh->data(fh).error = behavior.calc_error_on_face(
      local_mesh->data(fh).color,
      local_mesh->data(fh).pixels);
  }

  double color_error = std::accumulate(
    local_mesh->faces_begin(), local_mesh->faces_end(), 0.0,
    [&](double v, FaceHandle fh) {return v + local_mesh->data(fh).error;});
  return color_error;
}

double LEdgeFlipper::quality_after()
{
  return calc_max_cos(local_mesh.get());
}

void LEdgeFlipper::just_do_it()
{
  update_global_mesh();
}

/// @brief Find affected faces on mesh and affected pixels on image.
void LEdgeFlipper::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();

  ghs.a0 = mesh->halfedge_handle(flip_eh, 0);
  ghs.b0 = mesh->halfedge_handle(flip_eh, 1);
  ghs.a1 = mesh->next_halfedge_handle(ghs.a0);
  ghs.a2 = mesh->next_halfedge_handle(ghs.a1);
  ghs.b1 = mesh->next_halfedge_handle(ghs.b0);
  ghs.b2 = mesh->next_halfedge_handle(ghs.b1);
  ghs.va0 = mesh->to_vertex_handle(ghs.a0);
  ghs.va1 = mesh->to_vertex_handle(ghs.a1);
  ghs.vb0 = mesh->to_vertex_handle(ghs.b0);
  ghs.vb1 = mesh->to_vertex_handle(ghs.b1);
  ghs.fa = mesh->face_handle(ghs.a0);
  ghs.fb = mesh->face_handle(ghs.b0);

  global_affected_faces.insert(ghs.fa);
  global_affected_faces.insert(ghs.fb);
  for (FaceHandle fh : global_affected_faces)
  {
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void LEdgeFlipper::initialize_local_mesh()
{
  local_mesh = std::make_unique<LMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Intialize local handles
  lhs.a0 = local_mesh->find_halfedge(v2v.at(ghs.vb0), v2v.at(ghs.va0));
  lhs.b0 = local_mesh->opposite_halfedge_handle(lhs.a0);
  lhs.a1 = local_mesh->next_halfedge_handle(lhs.a0);
  lhs.a2 = local_mesh->next_halfedge_handle(lhs.a1);
  lhs.b1 = local_mesh->next_halfedge_handle(lhs.b0);
  lhs.b2 = local_mesh->next_halfedge_handle(lhs.b1);
  lhs.va0 = local_mesh->to_vertex_handle(lhs.a0);
  lhs.va1 = local_mesh->to_vertex_handle(lhs.a1);
  lhs.vb0 = local_mesh->to_vertex_handle(lhs.b0);
  lhs.vb1 = local_mesh->to_vertex_handle(lhs.b1);
  lhs.fa = local_mesh->face_handle(lhs.a0);
  lhs.fb = local_mesh->face_handle(lhs.b0);

  // Flip the local mesh
  EdgeHandle local_flip_eh = local_mesh->edge_handle(lhs.a0);
  local_mesh->flip(local_flip_eh);
}

void LEdgeFlipper::update_global_mesh()
{
  // update global mesh
  // do real flip operation.
  mesh->flip(flip_eh);
  // update pixels, color and error in global mesh.
  mesh->data(ghs.fa).pixels = local_mesh->data(lhs.fa).pixels;
  mesh->data(ghs.fb).pixels = local_mesh->data(lhs.fb).pixels;
  mesh->data(ghs.fa).color = local_mesh->data(lhs.fa).color;
  mesh->data(ghs.fb).color = local_mesh->data(lhs.fb).color;
  mesh->data(ghs.fa).error = local_mesh->data(lhs.fa).error;
  mesh->data(ghs.fb).error = local_mesh->data(lhs.fb).error;
}

}// namespace ImageTriSimp
}// namespace GCLF