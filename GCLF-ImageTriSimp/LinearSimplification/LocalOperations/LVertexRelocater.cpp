#include "LVertexRelocater.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

LVertexRelocater::LVertexRelocater(LinearSimplifier* _simplifier)
  :LLocalOperation(_simplifier)
{}

/// @brief Initialize neccessary data before relocating.
/// @param _v vertex to relocate.
/// @return currently always true.
bool LVertexRelocater::init(VertexHandle _v, size_t Np)
{
  clear();
  relocate_vh = _v;
  if (mesh->is_boundary(relocate_vh) && is_corner_vertex(relocate_vh))
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

void LVertexRelocater::clear()
{
  LLocalOperation::clear();
  relocate_vh.invalidate();
  halfedges.clear();
}

/// @brief Find affected faces on mesh and affected pixels on image.
void LVertexRelocater::find_affected_faces_pixels()
{
  halfedges.reserve(mesh->valence(relocate_vh));
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  for (HalfedgeHandle voh : mesh->voh_range(relocate_vh))
  {
    FaceHandle fh = mesh->face_handle(voh);
    if (mesh->is_boundary(voh))
      continue;
    halfedges.push_back(mesh->next_halfedge_handle(voh));
    global_affected_faces.insert(fh);
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void LVertexRelocater::initialize_local_mesh()
{
  local_mesh = std::make_unique<LMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Find control points to optimize.
  local_optimizing_vertex = v2v.at(relocate_vh);

  do_project_x = false;
  do_project_y = false;
  if (mesh->is_boundary(relocate_vh))
  {
    Vec3d& p = mesh->point(relocate_vh);
    double x_diff = std::min(fabs(p.x() - 0.0), fabs(p.x() - image->width));
    double y_diff = std::min(fabs(p.y() - 0.0), fabs(p.y() - image->height));
    if (x_diff < y_diff)
    {
      do_project_x = true;
      x_project_to = p.x();
    }
    else
    {
      do_project_y = true;
      y_project_to = p.y();
    }
  }
}

void LVertexRelocater::update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point)
{
  new_local_mesh.set_point(local_optimizing_vertex, new_point);
  if (do_project_x)
    new_local_mesh.point(local_optimizing_vertex).x() = x_project_to;
  else if (do_project_y)
    new_local_mesh.point(local_optimizing_vertex).y() = y_project_to;
}

void LVertexRelocater::update_global_mesh(const Vec3d& new_point)
{
  // update global mesh
  mesh->set_point(relocate_vh, local_mesh->point(local_optimizing_vertex));
  for (FaceHandle local_fh : local_mesh->vf_range(local_optimizing_vertex))
  {
    FaceHandle global_fh = rf2f[local_fh.idx()];
    // update pixels, color and error
    mesh->data(global_fh).pixels = local_mesh->data(local_fh).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fh).color;
    mesh->data(global_fh).error = local_mesh->data(local_fh).error;
  }
}

void LVertexRelocater::post_update()
{
#ifdef ONLY_RELOCATE_HIGH_ERROR_VERTEX
  // update error defines on vertex
  mesh->data(relocate_vh).error = 0.0;
  for (FaceHandle vf : mesh->vf_range(relocate_vh))
    mesh->data(relocate_vh).error += mesh->data(vf).error;

  for (VertexHandle vv : mesh->vv_range(relocate_vh))
  {
    mesh->data(vv).error = 0.0;
    for (FaceHandle vf : mesh->vf_range(vv))
      mesh->data(vv).error += mesh->data(vf).error;
  }
#endif
}

}// namespace ImageTriSimp
}// namespace GCLF