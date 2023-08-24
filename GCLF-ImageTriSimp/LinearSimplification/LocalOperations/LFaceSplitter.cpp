#include "LFaceSplitter.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

LFaceSplitter::LFaceSplitter(LinearSimplifier* _simplifier)
  :LLocalOperation(_simplifier)
{}

/// @brief Initialize neccessary data before flipping.
/// @return true if flip_ok.
bool LFaceSplitter::init(FaceHandle _f, size_t Np)
{
  clear();
  split_f = _f;

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

void LFaceSplitter::clear()
{
  LLocalOperation::clear();
  initialized = false;
}

bool LFaceSplitter::split_would_cause_over_valence(uint32_t max_valence)
{
  ASSERT(initialized, "splitter not initialized.");

  auto [v0, v1, v2] = face_vertices(*mesh, split_f);
  return mesh->valence(v0) < max_valence && mesh->valence(v1) < max_valence && mesh->valence(v2) < max_valence;
}

/// @brief Find affected faces on mesh and affected pixels on image.
void LFaceSplitter::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();

  global_affected_faces.insert(split_f);
  affected_pixel_coords.insert(affected_pixel_coords.end(),
    mesh->data(split_f).pixels.begin(), mesh->data(split_f).pixels.end());
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void LFaceSplitter::initialize_local_mesh()
{
  local_mesh = std::make_unique<LMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Split the local mesh
  local_split_f = FaceHandle(0);
  local_split_v = local_mesh->split(local_split_f, local_mesh->calc_centroid(local_split_f));
  local_optimizing_vertex = local_split_v;
}

void LFaceSplitter::update_global_mesh(const Vec3d& new_point)
{
  // update global mesh
  // do real split operation.
  split_v = mesh->split(split_f, new_point);
  // update pixels, color and error in global mesh.
  for (HalfedgeHandle voh : local_mesh->voh_range(local_split_v))
  {
    VertexHandle local_to_vh = local_mesh->to_vertex_handle(voh);
    FaceHandle local_fh = local_mesh->face_handle(voh);
    FaceHandle global_fh = mesh->face_handle(mesh->find_halfedge(split_v, rv2v[local_to_vh.idx()]));
    mesh->data(global_fh).pixels = local_mesh->data(local_fh).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fh).color;
    mesh->data(global_fh).error = local_mesh->data(local_fh).error;
  }
}

}// namespace ImageTriSimp
}// namespace GCLF